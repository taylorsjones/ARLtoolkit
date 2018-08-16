#R version of the ARL reader
#require('proj4')   #needed for CRS definitions
#require('raster')  #creates raster objects
#require('rts')     #creates raster time series (RTS) objects
#require('data.table')


ARL_read_main_header <- function(f){
	h <- list()
	fh <- readChar(f,50)
	h$met_type <- readChar(f,4)
	h$forecast_hr <- readChar(f,3)
	h$forecast_mn <- readChar(f,2)
	h$pole_lat    <- as.numeric( readChar(f,7) )
	h$pole_lon    <- as.numeric( readChar(f,7) )
	h$ref_lat     <- as.numeric( readChar(f,7) )
	h$ref_lon     <- as.numeric( readChar(f,7) )
	h$dx          <- as.numeric( readChar(f,7) )
	h$orientation <- readChar(f,7)
	h$cone_angle  <- as.numeric( readChar(f,7) )
	h$x_sync      <- as.numeric( readChar(f,7) )
	h$y_sync      <- as.numeric( readChar(f,7) )
	h$sync_lat    <- as.numeric( readChar(f,7) )
	h$sync_lon    <- as.numeric( readChar(f,7) )
	h$reserved    <- readChar(f,7)
	h$nx          <- as.numeric( readChar(f,3) )
	h$ny          <- as.numeric( readChar(f,3) )
	if( h$met_type == "HRRR"){
		h$nx <- h$nx + 1000
		h$ny <- h$ny + 1000
	}
	h$nz          <- as.numeric( readChar(f,3) )
	h$v_coord     <- readChar(f,2)
	h$index_length<- as.numeric( readChar(f,4) )
	d2_vars <- NULL
	z <- NULL
	for( j  in 1:h$nz){
		d3_vars <- NULL
		z <- c(z,readChar(f,6))
		num_of_vars <- as.numeric( readChar(f,2) )
		for( i in 1:num_of_vars){
			if( j == 1){
				d2_vars <- c(d2_vars, readChar(f,4) )
			}else{
				d3_vars <- c(d3_vars, readChar(f,4) )
			}
			readChar(f,3) #check sum...ignore it.
			readChar(f,1) #reserved
		}
	}
	h$z <- z
	h$d2_vars <- d2_vars
	h$d3_vars <- d3_vars

	#build the projection:
	lat_1 = acos( h$cone_angle /360.0 )*(180/pi)
	h$proj <- paste0("+proj=lcc +lat_1=",h$cone_angle," +lat_2=",h$cone_angle," +lat_0=",h$ref_lat," +lon_0=",h$ref_lon,' +units=m +e=0 +a=6371229')
	h$xx <- seq(-1*((h$nx*h$dx*1000)-1)/2,((h$nx*h$dx*1000)-1)/2,h$dx*1000) - 520.143 #offset for HRRR
	h$yy <- seq(-1*((h$ny*h$dx*1000)-1)/2,((h$ny*h$dx*1000)-1)/2,h$dx*1000) - 306.153 #offset for HRRR

	return(h)
}

ARL_read_variable <- function(f,h,x1,x2,y1,y2){
	output <- list()
	while( readChar(f,1) == ''){}
	seek(f,-1,"current")
	datetime <- readChar(f,10)
	dont_know   <- readChar(f,4)
	var_name    <- readChar(f,4)
	scale       <- as.numeric( readChar(f,4) )
	precision   <- as.numeric( readChar(f,14) )
	last_value  <- as.numeric( readChar(f,14) )
	var_map     <- matrix(0,1+y2-y1,1+x2-x1)
	for( i in 1:h$ny){
		if(i < y1){ #skip rows early rows
			last_value <- (readBin(f,"int",n=1,size=1,signed=F)-127.0)/(2^(7-scale)) + last_value
			seek(f,h$nx-1,"current")
		}else if( i > y2){ #leave after the rows we care about
			seek(f,h$nx*(1+h$ny-i),"current")
			break
		}else{
			first_col_last_value <- last_value
			for( j in 1:h$nx){
				value <- (readBin(f,"int",n=1,size=1,signed=F)-127.0)/(2^(7-scale)) + last_value
				if(j > x1){
        	var_map[i-y1,j-x1] <- value
				}
		    if(j > x2){
		    	seek(f,h$nx-j,"current")
		    	break
		    }
		    last_value <- value
			}
		last_value <- first_col_last_value
		}
	}
	var_map <- apply(var_map,2,rev)
	output$name <- var_name
	output$map <- var_map
	return(output)
}

#' Read an ARL file
#'
#' @param file_to_read the name of the file to read in
#' @param var_i_want the name of the variable to extract
#' @param ll the lat and lon of the lower-left corner of the area to extract
#' @param ur the lat and lon of the upper-right corner of the area to extract
#' @param latlon whether to convert the raster to lat-lon grid. Otherwise it stays in the native projection
#' @param verbose whether to output additional messages during extraction'
#' @return for 2d variables, a raster time series object is returned. For 3d variables, a list of RTS objects is returned, one for each altitude.
#' @examples
#' file_to_read <- 'hysplit.20161020.00z.hrrra'
#' lower_left <- c(-73,40)
#' upper_left <- c(-70,43)
#' pbl <- ARL_read(file_to_read,"PBLH",ll=lower_left,ur=upper_right)
#' plot(pbl)
#' @export
ARL_read <- function(file_to_read,var_i_want,ll,ur,latlon=TRUE, verbose=FALSE){
	f = file(file_to_read,'rb')
	raster_list <- list()
	levels <- NULL
	datetimes <- NULL
	for(step in seq(1,2000,1)){
		date <- readChar(f,10)
		level <- readChar(f,2)
		AA <- readChar(f,2)
		block_ID <- readChar(f,4)
		seek(f,-18,"current")
		if( length(block_ID) < 1 ){
			if(verbose){ message("...end of file") }
			break
		}
		if(block_ID == "INDX"){ #header detected...load in all the header info.
			h <- ARL_read_main_header(f)
			seek(f,h$nx*h$ny - h$index_length , "current")
			ll_xy <- project( ll, h$proj )
			x1 <- findInterval(ll_xy[1],h$xx)
			y1 <- findInterval(ll_xy[2],h$yy)
			ur_xy <- project( ur, h$proj )
			x2 <- findInterval(ur_xy[1],h$xx)
			y2 <- findInterval(ur_xy[2],h$yy)
			if(verbose){
				message("Array Dimensions: ")
				message(x2-x1)
				message(y2-y1)
			}
			if(var_i_want %in% h$d2_vars){
				if(verbose){message("2-d var detected")}
			}else if(var_i_want %in% h$d3_vars){
				if(verbose){message("3-d var detected")}
			}else{
				message("var not found! These vars are available:")
				message("2-dimensional:")
				print(h$d2_vars)
				message("3-dimensional:")
				print(h$d3_vars)
				break
			}

		}else if(block_ID == var_i_want){
			out <- ARL_read_variable(f,h,x1,x2,y1,y2)
			datetime <- paste0("20",gsub(" ","0",date))
			ras <- raster(out$map ,crs=h$proj, xmn=ll_xy[1], xmx=ur_xy[1], ymn=ll_xy[2], ymx=ur_xy[2] )
			if(latlon){
				ras <- projectRaster(ras,crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
			}
			raster_list <- c(raster_list,ras)
			datetimes <- c(datetimes,datetime)
			levels <- c(levels,level)
		}else{
			seek(f,50 + h$nx*h$ny,"current")
		}
	}

	close(f)

	unique_levels <- unique(levels)
	unique_datetimes <- unique(datetimes)

	if( length( unique_levels ) > 1){
		rts_list <- list()
  	output <- list()

		for( level in unique_levels){
			sub_list_of_rasters <- raster_list[ level == levels ]
			sub_list_of_datetimes <- datetimes[ level == levels ]
			rstack <- stack(sub_list_of_rasters)
			rts_object <- rts(rstack, as.POSIXct(sub_list_of_datetimes,format="%Y%m%d%H%M",tz="GMT"))
			rts_list <- c(rts_list,rts_object)
		}

		output$levels <- h$z
		output$rasters <- rts_list
		return(output)
	}else{
		rstack <- stack(raster_list)
		rts_object <- rts(rstack, as.POSIXct(datetimes,format="%Y%m%d%H%M",tz="GMT"))
		return(rts_object)
	}
}





#' Read an FSL formatted radiosonde file
#'
#' @param file_name name of the file to read
#' @param xlim left and right longitude limits to extract
#' @param ylim lower and upper latitude limits to extract
#' @return a data table is returned with launch timestamp and location, altitude, pressure, temperature, and U and V wind components
#' @examples
#' file_name <- 'roab_soundings43993.txt'
#' DT <- FSL_read(file_name,xlim=c(-125,-65),ylim=c(20,50))
#' sub_DT <- subset(DT, tstamp == DT$tstamp[1] )
#' plot(sub$lon, sub$lat)
#' map("worldHires",add=True)
#' @export
FSL_read <- function(file_name,xlim=c(-180,180),ylim=c(-90,90)){

  ts_strings <- NULL
  lats       <- NULL
  lons       <- NULL
  agls       <- NULL
  presses    <- NULL
  temps      <- NULL
  Us         <- NULL
  Vs         <- NULL
  SIDs       <- NULL

  con  <- file(file_name,open="r")

  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    #line <- readLines(con, n = 1, warn = FALSE)
    line_type <- substr(line,5,7)

    if(line_type == "254"){  #new sounding
      hour  <- substr(line,13,14)
      day   <- substr(line,20,21)
      month <- substr(line,28,30)
      year  <- substr(line,35,38)
      ts_string <- paste0(hour,day,month,year)
      ts_string <- gsub(" ","0",dt_string)
    }

    if(line_type == "  1"){  #station ID line
      WBAN <- substr(line,10,14)
      WMO  <- substr(line,17,21)
      lat  <- as.numeric(substr(line,22,28))
      lat_NS <- substr(line,29,29)
      if(lat_NS == "S"){
        lat <- -1*lat
      }
      lon  <- as.numeric(substr(line,30,35))
      lon_EW <- substr(line,36,36)
      if(lon_EW == "W"){
        lon <- -1*lon
      }
      elev  <- substr(line,37,42)
      rtime <- substr(line,43,50)
    }

    if(line_type == "  2"){                   #Sounding check lines
      n_of_l <- strtoi(substr(line,29,35))    #number of lines in this entry, including header(s)
      if( (lon>xlim[2]) | (lon<xlim[1]) | (lat>ylim[2]) | (lat<ylim[1]) ){
        dump <- readLines(con, n = n_of_l - 3) #skip entries outside our domain of interest
      }
    }

    if(line_type == "  3"){ #Station ID and other indicators
      SID      <- substr(line,18,21)     # Station name
      stype    <- substr(line,38,42)     # Type of sonde used
      ws_units <- substr(line,43,50)     # wind speed units (ms means tenths of m/s)
    }

    if( (line_type == "  4") | (line_type == "  5") | (line_type == "  7") | (line_type == "  9") ){ #data level
      dew_pt     <- substr(line,29,35)          # dewpoint in tenths of deg C
      wind_dir   <- strtoi(substr(line,36,42))  # wind direction in degrees
      wind_spd   <- strtoi(substr(line,43,49))  # wind speed check ws_units for units
      ts_strings <- c(ts_strings,ts_string)
      SIDs       <- c(SIDs,SID)
      lats       <- c(lats,lat)
      lons       <- c(lons,lon)
      agls       <- c(agls,strtoi(substr(line,15,21)) )
      presses    <- c(presses, strtoi(substr(line,8,14)) )
      temps      <- c(temps, strtoi(substr(line,22,28)) )
      Us         <- c(Us, cos(pi*wind_spd/180) )
      Vs         <- c(Vs, sin(pi*wind_spd/180) )
    }
  }

  close(con)
  tstamps <- as.POSIXct(ts_strings,tz="GMT",format="%H%d%b%Y")
  DT <- data.table(tstamp=tstamps,station=SIDs,lat=lats,lon=lons,agl=agls,press=presses,temp=temps,u=Us,v=Vs)
  return(DT)

}

semilog_func <- function(a,d,L){
  return( a*(1-exp(-1*d/L) ) )
}
#foo <- ARL_read(file,"V10M",ll=lower_left,ur=upper_right,latlon=FALSE)
#woot <- Variogram(foo[[1]])
#x <- woot@variogram$distance
#y <- woot@variogram$gamma
#zfit <- nls( y~semilog_func(a,x,L),start=list(a=1,L=3e4))


library("sp")
library("rgdal")
library("pbapply")

## Converts character strings to Title Case
titleCase <- function(x) { s <- strsplit(x, " ")[[1]]; paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ") }

## Given lat and long, and an sp::bbox object, return a matrix representing the points within the bounding box
subset_by_bbox <- function(lat, long, bbox, by_index = FALSE){
  long_bool <- long >= bbox[1, 1] & long <= bbox[1, 2]
  lat_bool <- lat >= bbox[2, 1] & lat <= bbox[2, 2]
  if (by_index){ return(which(long_bool & lat_bool)); }
  else { return(cbind(long = long[long_bool & lat_bool], lat = lat[long_bool & lat_bool])) }
}

## Given geocoded results, return address/message information
extractMessage <- function(res){
  if (res$status == "OK"){
    place <- res$results[[1]]
    place$types <- paste(strsplit(place$types, split = "_")[[1]], collapse = " ")
    place$geometry$location_type <- paste(strsplit(place$geometry$location_type, split = "_")[[1]], collapse = " ")
    sprintf("Address: %s \nPlace Type(s): %s \nLocation Method: %s", place$formatted_address, titleCase(place$types), place$geometry$location_type) 
  } else { "Reverse geocode lookup failed." }
}

## Extract Points of interest ids associated with points. Includes both the centroid and the original points that made up the POI.
extractPOIs <- function(x, cl){
  cl_ids <- unique(cl)
  cl_ids <- cl_ids[cl_ids != 0] ## exclude noise
  res <- pblapply(cl_ids, function(id) { 
    tmp <- x[which(cl == id),] 
    list(centroid=sp::Polygon(tmp)@labpt, ## sp::Polygon computes the center of mass, more representative of 'center' than mean centroid
         exemplars=tmp) 
  })
  names(res) <- as.character(cl_ids)
  res
}

## Extract clusters ids associated with points.
extractCentroids <- function(x, cl){
  cl_ids <- unique(cl)
  cl_ids <- cl_ids[cl_ids != 0] ## exclude noise
  cl_centroids <- pblapply(cl_ids, function(id) { 
    tmp <- x[which(cl == id),] 
    sp::Polygon(tmp)@labpt ## sp::Polygon computes the center of mass, more representative of 'center' than mean centroid
  })
  do.call(rbind, cl_centroids)
}

# Adapted from http://www.igorexchange.com/node/927
getUTMZone <- function(Long, Lat){
  LongTemp <- (Long+180)-floor((Long+180)/360)*360-180
  ZoneNumber <- floor((LongTemp + 180)/6) + 1
  if (Lat >= 56.0 && Lat < 64.0 && LongTemp >= 3.0 && LongTemp < 12.0 )
    ZoneNumber <- 32
  if( Lat >= 72.0 && Lat < 84.0 ){
    if  ( LongTemp >= 0.0  && LongTemp <  9.0 )
      ZoneNumber = 31
    else if( LongTemp >= 9.0  && LongTemp < 21.0 )
      ZoneNumber = 33
    else if(LongTemp >= 21.0 && LongTemp < 33.0 )
      ZoneNumber = 35
    else if(LongTemp >= 33.0 && LongTemp < 42.0 )
      ZoneNumber = 37
  }
  return(as.integer(ZoneNumber))
}

## Expects a matrix-coercible object of the form < lat, long >
## Zone specifies the UTM zone.
wgs84_to_meters <- function(lat, long){
  zone <- getUTMZone(Long = long, Lat = lat)
  if (length(unique(zone)) > 1){ stop("'wgs84_to_meters' expects the WGS84 coordinates passed to all be part of the same UTM zone.")} 
  zone <- unique(zone)[[1]]
  x <- cbind(long, lat)
  geo_pts <- sp::SpatialPoints(as.matrix(x), proj4string = CRS("+proj=longlat +ellps=WGS84"))
  metric_pts <- sp::spTransform(geo_pts, CRS(paste0("+proj=utm +zone=", zone, " +ellps=WGS84 +datum=WGS84 +units=m +no_defs")))
  coords <- metric_pts@coords
  new_coords <- cbind(coords[, 1] - min(coords[, 1]), coords[, 2] - min(coords[, 2]))
  offsets <- c(min(coords[, 1]), min(coords[, 2]))
  return(list(coords = new_coords, offset = offsets, utm_zone = zone))
}

## convert metric back to lat/long w/ the appropriate offset
meters_to_wgs84 <- function(metric_coords){
  if (!all(c("coords", "offset", "utm_zone") %in% names(metric_coords))){ stop("'meters_to_wgs84' expects a list of coordinates as returned by 'wgs84_to_meters'")}
  tmp_pts <- metric_coords$coords
  tmp_pts[, 1] <- tmp_pts[, 1] + metric_coords$offset[1]
  tmp_pts[, 2] <- tmp_pts[, 2] + metric_coords$offset[2]
  metric_pts <- sp::SpatialPoints(tmp_pts, proj4string = CRS(paste0("+proj=utm +zone=", metric_coords$utm_zone, " +ellps=WGS84 +datum=WGS84 +units=m +no_defs")))
  geo_pts <- sp::spTransform(metric_pts, CRS("+proj=longlat +ellps=WGS84"))
  colnames(geo_pts@coords) <- c("longitude", "latitude")
  return(geo_pts@coords)
}

## Convert two lat/long pairs to distance in meters
toMeters <- function(lat1, lon1, lat2, lon2){
  R = 6378.137 ## Radius of earth in KM
  dLat = lat2 * pi / 180 - lat1 * pi / 180
  dLon = lon2 * pi / 180 - lon1 * pi / 180
  a = sin(dLat/2) * sin(dLat/2) +
    cos(lat1 * pi / 180) * cos(lat2 * pi / 180) *
    sin(dLon/2) * sin(dLon/2)
  c = 2 * atan2(sqrt(a), sqrt(1-a))
  d = R * c
  return(d * 1000) #  meters
}


## Compile and load the stay point C++ function
if (length(find("stay_pts")) == 0){ Rcpp::sourceCpp("~/WaCS/Work/stay_pts.cpp", embeddedR = FALSE) }

## Parse through trajectories of <lat, lng> pairs, extracting stay points that satisfy the 
## given parameter settings 'd_thresh' and 't_thresh'. That is, any trajectory segment containing consecutive 
## coordinates each of which 1) are less than d_thresh apart (meters) and 2) collectively occur over a period of at least 
## 't_thresh' amount of time (units defined by the timesteps given), then the medoid is extracted as a 'stay point'
extractStayPoints <- function(lat, lng, timesteps, d_thresh, t_thresh){
  if (any(is.na(lat)) || any(is.na(lng)) || length(lat) < 3L){ return(NULL) }
  metric_coords <- try(wgs84_to_meters(lat = lat, long = lng), silent = TRUE)
  if (class(metric_coords) == "try-error") return(NULL)
  sps <- stay_pts(metric_coords$coords, # Metric coordinates reflect euclidean distance 1:1 w/ meters 
                  timesteps = timesteps, # timestep units determined by what's passed in
                  d_thresh = d_thresh, # within d_thresh meters 
                  t_thresh = t_thresh) # for at least d_thresh time
  if (nrow(sps) == 0){ return(NULL) }
  metric_coords$coords <- matrix(sps[, 1:2], ncol = 2)
  sps_lonlat <- try(meters_to_wgs84(metric_coords), silent = TRUE)
  if (class(sps_lonlat) == "try-error") return(NULL)
  return(sps_lonlat)
}

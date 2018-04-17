library("clv")
library("clues")
library("cluster")
library("dbscan")

## Davies Bouldin score (and Dunn Index; more computationally efficient when together)
DB <- function(x, cl, ...){
  if (all(cl == cl[1])) { return(-Inf) }
  noise_idx <- cl == 0
  valid_cl <- cl[!noise_idx]
  valid_x <- x[!noise_idx,]
  n <- length(valid_cl)
  if (n > 0){
    # cl_scatt_data <- clv::cls.scatt.data(data = valid_x, clust = valid_cl, dist = "euclidean")
    # db_index <- clv::clv.Davies.Bouldin(index.list = cl_scatt_data, intracls = c("complete"), intercls = c("single"))
    # dunn_index <- clv::clv.Dunn(index.list = cl_scatt_data, intracls = c("complete"), intercls = c("single"))
    db_index <- clusterSim::index.DB(x = valid_x, cl = valid_cl)
    return(list(DB=db_index, ...))
  } else {
    return(list(DB=NULL, ...))
  }
}

## Silhoutte global score
SIL <- function(x, cl, xdist, ...){
  if (all(cl == cl[1])) { return(-Inf) }
  noise_idx <- cl == 0
  valid_cl <- cl[!noise_idx]
  valid_x <- x[!noise_idx,]
  n <- length(valid_cl)
  sil_score <- cluster::silhouette(x = cl, dist = xdist)
  if (any(is.na(sil_score))){
    return(list(SIL=NA, ...))
  } else {
    sil_avgs <- summary(sil_score)$clus.avg.widths
    sil_global_index <- mean(sil_avgs[which(names(sil_avgs) != "0")])
    return(list(SIL=sil_global_index, ...))
  }
}

## Calinski-Harabasz score
CH <- function(x, cl, ...){
  if (all(cl == cl[1])) { return(-Inf) }
  noise_idx <- cl == 0
  valid_cl <- cl[!noise_idx]
  valid_x <- x[!noise_idx,]
  n <- length(valid_cl)
  ch_index <- fpc::calinhara(x = valid_x, clustering = valid_cl, cn = length(unique(valid_cl)))
  # ch_index <- clues::get_CH(y = valid_x, mem = valid_cl, disMethod = "Euclidean")
  list(CH=ch_index, ...)
}

## Density Based Clustering Validation
DBCV <- function(x, cl, xdist, ...){
  if (all(cl == cl[1])) { return(-Inf) }
  dbcv_index <- dbscan::dbcv(x = x, cl = cl, xdist = xdist)
  list(DBCV=dbcv_index, ...)
}


## Dunn Index
DUNN <- function(x, cl, ...){
  if (all(cl == cl[1])) { return(NULL) }
  noise_idx <- cl == 0
  valid_cl <- cl[!noise_idx]
  valid_x <- x[!noise_idx,]
  n <- length(valid_cl)
  cl_scatt_data <- clv::cls.scatt.data(data = valid_x, clust = valid_cl, dist = "euclidean")
  dunn_index <- clv::clv.Dunn(index.list = cl_scatt_data, intracls = c("complete"), intercls = c("single"))
  list(DUNN=dunn_index, ...)
}
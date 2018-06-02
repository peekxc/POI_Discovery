## run_clustering.R
## Simple function to run several clustering algorithms across a swath of parameter settings. 
## Each clustering algorithm is run within it's own local environment, with the results saved to 
## disk after finishing. This is done to minimize memory usage between algorithms.  
## Author: Matt Piekenbrock

## all_clusters
## ds_mat := Matrix of numeric observations in R^d
## folder_name := relative path of folder to save clustering results 
## core_num := number of cores to parallelize with
all_clusters <- function(ds_mat, folder_name, core_num, which_cl = c("average","kmeans", "single", "mcquitty", "optics", "ward", "dbscan", "rsl")){
  if (!dir.exists(folder_name)){ dir.create(folder_name) }
  library("pbapply")
  library("dbscan")
  library("clustertree")
  library("parallelDist")
  
  ## Setup parameter settings 
  k <- seq(5, min(c(300, nrow(ds_mat) - 1)), by = 5)
  ds_dist <- parallelDist::parallelDist(ds_mat, method = "euclidean")
  minPts <- seq(2, 40, by = 5)
  
  if ("kmeans" %in% which_cl){
    cat("Running K-means\n")
    with(new.env(), expr = {
      km_cl <- pbsapply(k, function(ki) { kmeans(x = ds_mat, centers = ki)$cluster }, cl = core_num)
      save(km_cl, file = file.path(paste0(folder_name, "/km_cl.rdata")))
    }); gc()
  }
  
  if ("single" %in% which_cl){
    cat("Running Single Linkage\n")
    with(new.env(), expr = {
      sl <- hclust(ds_dist, method = "single")
      eps <- quantile(sl$height, probs = seq(0.05, 0.95, by = .05))
      sl_cl <- pbsapply(eps, function(eps_h) { 
        cutree(sl, h = eps_h)
      }, cl = core_num)
      save(sl_cl, file = file.path(paste0(folder_name, "/sl_cl.rdata"))) 
    }); gc()
  }
  
  if ("average" %in% which_cl){
    cat("Running Average Linkage\n")
    with(new.env(), expr = {
      al <- hclust(ds_dist, method = "average")
      eps <- quantile(al$height, probs = seq(0.05, 0.95, by = .05))
      al_cl <- pbsapply(eps, function(eps_h) { 
        cutree(al, h = eps_h)
      }, cl = core_num)
      save(al_cl, file = file.path(paste0(folder_name, "/al_cl.rdata"))) 
    }); gc()
  }
  
  if ("mcquitty" %in% which_cl){
    cat("Running McQuitty Linkage\n")
    with(new.env(), expr = {
      mcq <- hclust(ds_dist, method = "mcquitty")
      eps <- quantile(mcq$height, probs = seq(0.05, 0.95, by = .05))
      mcq_cl <- pbsapply(eps, function(eps_h) {
        cutree(mcq, h = eps_h)
      }, cl = core_num)
      save(mcq_cl, file = file.path(paste0(folder_name, "/mcq_cl.rdata"))) 
    }); gc()
  }
  
  if ("ward" %in% which_cl){
    cat("Running Ward Linkage\n")
    with(new.env(), expr = {
      wrd <- hclust(ds_dist, method = "ward.D2")
      eps <- quantile(wrd$height, probs = seq(0.05, 0.95, by = .05))
      wrd_cl <- pbsapply(eps, function(eps_h) {
        cutree(wrd, h = eps_h)
      }, cl = core_num)
      save(wrd_cl, file = file.path(folder_name, "wrd_cl.rdata"))
    }); gc()
  }
  
  ## DBSCAN Tests
  if ("dbscan" %in% which_cl){
    cat("Running DBSCAN\n")
    with(new.env(), expr = {
      db_cl <- pblapply(minPts, function(mpts){
        ## HDBSCAN produces a hiearchy where cuts are equivalent to DBSCAN* clusters
        hdb <- dbscan::hdbscan(x = ds_mat, minPts = mpts, xdist = ds_dist)
        eps <- quantile(hdb$hc$height, probs = seq(0.05, 0.95, by = .05))
        do.call(cbind, lapply(eps, function(ep){ cutree(hdb$hc, h = ep) }))
      }, cl = core_num)
      db_cl <- do.call(cbind, db_cl)
      save(db_cl, file = file.path(paste0(folder_name, "/dbscan_cl.rdata"))) 
    }); gc()
  }
    
  ## OPTICS
  if ("optics" %in% which_cl){
    cat("Running OPTICS\n")
    with(new.env(), expr = {
      op_cl <- pblapply(minPts, function(mpts){
        op <- dbscan::optics(x = ds_mat, eps = Inf, minPts = mpts)
        op_eps <- quantile(op$reachdist, probs = seq(0.05, 0.95, by = .05))
        do.call(cbind, lapply(op_eps, function(ep){
          dbscan::extractDBSCAN(object = op, eps_cl = ep)$cluster
        }))
      }, cl = core_num)
      op_cl <- do.call(cbind, op_cl)
      save(op_cl, file = file.path(paste0(folder_name, "/optics_cl.rdata"))) 
    }); gc()
  }
  
  
  ## RSL / Cluster Tree 
  if ("rsl" %in% which_cl){
    cat("Running RSL\n")
    with(new.env(), expr = {
      k_init <- floor(ncol(ds_mat)*log(nrow(ds_mat))) ## Optimal setting as given by Chaudhuri's analysis. 
      k <- seq(k_init, 2*k_init, by = 2) ## Practically speaking, noisier data sets may require a slightly higher setting
      rsl_cl <- pblapply(k, function(k_){
        cltree <- clustertree::clustertree(x = ds_dist, k = k_, alpha = sqrt(2), estimator = "RSL")
        minpts <- seq(5, k_, by = 5L)
        do.call(cbind, lapply(minpts, function(min_size){ dbscan::extractFOSC(cltree$hc, minPts = min_size)$cluster }))
      })
      rsl_cl <- do.call(cbind, rsl_cl)
      save(rsl_cl, file = file.path(paste0(folder_name, "/rsl_cl.rdata"))) 
    }); gc()
  }
}
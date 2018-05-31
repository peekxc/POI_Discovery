## validate.R
## Simple function to load run multiple validation scores for each 
## Author: Matt Piekenbrock

validate <- function(ds_mat, folder_name, core_num, 
                     which_cl = c("al_cl","dbscan_cl", "km_cl", "mcq_cl", "optics_cl", "sl_cl", "wrd_cl", "rsl_cl"), 
                     root_dir = "~/ClusterTree_POIs"){
  if (!dir.exists(folder_name)){ stop(paste0("Folder path ", folder_name, " not found.")) }
  if (is.null(dim(ds_mat))){ stop("Stopping. 'ds_mat' must be a matrix-coercible object.") }
  
  ## Check all files exist
  cluster_path <- file.path(root_dir, folder_name)
  cls_to_validate <- paste0(which_cl, ".rdata")
  if (any(!cls_to_validate %in% list.files(cluster_path))){
    browser()
    stop("Not all cluster configurations found.")
  }
  
  ## Cretae new score directory if it doesn't exist
  score_dir <- file.path(root_dir, folder_name, "scores")
  if (!dir.exists(score_dir)){ dir.create(score_dir) }
  
  ## Compute distance matrix
  cat("Computing distances...\n")
  ds_dist <- parallelDist::parallelDist(x = ds_mat, method = "euclidean")
  
  ## Run validation
  for (cl_file in cls_to_validate){
    
    ## Load the data set as 'x_cl'; expected to be a matrix of clusters
    cl_env <- new.env(parent = .BaseNamespaceEnv)
    load(file.path(cluster_path, cl_file), envir = cl_env)
    cl_name <- ls(cl_env)
    x_cl <- cl_env[[cl_name]]
    if (is.null(dim(x_cl)) || is.list(x_cl)) { stop("'x_cl' not a matrix-like object.") }
    
    ## If the score file already exists, don't re-run validation
    score_file <- file.path(score_dir, paste0(cl_name, ".rdata"))
    if (file.exists(score_file) && missing(which_cl)){ print(sprintf("Skipping %s; score file already exists.", cl_name)); next; }
    
    ## Run validation scores
    cat("Running DBCV score...\n")
    x_dbcv <- pbapply::pbapply(x_cl, 2, function(cl){ DBCV(x = ds_mat, cl = cl, xdist = ds_dist) },  cl = core_num)
    
    # cat("Running Calinski-Harabasz score...\n")
    # x_ch <- pbapply::pbapply(x_cl, 2, function(cl){ CH(x = ds_mat, cl = cl) }, cl = core_num)
    
    # cat("Running Silhouette score...\n")
    # x_sil <- pbapply::pbapply(x_cl, 2, function(cl){ SIL(x = ds_mat, cl = cl, xdist = ds_dist) }, cl = core_num)
    
    # { rm(ds_dist); gc() } ## cleanup 
    
    # cat("Running Davies Bouldin and Dunn Index scores...\n")
    # x_db <- pbapply::pbapply(x_cl, 2, function(cl){ DB(x = ds_mat, cl = cl) },  cl = core_num)
    
    ## Save scores 
    cat(sprintf("Saving scores for %s...\n", cl_file))
    x_scores <- list(x_dbcv) # x_db, x_ch, x_sil, 
    save(x_scores, file = score_file)
  }
}
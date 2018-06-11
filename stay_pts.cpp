#include <Rcpp.h>
using namespace Rcpp;

NumericVector col_means( NumericMatrix& X ) {
  int nCols = X.ncol();
  NumericVector out = no_init(nCols);
  for( int j=0; j < nCols; j++ ) {
    NumericMatrix::Column tmp = X(_, j);
    out[j] = mean(tmp);
  }
  return out;

}

// Returns squared euclidean distance
// [[Rcpp::export]]
static inline double euc_dist(NumericVector x, NumericVector y){
  return sum((x - y) * (x - y));
}


// Expects: 
// x := Numeric Matrix of < x, y > \in R^2.
// timesteps := Numeric vector of the time steps associated with each point. Length should be equal to the number of rows of x. 
// d_thresh := distance threshold. Consecutive points must be below this threshold to be considered a segment. 
// t_thresh := time threshold. Consecutive points must be below this threshold to be considered a segment. 
// medoid := boolean. If TRUE, returns the point in x, for each segment, which has the smallest distance to the next point. If FALSE, the mean of the segment is used. 
// returnIndex := (Future Work)
// Description := From the paper: "In this experiment, we set t_thresh to 20 minutes and d_thresh to 200 meters for stay point detection.
// In other words, if an individual stays over 20 minutes within a distance of 200 meters, a stay point is detected.
// [[Rcpp::export]]
NumericMatrix stay_pts(NumericMatrix& x, const NumericVector& timesteps, double d_thresh, double t_thresh, 
                       bool medoid = true, bool returnIndex = false) {

  int n = x.nrow();
  //Rcout << "here" << std::endl;
  // Don't consider less than 3 points to make a stay point
  if (n < 3){
    NumericMatrix res = no_init_matrix(0, 3);
    return(res);
  }

  // Compute distance from each point to the next point
  NumericVector next_dist = Rcpp::no_init(n - 1);
  for (int i = 0; i < n - 2; ++i){
    next_dist.at(i) = euc_dist(x(i, _), x(i+1, _));
  }

  // Check to see which distances are greater than the given threshold
  LogicalVector dist_checks = next_dist < (d_thresh * d_thresh); // use squared distance

  // Get Run time length encoding
  Function rle("rle");
  List rle_res = rle(dist_checks);

  // Retrieve rle lengths, both regular and cumulative
  IntegerVector rle_len = rle_res["lengths"];
  IntegerVector cum_idx = cumsum(rle_len);

  // Iterate through and collect stay points
  int idx = 0;
  std::vector<double> sp_x = std::vector<double>();
  std::vector<double> sp_y = std::vector<double>();
  for (IntegerVector::iterator slen = rle_len.begin(); slen != rle_len.end(); ++slen, ++idx){
    int begin = idx == 0 ? 0 : cum_idx.at(idx - 1), end = cum_idx.at(idx) - 1; // inclusive start and end
    int v_len = (*slen);
    if (end == 0 || end <= begin) continue;
    if (begin < 0 || end >= next_dist.size() - 1) continue;
    if (v_len > 1 && dist_checks.at(begin) && timesteps.at(end) - timesteps.at(begin) >= t_thresh){
      //Rcout << "TRUE: begin: " << begin << ", end: " << end << std::endl;
      NumericMatrix subset = x(Rcpp::Range(begin, end), _); // -1 for 0-based inclusive
      // NumericVector stay_pt = col_means(subset); stay point defined by the mean of all points
      NumericVector stay_pt = NumericVector();
      if (medoid){ stay_pt = subset(which_min(next_dist[Range(begin, end)]), _); }
      else { stay_pt = col_means(subset); }
      if (stay_pt.size() == 2){
        sp_x.push_back(stay_pt.at(0));
        sp_y.push_back(stay_pt.at(1));
      }
    }
  }
  const int sp_n = sp_x.size();
  NumericMatrix stay_pts_res;
  if (returnIndex){
    stay_pts_res = no_init_matrix(sp_n, 3);
  } else {
    stay_pts_res = no_init_matrix(sp_n, 4);
  }
  for (int i = 0; i < sp_n; ++i){
    stay_pts_res(i, _) = returnIndex ? NumericVector::create(sp_x.at(i), sp_y.at(i), timesteps.at(i)) : NumericVector::create(sp_x.at(i), sp_y.at(i), timesteps.at(i), i);
  }
  return(stay_pts_res);
}

/*** R
# traj_data <- track_data[order(timestep), .(id, x, y, timestep)]
#
# ## If a car moves less than 10m each timestep for at least 2 seconds
# wut <- stay_pts(as.matrix(traj_data[id==0]), d_thresh = 5, t_thresh = 2)
#
#
# plot(traj_data[id == 0][, .(x, y)])
# points(as.matrix(data.table::rbindlist(lapply(wut, as.list))), col= "red")

# sp_list <- lapply(unique(sort(traj_data$id)), function(cid){
#   sp <- stay_pts(as.matrix(traj_data[id==cid]), d_thresh = 1000, t_thresh = 0.00001)
#   as.matrix(data.table::rbindlist(lapply(sp, as.list)))
# })
#
# traj_data[, .(sp = stay_pts(as.matrix(.SD), d_thresh = 1000, t_thresh = 0.00001)), by = id]
# traj_data[, .(stay_pts(as.matrix(.SD), d_thresh = 1000, t_thresh = 0.00001)), by = id]
# sum(sapply(sp_list, length))
#
# traj_data[id==0]
# plot(traj_data[, 2:3])
# points(, col = "red", cex=2)

Rcpp::sourceCpp('~/WaCS/net-analysis/stay_pts.cpp', embeddedR = FALSE)
  library("data.table")
  load(file = "veh_data.rdata")

   #, t1, 5, 1
  veh_data <- CLIF[, .(veh_id = trackID, V1 = lat, V2 = long, time = frame)]

  test <- CLIF[sample(1:nrow(CLIF), size = 1000), .(lat, long)]



  library(leaflet)
  poi_centroids_latlong <- sumo_net$projectCRS(poi_centroids)@coords
  m <- leaflet() %>%
  addTiles() %>%  # Add default OpenStreetMap map tiles
    addCircles(lng = poi_centroids_latlong[, 1], poi_centroids_latlong[, 2])
  # addMarkers(lng=test$long, lat=test$lat, popup="The birthplace of R")


  coords <- wgs84_to_meters(veh_data[, .(V1, V2)])
  veh_data$V1 <- coords[, 1] ## Replace latitude with metric space
  veh_data$V2 <- coords[, 2] ## Replace longitude with metric space
  sps <- extractSPs(veh_data, d_thresh = 5L, t_thresh = 1L)



  inp <- as.matrix(all_sps[, .(X1, X2)])
  # inp <- inp + runif(nrow(inp), max = 3) ## add up to 3 meters of random noise
  poi_model <- clustertree::clustertree(inp, k = 5L) ## Create a cluster tree model of the points
  # poi_cl <- dbscan::extractFOSC(poi_model$hc, minPts = 2L) ## Extract stability-based clusters
  poi_cl <- cutree(poi_model$hc, h = 25)

  table(poi_cl)
  length(unique(poi_cl))
  # plot(inp, col = poi_cl$cluster+1, asp = 1)

  cl_ids <- sort(unique(poi_cl))[-1]
  POIs <- lapply(cl_ids, function(poi_id){
    matrix(inp[poi_cl == poi_id,], ncol = 2)
  })
  POIs <- lapply(POIs, function(poi) poi[chull(poi),])

  clif_net <- createGeoNet(veh_data, POIs)

  ## Creates the dynamic network of the geospace
  ## Expects as arguments:
  ## 1) Trajectory data of the form < agent_id, x, y, timestep>
  ## 2) A list of point geometries representing locations to track
  createGeoNet <- function(tracks, POIs){
    if (ncol(tracks) != 4) { stop("Expecting a 4 column matrix-coercible object of track data of the form of the form < agent_id, x, y, time>") }
    if (!is(POIs, "list") || length(POIs) == 0) { stop("Expecting POIs to be a list of point geometries.") }
    POIs <- lapply(POIs, function(poi) matrix(poi, ncol =2))
    colnames(tracks) <- c("agent_id", "x", "y", "timestep")
    tracks <- data.table::data.table(tracks)
    track_ids <- sort(unique(tracks$agent_id))
    n <- length(POIs)
    poi_ids <- seq(n)
    offset <- min(tracks$timestep)-1
    tracks$timestep <- tracks$timestep - offset

    ## Get the record of visited locations for each 'agent'
    cat("\rTracing agent visitation records...            ")
    pois_visited <- vector(mode = "list", length = length(track_ids))
    pb <- txtProgressBar(max = length(track_ids), style = 3)
    for (i in seq(length(track_ids))){
      id <- track_ids[i]
      base_traj <- tracks[agent_id == id][order(timestep)]
      traj <- as.matrix(base_traj[, .(x, y)])
      visited <- lapply(poi_ids, function(poi_id){
        poi <- matrix(POIs[[poi_id]], ncol = 2)

        ## First check: if the point ever inside the POI, record the earliest time visited
        pip_res <- sp::point.in.polygon(point.x = traj[, 1], point.y = traj[, 2], pol.x = poi[, 1],  pol.y = poi[, 2])
        if(any(pip_res > 0)){ return(data.frame(agent_id = id, time_index = base_traj$timestep[which(pip_res > 0)[1]], poi = poi_id)) }

        ## Second check: if the trajectory ever gets within m meters of the centroid of the POI, record the earliest time visited
        centroid <- apply(poi, 2, mean)
        distances <- apply(traj, 1, function(pt) euc_dist(pt, centroid))
        if (any(distances < 25)){ return(data.frame(agent_id = id, time_index = base_traj$timestep[which(distances < 25)[1]], poi = poi_id)) }
        return(NULL)
      })
      pois_visited[[i]] <- data.table::rbindlist(visited)
      i <- i + 1
      setTxtProgressBar(pb, value = i)
    }
    close(pb)
    POIs_visited <- Filter(function(dt) nrow(dt) != 0, pois_visited) ## This is necessary to prevent a memory access violation!
    POIs_visited <- data.table::rbindlist(POIs_visited)
    POIs_visited <- POIs_visited[order(agent_id, time_index)]

    ## Convert the agents visitation sequence into an edge list. Save times according to the earliest arrival.
    cat("\rConverting visit sequences into an edge list...            ")
    all_visited <- lapply(unique(POIs_visited$agent_id), function(id){
      visit_seq <- POIs_visited[agent_id == id]$poi
      visit_time <- POIs_visited[agent_id == id]$time_index
      if (length(visit_seq) > 1){
        edge <- cbind(from = visit_seq[1:length(visit_seq)-1], to = visit_seq[2:length(visit_seq)])
        edge_time <- cbind(depart = visit_time[1:length(visit_seq)-1], arrive = visit_time[2:length(visit_seq)])
        return(data.frame(cbind(edge, edge_time)))
      }
      else return(NULL)
    })
    visit_el <- data.table::rbindlist(all_visited)

    ## Create the edge 'spells' (see ?networkDynamic)
    edge_spells <- as.matrix(visit_el[, .(onset=depart, terminus=arrive, tail=from, head = to)])
    el <- data.table::data.table(edge_spells)
    el <- el[order(onset, terminus)]

    ## Assume undirected discrete network
    minmax <- function(x, y) { c(min(x), max(y)) }
    keys <- apply(apply(combn(1:n, 2), 2, as.character), 2, function(key) paste0(key, collapse = ","))

    ## Fill in active times to create edge weights
    cat("\rCalculating edge weights...              ")
    time_rng <- range(c(el$onset, el$terminus))
    ew <- structure(vector(mode = "list", length = choose(n, 2)), names = keys)
    for (key in keys){ ew[[key]] <- rep(0L, (time_rng[2] - time_rng[1]) + 1L) } ## Initialize edge weights
    for (i in 1:nrow(el)){
      edge <- el[i, ]
      edge_key <- paste0(as.character(minmax(edge$tail, edge$head)), collapse = ",")
      edge_activity <- ew[[edge_key]]
      edge_activity[edge$onset:edge$terminus] <- edge_activity[edge$onset:edge$terminus] + 1L
      ew[[edge_key]] <- edge_activity
    }

    ## Convert the edge list into a set of edge toggles
    ## Assumes an undirected discrete network. This enables a smaller, compressed representation of the
    ## dynamic network that's easier to post-process.
    message("Creating edge state changes...")
    edge_toggles <- structure(vector(mode = "list", length = length(keys)), names = keys)
    for (key in keys){
      edge_ea <- rle(ew[[key]])
      v_ids <- as.integer(unlist(strsplit(key, ",")))
      if (all(edge_ea$values == 0)) next ## Don't consider inactive edges
      if (all(edge_ea$values == 1)) {
        edge_toggles[[key]] <- data.frame(time = time_rng[1], tail = v_ids[1], head = v_ids[2], direction = 1)
        next
      }
      idx <- cumsum(edge_ea$lengths) ## The start of every activation or deactivation
      toggles <- t(mapply(function(start, end) c(start, ew[[key]][start+1] > 0), c(1, idx[-length(idx)]), idx))
      # active_idx <- idx[which(edge_ea$values != 0)] ## The start of every activation
      # deactive_idx <- idx[which(edge_ea$values == 0)] ## The start of every deactivation
      #
      # deactive_start <- if (deactive_idx[1] == 1) c(1, act_idx[deactive_idx[-1] - 1]) else act_idx[deactive_idx]
      # deactive_end <- act_idx[deactive_idx]
      # toggle_changes <- matrix(rbind(cbind(deactive_start, rep(0L, length(deactive_start))),
      #                                cbind(deactive_end, rep(1L, length(deactive_end)))), ncol = 2)
      # toggle_changes <- matrix(toggle_changes[order(toggle_changes[,1]),], ncol = 2)

      ## Fix to make sure deactivations don't end when the network ends
      # if (toggle_changes[nrow(toggle_changes), 1] == time_rng[2] && nrow(toggle_changes) > 1){
      #   toggle_changes <- matrix(toggle_changes[1:(nrow(toggle_changes)-1),], ncol=2)
      # }


      ## Create the edge changes
      edge_act <- data.frame(time = toggles[, 1],                 ## The time increment the edge changed state
                             tail = rep(v_ids[1], nrow(toggles)), ## The tail vertex id
                             head = rep(v_ids[2], nrow(toggles)), ## The head vertex id
                             direction = toggles[, 2]             ## The direction/state (1 = active, 0 = inactive)
      )
      ## Fix the start to be active if it started active
      #if (deactive_start[1] > 1){  edge_act <- rbind(c(1L, v_ids[1], v_ids[2], 1L), edge_act) }
      edge_toggles[[key]] = edge_act
    }
    el_toggles <- as.matrix(data.table::rbindlist(edge_toggles))


    ## Final network
    message("Creating network...")
    base_net <- network::network.initialize(length(POIs), directed = F)
    final_net <- networkDynamic::networkDynamic(base_net, edge.changes = el_toggles)

    message("Adding dynamic edge weights...")
    for (key in keys){
      edge_ea <- rle(ew[[key]]) ## edge activity
      if (all(edge_ea$values == 0)) next ## Don't consider inactive edges
        v_ids <- as.integer(unlist(strsplit(key, ","))) ## vertex ids
        termini <- cumsum(edge_ea$lengths)
        onset <- c(1L, termini[-length(termini)])
        act_idx <- which(edge_ea$values != 0)
        networkDynamic::activate.edge.attribute(final_net, prefix = "weight",
                                                value = edge_ea$values[act_idx], ## The edge weight
                                                onset = onset[act_idx],          ## The onset of a new edge weight
                                                terminus = termini[act_idx]      ## The terminus of the edge weight
        )
    }
    return(final_net)
  }


  poi_centroids <- do.call(rbind, lapply(POIs, function(poi) sp::Polygon(poi)@labpt))
  final_net %v% "x" <- poi_centroids[, 1]
  final_net %v% "y" <- poi_centroids[, 2]














  change_times <- networkDynamic::get.change.times(final_net)
  mins <- c(5, 10, 15, 30, 1)
  for (n_min in mins){
    n_seq <- floor(max(change_times)/(5*60))
    switches <- seq(min(change_times), max(change_times), length.out = n_seq)
    onsets <- switches[1:length(switches)-1L]
    termini <- switches[2:length(switches)]
    nets <- mapply(function(start, end) { networkDynamic::network.collapse(final_net, onset = start, terminus = end) }, onsets, termini, SIMPLIFY = F)
  }



  ndtv::compute.animation(final_net, animation.mode = "kamadakawai",
                    slice.par=list(start=min(change_times), end=max(change_times), interval=1, aggregate.dur=1, rule='any'))
  ndtv::render.d3movie(final_net,
    edge.lwd = function(slice) {
      edge_weight <- (slice %e% "weight")
      if (length(edge_weight) == 0) return(0L)
      if (is.na(edge_weight)){ return(rep(1L, network::network.edgecount(slice))); }
      else {return(edge_weight/3) }
    })

  coords <- cbind(runif(53), runif(53))
  l <- network::network.layout.circle(nets[[1]])

  lapply()plot(nets[[1]], coord = coords)


  ## Populate the weight attribute in the network
  get.edge.activity(final_net, as.spellList = T)

  networkDynamic::activate.edge.attribute()

  data.table::data.table(edge_spells)[order(onset, terminus)]
  networkDynamic::network.extract(final_net, onset = 1, terminus = 10)


  networkDynamic::activate.edge.attribute(final_net, prefix = "weight",)
  load(file = "final_net.rdata")
  library("networkDynamic")
  library("ndtv")
  ndtv::render.d3movie(final_net, filename = "netDy.html")
  # plot(inp, col = poi_cl$cluster+1, asp = 1)

  # load("/Users/mpiekenbrock/WaCS/spatnet/data/CLIF.rda")
  # track_data <- CLIF[, .(lat, long, frame), by = trackID]

*/

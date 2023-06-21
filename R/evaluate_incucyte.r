bcm3.incucyte.get.likelihood <- function(bcm3, param.values) {
  res <- .C("bcm3_rbridge_incucyte_get_log_likelihood", bcm3$.cpp, param.values, 0, as.integer(0), PACKAGE="bcmrbridge")
  if (res[[4]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[4]]))
  }
  return(res[[3]])
}

bcm3.incucyte.get.simulated.trajectories <- function(bcm3, param.values, experiment_ix = 1) {
  traj_buffer <- rep(0.0, 5 * 100 * 11)
  res <- .C("bcm3_rbridge_incucyte_get_simulated_trajectories", bcm3$.cpp, as.numeric(param.values), traj_buffer, as.integer(experiment_ix), as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[6]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[6]]))
  }
  num_timepoints <- res[[5]]
  arr <- array(res[[3]], dim=c(num_timepoints,11,5))
  result <- list()
  result$cell_count <- arr[,,1]
  result$apoptotic_cell_count <- arr[,,2]
  result$debris <- arr[,,3]
  result$confluence <- arr[,,4]
  result$apoptosis_marker <- arr[,,5]
  return(result)
}

bcm3.incucyte.get.simulated.ctb <- function(bcm3, param.values, experiment_ix = 1) {
  ctb_buffer <- rep(0.0, 9)
  res <- .C("bcm3_rbridge_incucyte_get_simulated_ctb", bcm3$.cpp, as.numeric(param.values), ctb_buffer, as.integer(experiment_ix), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[5]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[5]]))
  }
  return(res[[3]])
}

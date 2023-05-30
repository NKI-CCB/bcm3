bcm3.ODE.get.likelihood <- function(bcm3, param.values) {
  res <- .C("bcm3_rbridge_ODE_get_log_likelihood", bcm3$.cpp, param.values, 0, as.integer(0), PACKAGE="bcmrbridge")
  if (res[[4]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[4]]))
  }
  return(res[[3]])
}

bcm3.ODE.get.simulated.trajectories <- function(bcm3, param.values) {
  traj_buffer <- rep(0.0, 33 * 6)
  res <- .C("bcm3_rbridge_ODE_get_simulated_trajectories", bcm3$.cpp, as.numeric(param.values), traj_buffer, as.integer(0), PACKAGE="bcmrbridge")
  if (res[[4]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[4]]))
  }
  return(res[[3]])
}

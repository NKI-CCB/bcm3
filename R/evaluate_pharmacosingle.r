
bcm3.pharmacosingle.get.likelihood <- function(bcm3, param.values) {
  res <- .C("bcm3_rbridge_pharmacosingle_get_log_likelihood", bcm3$.cpp, param.values, 0, as.integer(0), PACKAGE="bcmrbridge")
  if (res[[4]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[4]]))
  }
  return(res[[3]])
}

bcm3.pharmacosingle.get.observed.data <- function(bcm3) {
  retval <- list()
  max_nt <- 100
  traj_buffer <- rep(0.0, max_nt)
  time_buffer <- rep(0.0, max_nt)
  
  res <- .C("bcm3_rbridge_pharmacosingle_get_observed_data", bcm3$.cpp, traj_buffer, time_buffer,
            as.integer(max_nt), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[5]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[5]]))
  }
  
  ntimepoints <- res[[4]]
  
  retval$time <- res[[3]][1:ntimepoints]
  retval$data <- res[[2]][1:ntimepoints]
  return(retval)
}

bcm3.pharmacosingle.get.simulated.data <- function(bcm3, param.values) {
  retval <- list()
  max_nt <- 100
  traj_buffer <- rep(0.0, max_nt)
  time_buffer <- rep(0.0, max_nt)
  
  res <- .C("bcm3_rbridge_pharmacosingle_get_simulated_data", bcm3$.cpp, param.values, traj_buffer, time_buffer,
            as.integer(max_nt), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[6]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[6]]))
  }
  
  ntimepoints <- res[[5]]
  
  retval$time <- res[[4]][1:ntimepoints]
  retval$data <- res[[3]][1:ntimepoints]
  return(retval)
}

bcm3.pharmacosingle.get.simulated.trajectory <- function(bcm3, param.values, timepoints) {
  retval <- list()
  nt <- as.integer(length(timepoints))
  max_compartments <- 10
  traj_buffer <- rep(0.0, nt * max_compartments)
  conc_buffer <- rep(0.0, nt)
  
  res <- .C("bcm3_rbridge_pharmacosingle_get_simulated_trajectory", bcm3$.cpp, param.values, timepoints, nt, traj_buffer, conc_buffer,
            as.integer(max_compartments), as.integer(0), package="bcmrbridge")
  if (res[[8]] != 0) {
    stop(paste("bcm3 c++ bridge error:", res[[8]]))
  }
  
  ncompartments <- res[[7]]
  
  retval$central_compartment_concentrations <- res[[6]][1:nt]
  retval$trajectories <- array(res[[5]], c(ncompartments, nt))
  return(retval)
}
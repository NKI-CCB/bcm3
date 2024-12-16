
bcm3.PK.get.likelihood <- function(bcm3, param.values) {
  res <- .C("bcm3_rbridge_PK_get_log_likelihood", bcm3$.cpp, param.values, 0, as.integer(0), PACKAGE="bcmrbridge")
  if (res[[4]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[4]]))
  }
  return(res[[3]])
}

bcm3.PK.get.observed.data <- function(bcm3) {
  retval <- list()
  max_nt <- 10000
  traj_buffer <- rep(0.0, max_nt)
  time_buffer <- rep(0.0, max_nt)
  
  res <- .C("bcm3_rbridge_PK_get_observed_data", bcm3$.cpp, traj_buffer, time_buffer,
            as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[5]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[5]]))
  }
  
  ntimepoints <- res[[4]]
  
  retval$time <- res[[3]][1:ntimepoints]
  retval$data <- res[[2]][1:ntimepoints]
  return(retval)
}

bcm3.PK.get.simulated.data <- function(bcm3, param.values) {
  retval <- list()
  max_nt <- 10000
  max_compartments <- 5
  traj_buffer <- rep(0.0, max_nt * max_compartments)
  conc_buffer <- rep(0.0, max_nt)
  time_buffer <- rep(0.0, max_nt)
  
  res <- .C("bcm3_rbridge_PK_get_simulated_data", bcm3$.cpp, param.values, traj_buffer, conc_buffer, time_buffer,
            as.integer(0), as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[8]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[8]]))
  }
  
  ntimepoints <- res[[6]]
  
  retval$time <- res[[5]][1:ntimepoints]
  retval$data <- res[[4]][1:ntimepoints]
  return(retval)
}

bcm3.PK.get.simulated.trajectories <- function(bcm3, param.values) {
  retval <- list()
  max_nt <- 10000
  max_compartments <- 5
  traj_buffer <- rep(0.0, max_nt * max_compartments)
  conc_buffer <- rep(0.0, max_nt)
  time_buffer <- rep(0.0, max_nt)
  
  res <- .C("bcm3_rbridge_PK_get_simulated_data", bcm3$.cpp, param.values, traj_buffer, conc_buffer, time_buffer,
            as.integer(0), as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[8]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[8]]))
  }
  
  ntimepoints <- res[[6]]
  ncompartments <- res[[7]]
  
  retval$time <- res[[5]][1:ntimepoints]
  retval$trajectories <- array(res[[3]], c(ntimepoints, ncompartments))
  return(retval)
}
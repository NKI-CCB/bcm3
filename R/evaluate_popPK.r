
bcm3.popPK.get.likelihood <- function(bcm3, param.values) {
  res <- .C("bcm3_rbridge_popPK_get_log_likelihood", bcm3$.cpp, param.values, 0, as.integer(0), PACKAGE="bcmrbridge")
  if (res[[4]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[4]]))
  }
  return(res[[3]])
}

bcm3.popPK.get.observed.data <- function(bcm3) {
  retval <- list()
  max_nt <- 10000
  max_patients <- 100
  traj_buffer <- rep(0.0, max_nt * max_patients)
  time_buffer <- rep(0.0, max_nt)
  
  res <- .C("bcm3_rbridge_popPK_get_observed_data", bcm3$.cpp, traj_buffer, time_buffer,
            as.integer(0), as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[6]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[6]]))
  }
  
  npatients <- res[[4]]
  ntimepoints <- res[[5]]
  
  retval$time <- res[[3]][1:ntimepoints]
  retval$data <- array(res[[2]], c(ntimepoints, npatients))
  return(retval)
}

bcm3.popPK.get.simulated.data <- function(bcm3, param.values) {
  retval <- list()
  max_nt <- 10000
  max_patients <- 100
  max_compartments <- 5
  traj_buffer <- rep(0.0, max_nt * max_patients * max_compartments)
  conc_buffer <- rep(0.0, max_nt * max_patients)
  time_buffer <- rep(0.0, max_nt)
  
  res <- .C("bcm3_rbridge_popPK_get_simulated_data", bcm3$.cpp, param.values, traj_buffer, conc_buffer, time_buffer,
            as.integer(0), as.integer(0), as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[9]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[9]]))
  }
  
  npatients <- res[[6]]
  ntimepoints <- res[[7]]
  
  retval$time <- res[[5]][1:ntimepoints]
  retval$data <- array(res[[4]], c(ntimepoints, npatients))
  return(retval)
}

bcm3.popPK.get.simulated.trajectories <- function(bcm3, param.values) {
  retval <- list()
  max_nt <- 10000
  max_patients <- 100
  max_compartments <- 5
  traj_buffer <- rep(0.0, max_nt * max_patients * max_compartments)
  conc_buffer <- rep(0.0, max_nt * max_patients)
  time_buffer <- rep(0.0, max_nt)
  
  res <- .C("bcm3_rbridge_popPK_get_simulated_data", bcm3$.cpp, param.values, traj_buffer, conc_buffer, time_buffer,
            as.integer(0), as.integer(0), as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[9]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[9]]))
  }
  
  npatients <- res[[6]]
  ntimepoints <- res[[7]]
  ncompartments <- res[[8]]
  
  retval$time <- res[[5]][1:ntimepoints]
  retval$trajectories <- array(res[[3]], c(ntimepoints, ncompartments, npatients))
  return(retval)
}
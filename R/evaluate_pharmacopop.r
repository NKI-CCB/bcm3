
bcm3.pharmacopop.get.likelihood <- function(bcm3, param.values) {
  res <- .C("bcm3_rbridge_pharmacopop_get_log_likelihood", bcm3$.cpp, param.values, 0, as.integer(0), PACKAGE="bcmrbridge")
  if (res[[4]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[4]]))
  }
  return(res[[3]])
}

bcm3.pharmacopop.get.observed.data <- function(bcm3, patient_ix) {
  retval <- list()
  max_nt <- 100
  traj_buffer <- rep(0.0, max_nt)
  time_buffer <- rep(0.0, max_nt)
  
  res <- .C("bcm3_rbridge_pharmacopop_get_observed_data", bcm3$.cpp, as.integer(patient_ix), traj_buffer, time_buffer,
            as.integer(max_nt), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[6]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[6]]))
  }
  
  ntimepoints <- res[[5]]
  
  retval$time <- res[[4]][1:ntimepoints]
  retval$data <- res[[3]][1:ntimepoints]
  return(retval)
}

bcm3.pharmacopop.get.num.patients <- function(bcm3) {
  retval <- list()
  max_nt <- 100
  traj_buffer <- rep(0.0, max_nt)
  time_buffer <- rep(0.0, max_nt)
  
  res <- .C("bcm3_rbridge_pharmacopop_get_num_patients", bcm3$.cpp, as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[3]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[6]]))
  }
  num_patients <- res[[2]]
  return(num_patients)
}

bcm3.pharmacopop.get.simulated.trajectory <- function(bcm3, param.values, timepoints, patient_ix) {
  retval <- list()
  nt <- as.integer(length(timepoints))
  max_compartments <- 10
  traj_buffer <- rep(0.0, nt * max_compartments)
  conc_buffer <- rep(0.0, nt)
  
  res <- .C("bcm3_rbridge_pharmacopop_get_simulated_trajectory", bcm3$.cpp, param.values, timepoints, nt,
            as.integer(patient_ix), traj_buffer, conc_buffer, as.integer(max_compartments), as.integer(0), package="bcmrbridge")
  if (res[[9]] != 0) {
    stop(paste("bcm3 c++ bridge error:", res[[9]]))
  }
  
  ncompartments <- res[[8]]
  
  retval$central_compartment_concentrations <- res[[7]][1:nt]
  retval$trajectories <- array(res[[5]], c(ncompartments, nt))
  return(retval)
}

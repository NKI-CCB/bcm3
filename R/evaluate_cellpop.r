bcm3.cellpop.get.likelihood <- function(bcm3, experiment, param.values) {
  res <- .C("bcm3_rbridge_cellpop_get_log_likelihood", bcm3$.cpp, param.values, 0, as.integer(0), PACKAGE="bcmrbridge")
  if (res[[4]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[4]]))
  }
  return(res[[3]])
}

bcm3.cellpop.get.num.species <- function(bcm3, experiment) {
  res <- .C("bcm3_rbridge_cellpop_get_num_species", bcm3$.cpp, as.character(experiment), as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[4]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[4]]))
  }
  return(res[[3]])
}

bcm3.cellpop.get.species.name <- function(bcm3, experiment, species_ix) {
  res <- .C("bcm3_rbridge_cellpop_get_species_name", bcm3$.cpp, as.character(experiment), as.integer(species_ix-1), "", as.integer(0), PACKAGE="bcmrbridge")
  if (res[[5]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[5]]))
  }
  return(res[[4]])
}

bcm3.cellpop.get.simulated.trajectories <- function(bcm3, experiment, param.values) {
  ns <- bcm3.cellpop.get.num.species(bcm3, experiment)
  max_nt <- 500
  max_cells <- 500
  traj_buffer <- rep(0.0, ns * max_nt * max_cells)
  time_buffer <- rep(0.0, max_nt)
  res <- .C("bcm3_rbridge_cellpop_get_simulated_trajectories", bcm3$.cpp, as.character(experiment), as.numeric(param.values), traj_buffer, time_buffer, as.integer(0), as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[8]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[8]]))
  }

  ntimepoints <- res[[7]]
  ncells <- res[[6]]

  retval <- list()
  retval$time <- res[[5]][1:ntimepoints]
  retval$cells <- array(res[[4]], c(ns, ntimepoints, ncells))
  return(retval)
}

bcm3.cellpop.get.observed.data <- function(bcm3, experiment) {
  res <- .C("bcm3_rbridge_cellpop_get_num_data", bcm3$.cpp, as.character(experiment), as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[4]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[4]]))
  }
  nd <- res[[3]]
  
  retval <- list()
  for (i in 1:nd) {
    retval[[i]] <- list()
    max_nt <- 500
    max_cells <- 500
    max_markers <- 5
    traj_buffer <- rep(0.0, max_nt * max_cells * max_markers)
    time_buffer <- rep(0.0, max_nt)
    
    res <- .C("bcm3_rbridge_cellpop_get_observed_data", bcm3$.cpp, as.character(experiment), as.integer(i-1), traj_buffer, time_buffer,
              as.integer(0), as.integer(0), as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
    if (res[[9]] != 0) {
      stop(paste("BCM3 C++ bridge error:", res[[9]]))
    }
    
    ncells <- res[[6]]
    ntimepoints <- res[[7]]
    nmarkers <- res[[8]]
    
    retval[[i]]$time <- res[[5]][1:ntimepoints]
    retval[[i]]$data <- array(res[[4]], c(ntimepoints, ncells, nmarkers))
  }
  return(retval)
}

bcm3.cellpop.get.simulated.data <- function(bcm3, experiment, param.values) {
  res <- .C("bcm3_rbridge_cellpop_get_num_data", bcm3$.cpp, as.character(experiment), as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[4]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[4]]))
  }
  nd <- res[[3]]
  
  retval <- list()
  for (i in 1:nd) {
    retval[[i]] <- list()
    max_nt <- 500
    max_cells <- 500
    max_markers <- 5
    traj_buffer <- rep(0.0, max_nt * max_cells * max_markers)
    time_buffer <- rep(0.0, max_nt)
    
    res <- .C("bcm3_rbridge_cellpop_get_simulated_data", bcm3$.cpp, as.character(experiment), as.numeric(param.values), as.integer(i-1), traj_buffer, time_buffer,
              as.integer(0), as.integer(0), as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
    if (res[[10]] != 0) {
      stop(paste("BCM3 C++ bridge error:", res[[10]]))
    }
    
    ncells <- res[[7]]
    ntimepoints <- res[[8]]
    nmarkers <- res[[9]]
    
    retval[[i]]$time <- res[[6]][1:ntimepoints]
    retval[[i]]$data <- array(res[[5]], c(ntimepoints, ncells, nmarkers))
  }
  return(retval)
}

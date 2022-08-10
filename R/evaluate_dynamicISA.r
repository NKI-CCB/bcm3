bcm3.dynamicISA.get.num.data <- function(bcm3, experiment) {
  res <- .C("bcm3_rbridge_dynamicISA_get_num_data", bcm3$.cpp, experiment, as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[4]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[4]]))
  }
  return(res[[3]])
}

bcm3.dynamicISA.get.num.conditions <- function(bcm3, experiment) {
  res <- .C("bcm3_rbridge_dynamicISA_get_num_conditions", bcm3$.cpp, experiment, as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[4]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[4]]))
  }
  return(res[[3]])
}

bcm3.dynamicISA.get.observed.data <- function(bcm3, experiment, data.ix) {
  nc <- bcm3.dynamicISA.get.num.conditions(bcm3, experiment)
  
  res <- .C("bcm3_rbridge_dynamicISA_get_num_replicates", bcm3$.cpp, experiment, as.integer(data.ix-1), as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[5]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[4]]))
  }
  nr <- res[[4]]

  buffer <- rep(0.0, nr*nc)
  res <- .C("bcm3_rbridge_dynamicISA_get_observed_data", bcm3$.cpp, experiment, as.integer(data.ix-1), buffer, as.integer(0), PACKAGE="bcmrbridge")
  if (res[[5]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[5]]))
  }
  return(matrix(res[[4]], nr, nc))
}

bcm3.dynamicISA.get.modeled.data <- function(bcm3, experiment, data.ix, param.values) {
  nc <- bcm3.dynamicISA.get.num.conditions(bcm3, experiment)
  
  buffer <- rep(0.0, nc)
  res <- .C("bcm3_rbridge_dynamicISA_get_modeled_data", bcm3$.cpp, experiment, param.values, as.integer(data.ix-1), buffer, as.integer(0), PACKAGE="bcmrbridge")
  if (res[[6]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[6]]))
  }
  return(res[[5]])
}

bcm3.dynamicISA.get.modeled.activities <- function(bcm3, experiment, param.values) {
  nc <- bcm3.dynamicISA.get.num.conditions(bcm3, experiment)
  
  res <- .C("bcm3_rbridge_dynamicISA_get_num_signaling_molecules", bcm3$.cpp, as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[3]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[3]]))
  }
  nsm <- res[[2]]
  
  buffer <- rep(0.0, nc * nsm)
  res <- .C("bcm3_rbridge_dynamicISA_get_modeled_activities", bcm3$.cpp, experiment, param.values, buffer, as.integer(0), PACKAGE="bcmrbridge")
  if (res[[5]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[5]]))
  }
  return(matrix(res[[4]], nsm, nc))
}

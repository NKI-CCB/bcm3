bcm3.fISA.get.num.data <- function(bcm3, experiment) {
  res <- .C("bcm3_rbridge_fISA_get_num_data", bcm3$.cpp, experiment, as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[4]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[4]]))
  }
  return(res[[3]])
}

bcm3.fISA.get.num.cell.lines <- function(bcm3, experiment) {
  res <- .C("bcm3_rbridge_fISA_get_num_cell_lines", bcm3$.cpp, experiment, as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[4]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[4]]))
  }
  return(res[[3]])
}

bcm3.fISA.get.cell.line.names <- function(bcm3, experiment) {
  nc <- bcm3.fISA.get.num.cell.lines(bcm3, experiment)
  cell_lines <- list()
  for (i in 1:nc) {
    res <- .C("bcm3_rbridge_fISA_get_cell_line_name", bcm3$.cpp, experiment, as.integer(i-1), "", as.integer(0), PACKAGE="bcmrbridge")
    if (res[[5]] != 0) {
      stop(paste("BCM3 C++ bridge error:", res[[5]]))
    }
    cell_lines[[i]] <- res[[4]]
  }
  return(unlist(cell_lines))
}

bcm3.fISA.get.observed.data <- function(bcm3, experiment, data.ix) {
  nc <- bcm3.fISA.get.num.cell.lines(bcm3, experiment)
  
  res <- .C("bcm3_rbridge_fISA_get_num_replicates", bcm3$.cpp, experiment, as.integer(data.ix-1), as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[5]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[5]]))
  }
  nr <- res[[4]]

  buffer <- rep(0.0, nr*nc)
  res <- .C("bcm3_rbridge_fISA_get_observed_data", bcm3$.cpp, experiment, as.integer(data.ix-1), buffer, as.integer(0), PACKAGE="bcmrbridge")
  if (res[[5]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[5]]))
  }
  return(matrix(res[[4]], nr, nc))
}

bcm3.fISA.get.likelihood <- function(bcm3, param.values) {
  res <- .C("bcm3_rbridge_fISA_get_log_likelihood", bcm3$.cpp, param.values, 0, as.integer(0), PACKAGE="bcmrbridge")
  if (res[[4]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[4]]))
  }
  return(res[[3]])
}

bcm3.fISA.get.modeled.data <- function(bcm3, experiment, data.ix, param.values) {
  nc <- bcm3.fISA.get.num.cell.lines(bcm3, experiment)
  
  buffer <- rep(0.0, nc)
  res <- .C("bcm3_rbridge_fISA_get_modeled_data", bcm3$.cpp, experiment, param.values, as.integer(data.ix-1), buffer, as.integer(0), PACKAGE="bcmrbridge")
  if (res[[6]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[6]]))
  }
  return(res[[5]])
}

bcm3.fISA.get.modeled.activities <- function(bcm3, experiment, param.values) {
  nc <- bcm3.fISA.get.num.cell.lines(bcm3, experiment)
  
  res <- .C("bcm3_rbridge_fISA_get_num_signaling_molecules", bcm3$.cpp, as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[3]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[3]]))
  }
  nsm <- res[[2]]
  
  buffer <- rep(0.0, nc * nsm)
  res <- .C("bcm3_rbridge_fISA_get_modeled_activities", bcm3$.cpp, experiment, param.values, buffer, as.integer(0), PACKAGE="bcmrbridge")
  if (res[[5]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[5]]))
  }
  return(matrix(res[[4]], nsm, nc))
}

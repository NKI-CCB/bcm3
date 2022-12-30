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
  parent_buffer <- as.integer(rep(-1, max_cells))
  res <- .C("bcm3_rbridge_cellpop_get_simulated_trajectories", bcm3$.cpp, as.character(experiment), as.numeric(param.values), traj_buffer, time_buffer, parent_buffer, as.integer(0), as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[9]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[9]]))
  }

  ntimepoints <- res[[8]]
  ncells <- res[[7]]

  retval <- list()
  retval$time <- res[[5]][1:ntimepoints]
  retval$cells <- array(res[[4]], c(ns, ntimepoints, ncells))
  retval$parents <- res[[6]][1:ncells]
  
  dimnames(retval$cells)[[1]] <- 1:ns
  for (i in 1:ns) {
    dimnames(retval$cells)[[1]][i] <- bcm3.cellpop.get.species.name(bcm3, experiment, i)
  }
  
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

bcm3.cellpop.get.matched.simulation <- function(bcm3, experiment, param.values) {
  res <- .C("bcm3_rbridge_cellpop_get_num_data", bcm3$.cpp, as.character(experiment), as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[4]] != 0) {
    stop(paste("BCM3 C++ bridge error:", res[[4]]))
  }
  nd <- res[[3]]
  ns <- bcm3.cellpop.get.num.species(bcm3, experiment)
  
  retval <- list()
  for (i in 1:nd) {
    retval[[i]] <- list()
    max_nt <- 500
    max_cells <- 500
    traj_buffer <- rep(0.0, ns * max_nt * max_cells)
    time_buffer <- rep(0.0, max_nt)
    
    res <- .C("bcm3_rbridge_cellpop_get_matched_simulation", bcm3$.cpp, as.character(experiment), as.numeric(param.values), as.integer(i-1), traj_buffer, time_buffer,
              as.integer(0), as.integer(0), as.integer(0), PACKAGE="bcmrbridge")
    if (res[[9]] != 0) {
      stop(paste("BCM3 C++ bridge error:", res[[10]]))
    }
    
    ncells <- res[[7]]
    ntimepoints <- res[[8]]
    
    retval[[i]]$time <- res[[6]][1:ntimepoints]
    retval[[i]]$cells <- array(res[[5]], c(ns, ntimepoints, ncells))
  }
  return(retval)
}

bcm3.cellpop.posterior.predictive.plot <- function(bcm3, experiment_ix, temperature_ix, numppdsamples, pdffilename, pdfwidth, pdfheight, xlim=NULL, time_unit="minutes", pch=20)
{
  model <- bcm3
  sample <- model$posterior$samples[,temperature_ix,dim(model$posterior$samples)[3]]
  names(sample) <- model$variables
  obsdata <- bcm3.cellpop.get.observed.data(model, model$likelihood$experiments[[experiment_ix]]$name)
  
  sample_ix <- (dim(model$posterior$samples)[3]/2+1):(dim(model$posterior$samples)[3])
  ppd_sample_ix <- sample(sample_ix, numppdsamples)
  
  merged <- list()
  for (j in 1:length(obsdata)) {
    merged[[j]] <- list()
    for (l in 1:dim(obsdata[[j]]$data)[2]) {
      merged[[j]][[l]] <- matrix(NA, numppdsamples, length(obsdata[[j]]$time))
    }
  }
  for (i in 1:numppdsamples) {
    sample <- model$posterior$samples[,temperature_ix,ppd_sample_ix[i]]
    names(sample) <- model$variables
    simdata <- bcm3.cellpop.get.simulated.data(model, model$likelihood$experiments[[experiment_ix]]$name, sample)
    
    for (j in 1:length(obsdata)) {
      for (l in 1:dim(obsdata[[j]]$data)[2]) {
        vals <- simdata[[j]]$data[,l,1]
        vals[vals == -Inf] <- -2
        merged[[j]][[l]][i,] <- vals
      }
    }
  }
  
  plot_count <- 0
  for (k in 1:length(obsdata)) {
    if (model$likelihood$experiments[[experiment_ix]]$data[[k]]$type == "duration") {
      if (sum(!is.na(obsdata[[k]]$data[1,,1])) > 0) {
        plot_count <- plot_count + 1
      }
    } else {
      plot_count <- plot_count + dim(obsdata[[k]]$data)[2]
    }
  }
  
  pdf_tile(pdffilename, pdfwidth, pdfheight, plot_count)
  
  for (k in 1:length(obsdata)) {
    stdev_ix <- match(model$likelihood$experiments[[experiment_ix]]$data[[k]]$stdev, model$variables)
    if (is.na(stdev_ix)) {
      sds <- as.numeric(model$likelihood$experiments[[experiment_ix]]$data[[k]]$stdev)
    } else {
      sds <- 10^as.numeric(model$posterior$samples[stdev_ix, temperature_ix, ppd_sample_ix])
    }
    if (model$likelihood$experiments[[experiment_ix]]$data[[k]]$type == "duration") {
      if (sum(!is.na(obsdata[[k]]$data[1,,1])) > 0) {
        tmp <- matrix(unlist(merged[[k]]), nrow=numppdsamples)
        ppd_barplot(model, tmp, obsdata[[k]]$data[1,,1], labels=1:ncol(tmp),
                    error_model = "normal", sd_samples = sds * model$likelihood$experiments[[experiment_ix]]$data[[k]]$stdev_multiplication_factor,
                    xlab="Cells", ylab="NEBD to anaphase time (seconds)", ylim=c(0,6000), pointsize=1.0)
      }
    } else {
      if (time_unit == "minutes") {
        t <- obsdata[[k]]$time / 60
      } else if (time_unit == "hours") {
        t <- obsdata[[k]]$time / 3600
      } else {
        t <- obsdata[[k]]$time
      }
      for (cell_ix in 1:dim(obsdata[[k]]$data)[2]) {
        lower <- rep(NA, length(obsdata[[k]]$time))
        upper <- rep(NA, length(obsdata[[k]]$time))
        for (i in 1:length(obsdata[[k]]$time)) {
          x <- rnorm(length(merged[[k]][[cell_ix]][,i]) * 100, merged[[k]][[cell_ix]][,i], sds * model$likelihood$experiments[[experiment_ix]]$data[[k]]$stdev_multiplication_factor)
          lower[i] <- quantile(x, 0.05, na.rm=T)
          upper[i] <- quantile(x, 0.95, na.rm=T)
        }
        ylim <- c(min(-0.1, min(obsdata[[k]]$data-0.1, lower, na.rm=T)),
                  max(1.1, max(obsdata[[k]]$data, upper, na.rm=T)+0.1, lower, na.rm=T))
        if (is.null(xlim)) {
          use_xlim <- range(t)
        } else {
          use_xlim <- xlim
        }
        plot(t, obsdata[[k]]$data[,cell_ix,1], pch=pch, ylab=model$likelihood$experiments[[experiment_ix]]$data[[k]]$data_name, xlab="Time (seconds)", xlim=xlim, ylim=ylim, yaxt='n', type='n')
        axis(2, at=seq(0,1,by=0.5))
        rug(seq(-0.1,1.1,by=0.1), side=2, ticksize=-0.01)
        ppd_color <- rgb(43,131,186, max=255)
        have_plot <- which(!is.na(lower) & !is.na(upper))
        polygon(c(t[have_plot], rev(t[have_plot])), c(lower[have_plot], rev(upper[have_plot])), col=adjustcolor(ppd_color, alpha.f = 0.4), border=F)
        lines(t[have_plot], lower[have_plot], col=ppd_color, lwd=2)
        lines(t[have_plot], upper[have_plot], col=ppd_color, lwd=2)
        points(t, obsdata[[k]]$data[,cell_ix,1], pch=pch)
      }
    }
  }
  
  res <- dev.off()
}
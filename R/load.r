library(XML)
library(hdf5r)

bcm3.load.results <- function(base_folder, output_folder, prior_file="prior.xml", likelihood_file="likelihood.xml", output_filename="output.nc", load_sampler_adaptation = T)
{
  model <- bcm3.load(base_folder, prior_file, likelihood_file)
  model$output_folder <- paste(base_folder, "/", output_folder, sep="")
  
  # Posterior
  output_file <- H5File$new(paste(model$output_folder, "/", output_filename, sep=""), 'r')
  
  model$posterior <- list()
  model$posterior$temperatures <- output_file[["samples/temperature"]][]
  model$posterior$samples <- output_file[["samples/variable_values"]][,,,drop=F]
  if ("weights" %in% names(output_file[["samples"]])) {
    model$posterior$weights <- output_file[["samples/weights"]][,,drop=F]
  } else {
    model$posterior$weights <- matrix(1, dim(model$posterior$samples)[2], dim(model$posterior$samples)[3])
  }
  if (length(output_file[["samples/log_prior"]]$dims) == 1) {
    ntemps <- dim(model$posterior$samples)[2]
    nsamples <- dim(model$posterior$samples)[3]
    model$posterior$lprior <- matrix(NA, ntemps, nsamples)
    model$posterior$llikelihood <- matrix(NA, ntemps, nsamples)
    model$posterior$lprior[ntemps,] <- output_file[["samples/log_prior"]][]
    model$posterior$llikelihood[ntemps,] <- output_file[["samples/log_likelihood"]][]
  } else {
    model$posterior$lprior <- output_file[["samples/log_prior"]][,,drop=F]
    model$posterior$llikelihood <- output_file[["samples/log_likelihood"]][,,drop=F]
  }
  
  fill_value <- output_file[["samples/variable_values"]]$get_fill_value()
  model$posterior$samples[model$posterior$samples == fill_value] <- NA
  model$posterior$weights[model$posterior$weights == fill_value] <- NA
  model$posterior$lprior[model$posterior$lprior == fill_value] <- NA
  model$posterior$llikelihood[model$posterior$llikelihood == fill_value] <- NA
  
  model$posterior$lposterior <- model$posterior$lprior + model$posterior$llikelihood
  model$posterior$lfracposterior <- matrix(NA, nrow(model$posterior$lprior), ncol(model$posterior$lprior))
  for (i in 1:length(model$posterior$temperatures)) {
    model$posterior$lfracposterior[i,] <- model$posterior$lprior[i,] + model$posterior$temperatures[i] * model$posterior$llikelihood[i,]
  }
  
  output_file$close_all()
  
  if (load_sampler_adaptation) {
    sampler_adaptation_fn <- paste(model$output_folder, "/sampler_adaptation.nc", sep="")
    if (file.exists(sampler_adaptation_fn)) {
      model$sampler_adaptation <- load.netcdf.bundler.data(sampler_adaptation_fn)
    } else {
      model$sampler_adaptation <- NULL
    }
  } else {
    model$sampler_adaptation <- NULL
  }
  
  model$AIC <- 2 * model$nvar - 2 * max(model$posterior$llikelihood, na.rm=T)
  
  return(structure(model, class = "bcm3"))
}

bcm3.load <- function(base_folder, prior_file="prior.xml", likelihood_file="likelihood.xml")
{
  model <- list()
  model$base_folder = base_folder
  
  # Prior
  model$prior <- list()
  model$prior$file_name <- prior_file
  model$prior$variable_attrs <- list()
  varxml <- xmlTreeParse(paste(base_folder, "/", prior_file, sep=""))
  variables <- xmlRoot(varxml)["variable", all=T]
  variable_attrs <- lapply(variables, xmlAttrs)
  model$variables <- c()
  var_ix <- 1
  for (i in 1:length(variables)) {
    if ("repeat" %in% names(variable_attrs[[i]])) {
      repeats = as.numeric(variable_attrs[[i]]["repeat"])
      for (k in 1:repeats) {
        model$prior$variable_attrs[[var_ix]] <- variable_attrs[[i]]
        model$variables[var_ix] <- paste(variable_attrs[[i]]["name"], k, sep="_")
        var_ix <- var_ix + 1
      }
    } else {
      model$prior$variable_attrs[[var_ix]] <- variable_attrs[[i]]
      model$variables[var_ix] <- variable_attrs[[i]]["name"]
      var_ix <- var_ix + 1
    }
  }
  model$nvar <- length(model$variables)
  
  # Likelihood
  model$likelihood <- list()
  model$likelihood$file_name <- likelihood_file
  likelihoodxml <- xmlTreeParse(paste(base_folder, "/", likelihood_file, sep=""))
  model$likelihood$type <- xmlAttrs(xmlRoot(likelihoodxml))["type"]
  experiments <- xmlRoot(likelihoodxml)["experiment", all=T]
  if (length(experiments) > 0) {
    model$likelihood$experiments <- list()
    for (i in 1:length(experiments)) {
      model$likelihood$experiments[[i]] <- list()
      model$likelihood$experiments[[i]]$name <- xmlAttrs(experiments[[i]])["name"]
      data <- experiments[[i]]["data", all=T]
      if (length(data) > 0) {
        model$likelihood$experiments[[i]]$data <- list()
        for (j in 1:length(data)) {
          model$likelihood$experiments[[i]]$data[[j]] <- list()
          model$likelihood$experiments[[i]]$data[[j]]$type <- xmlAttrs(data[[j]])["type"]
          if (is.na(model$likelihood$experiments[[i]]$data[[j]]$type)) {
            model$likelihood$experiments[[i]]$data[[j]]$type <- "time_course"
          }
          model$likelihood$experiments[[i]]$data[[j]]$species_name <- xmlAttrs(data[[j]])["species_name"]
          model$likelihood$experiments[[i]]$data[[j]]$data_name <- xmlAttrs(data[[j]])["data_name"]
          model$likelihood$experiments[[i]]$data[[j]]$stdev <- xmlAttrs(data[[j]])["stdev"]
          model$likelihood$experiments[[i]]$data[[j]]$stdev_multiplication_factor <- as.numeric(xmlAttrs(data[[j]])["stdev_multiplication_factor"])
          model$likelihood$experiments[[i]]$data[[j]]$base_scale_sd_suffix <- xmlAttrs(data[[j]])["base_scale_sd_suffix"]
          if (is.na(model$likelihood$experiments[[i]]$data[[j]]$stdev_multiplication_factor)) {
            model$likelihood$experiments[[i]]$data[[j]]$stdev_multiplication_factor <- 1.0
          }
        }
      }
    }
  }
  
  return(structure(model, class = "bcm3"))
}

load.netcdf.bundler.data <- function(netcdf_bundler_filename) {
  file <- H5File$new(netcdf_bundler_filename, 'r')
  res <- list()
  for (i in 1:length(names(file))) {
    groupname <- names(file)[i]
    res[[groupname]] <- load.netcdf.bundler.group(file[[groupname]])
  }
  file$close_all()
  return(res)
}

load.netcdf.bundler.group <- function(group) {
  result <- list()
  for (j in 1:length(names(group))) {
    elementname <- names(group)[j]
    
    if (class(group[[elementname]])[1] == "H5Group") {
      result[[elementname]] <- load.netcdf.bundler.group(group[[elementname]])
    } else {
      varname <- elementname
      if (substr(varname, nchar(varname)-3,nchar(varname)) != "dim1" &&
          substr(varname, nchar(varname)-3,nchar(varname)) != "dim2") {
        if (length(group[[varname]]$dims) == 1) {
          result[[varname]] <- group[[varname]][]
        } else if (length(group[[varname]]$dims) == 2) {
          result[[varname]] <- group[[varname]][,]
        } else {
          stop()
        }
      }
    }
  }
  return(result)
}

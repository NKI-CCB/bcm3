library(XML)
library(hdf5r)

bcm3.load.results <- function(base_folder, output_folder, prior_file="prior.xml", likelihood_file="likelihood.xml", output_filename="output.nc", load_sampler_adaptation = T)
{
  model <- bcm3.load(base_folder, prior_file, likelihood_file)
  model$output_folder <- paste(base_folder, "/", output_folder, sep="")
  
  # Posterior
  output_file <- H5File$new(paste(model$output_folder, "/", output_filename, sep=""), 'r')
  
  model$posterior <- list()
  model$posterior$samples <- output_file[["samples/variable_values"]][,,]
  if (length(output_file[["samples/log_prior"]]$dims) == 1) {
    ntemps <- dim(model$posterior$samples)[2]
    nsamples <- dim(model$posterior$samples)[3]
    model$posterior$lprior <- matrix(NA, ntemps, nsamples)
    model$posterior$llikelihood <- matrix(NA, ntemps, nsamples)
    model$posterior$lprior[ntemps,] <- output_file[["samples/log_prior"]][]
    model$posterior$llikelihood[ntemps,] <- output_file[["samples/log_likelihood"]][]
  } else {
    model$posterior$lprior <- output_file[["samples/log_prior"]][,]
    model$posterior$llikelihood <- output_file[["samples/log_likelihood"]][,]
  }
  #colnames(model$posterior$samples) <- output_file[["samples/variable"]][]
  model$posterior$temperatures <- output_file[["samples/temperature"]][]
  
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
  varxml <- xmlTreeParse(paste(base_folder, "/", prior_file, sep=""))
  variables <- xmlRoot(varxml)["variable", all=T]
  model$prior$variable_attrs <- lapply(variables, xmlAttrs)
  model$variables <- rep("", length(variables))
  for (i in 1:length(variables)) {
    model$variables[i] <- model$prior$variable_attrs[[i]]["name"]
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
    group <- names(file)[i]
    res[[group]] <- list()
    
    for (j in 1:length(names(file[[group]]))) {
      varname <- names(file[[group]])[j]
      if (substr(varname, nchar(varname)-3,nchar(varname)) != "dim1" &&
          substr(varname, nchar(varname)-3,nchar(varname)) != "dim2") {
        if (length(file[[group]][[varname]]$dims) == 1) {
          res[[group]][[varname]] <- file[[group]][[varname]][]
        } else if (length(file[[group]][[varname]]$dims) == 2) {
          res[[group]][[varname]] <- file[[group]][[varname]][,]
        } else {
          stop()
        }
      }
    }
  }
  file$close_all()
  return(res)
}
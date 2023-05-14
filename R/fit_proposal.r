suppressPackageStartupMessages(library(hdf5r))
suppressPackageStartupMessages(library(EMMIXmfa))
suppressPackageStartupMessages(library(mclust))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1)
tmpfile <- args[1]
if (length(args) > 1) {
  log_info <- T
} else {
  log_info <- F
}

file <- h5file(tmpfile)
samples <- t(file[["history"]][["samples"]][,])
h5close(file)

ndim <- ncol(samples)

ncomponents <- c(1,2,3,5,8)
ncomponents <- ncomponents[ncomponents < sqrt(nrow(samples))]

if (ndim > 1) {
  nfactors <- c(1, 1)
  for (i in 1:ndim) {
    nfactors <- c(nfactors, nfactors[i] + nfactors[i+1])
  }
  nfactors <- unique(nfactors[nfactors <= (ndim-1)])
} else {
  nfactors <- c(1)
}

mtfaress <- list()
BICs <- matrix(NA, length(ncomponents), length(nfactors))
minbic <- Inf
mini <- NA
minj <- NA
for (i in 1:length(ncomponents)) {
  mtfaress[[i]] <- list()
  for (j in 1:length(nfactors)) {
    mtfa <- NULL
    output <- capture.output(
      result <- try ({
        mtfa <- EMMIXmfa::mtfa(samples, ncomponents[i], nfactors[j], sigma_type = "unique", D_type = "common", tol=1e-4, nkmeans = 5, nrandom = 5, conv_measure='ratio')
      }, silent = T)
    )
    
    if (!is.null(mtfa)) {
      mtfa <- result
      mtfaress[[i]][[j]] <- mtfa
      BICs[i,j] <- mtfaress[[i]][[j]]$BIC
      
      if (!is.na(mtfaress[[i]][[j]]$BIC) && mtfaress[[i]][[j]]$BIC < minbic) {
        minbic <- mtfaress[[i]][[j]]$BIC
        mini <- i
        minj <- j
      }
    }
  }
}

mc <- Mclust(samples, G=ncomponents)

mtfa <- mtfaress[[mini]][[minj]]

if (log_info) {
  cat("mtfa BICs:\n")
  rownames(BICs) <- ncomponents
  colnames(BICs) <- nfactors
  print(BICs)
  
  cat("\nmclust BICs:\n")
  print(mc$BIC)
}

unlink(tmpfile)
file <- h5file(tmpfile)
gr <- createGroup(file, "fitted_proposal")

if (-max(mc$BIC, na.rm=T) < minbic) {
  if (log_info) {
    cat("\nUsing mclust fit\n\n")
  }
  gr[["weights"]] <- mc$parameters$pro
  for (k in 1:length(mc$parameters$pro)) {
    mean <- mc$parameters$mean[,k]
    cov <- mc$parameters$variance$sigma[,,k]
    gr[[paste("mean", k, sep="")]] <- mean
    gr[[paste("covariance", k, sep="")]] <- cov
  }
} else {
  if (log_info) {
    cat("\nUsing mtfa fit\n\n")
  }
  gr[["weights"]] <- mtfa$pivec[1,]
  for (k in 1:mtfa$g) {
    mean <- mtfa$mu[,k]
    cov <- mtfa$B[,,k] %*% t(mtfa$B[,,k]) + mtfa$D
    gr[[paste("mean", k, sep="")]] <- mean
    gr[[paste("covariance", k, sep="")]] <- cov
  }
}
h5close(file)

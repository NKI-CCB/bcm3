bcm3.init.cpp <- function(bcm3, clparam = "", threads = NA) {
  if (.Platform$OS.type == "windows") {
    # Need to go into the Release directory to resolve dependency dlls
    cwd <- getwd()
    setwd(paste(Sys.getenv("BCM3_ROOT"), "/bin/Release/", sep=""))
    dyn.load("bcmrbridge.dll")
    setwd(cwd)
    bcm3$.cppdllfn <- paste(Sys.getenv("BCM3_ROOT"), "/bin/Release/bcmrbridge.dll", sep="")
  } else {
    bcm3$.cppdllfn <- paste(Sys.getenv("BCM3_ROOT"), "/bin/bcmrbridge.so", sep="")
    dyn.load(bcm3$.cppdllfn)
  }
  
  if (is.na(threads)) {
    threads <- as.integer(-1)
  }
  
  res <- .C("bcm3_rbridge_init", "", bcm3$base_folder, bcm3$prior$file_name, bcm3$likelihood$file_name, clparam, as.integer(threads), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[7]] != 0) {
    dyn.unload(bcm3$.cppdllfn)
    stop(paste("BCM3 C++ bridge init failed:", res[[7]]))
  }
  
  bcm3$.cpp <- res[[1]]
  return(bcm3)
}

bcm3.reinit.cpp <- function(bcm3, clparam = "", threads = NA) {
  res <- .C("bcm3_rbridge_cleanup", bcm3$.cpp, as.integer(0), PACKAGE="bcmrbridge")
  if (res[[2]] != 0) {
    stop(paste("BCM3 C++ bridge cleanup failed:", res[[2]]))
  }
  
  if (is.na(threads)) {
    threads <- as.integer(-1)
  }
  
  res <- .C("bcm3_rbridge_init", "", bcm3$base_folder, bcm3$prior$file_name, bcm3$likelihood$file_name, clparam, as.integer(threads), as.integer(0), PACKAGE="bcmrbridge")
  if (res[[7]] != 0) {
    dyn.unload(bcm3$.cppdllfn)
    stop(paste("BCM3 C++ bridge init failed:", res[[7]]))
  }
  
  bcm3$.cpp <- res[[1]]
  return(bcm3)
}

bcm3.release.cpp <- function(bcm3) {
  dyn.unload(bcm3$.cppdllfn)
  bcm3$.cpp <- ""
  bcm3$.cppdllfn <- ""
  return(bcm3)
}

source(paste(Sys.getenv("BCM3_ROOT"), "/R/evaluate_cellpop.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/evaluate_dynamicISA.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/evaluate_fISA.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/evaluate_popPK.r", sep=""))

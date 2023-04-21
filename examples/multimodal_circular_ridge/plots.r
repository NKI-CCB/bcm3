source(paste(Sys.getenv("BCM3_ROOT"), "/R/load.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/plots_functions.r", sep=""))

model <- bcm3.load.results(".", "output_t1_n1_e1_clustered")

sample_ix <- (dim(model$posterior$samples)[3]/2+1):(dim(model$posterior$samples)[3])
temperature_ix <- dim(model$posterior$samples)[2]
plot(t(model$posterior$samples[c(1,2),temperature_ix,sample_ix]), xlim=c(-6,6), ylim=c(-2.5,2.5), col=rgb(0,0,0,0.05), pch=19)

marginal_likelihood(model)

plot_trace(model, var_ix=1)
plot_trace(model, var_ix=2)

adaptation_iter <- 1
i <- 1
j <- 2
x <- model$sampler_adaptation[[adaptation_iter]]$clustering_input_samples[i,] * model$sampler_adaptation[[adaptation_iter]]$clustering_input_sample_scaling[i]
y <- model$sampler_adaptation[[adaptation_iter]]$clustering_input_samples[j,] * model$sampler_adaptation[[adaptation_iter]]$clustering_input_sample_scaling[j]
plot(x, y, col=model$sampler_adaptation[[adaptation_iter]]$assignment+1, pch=20)

nclusters <- length(unique(model$sampler_adaptation[[adaptation_iter]]$assignment))

library(ellipse)
for (clusti in 1:nclusters) {
  clust_sample_ix <- which(model$sampler_adaptation[[adaptation_iter]]$assignment == clusti-1)
  clustx <- x[clust_sample_ix]
  clusty <- y[clust_sample_ix]
  name <- paste("cluster", clusti-1, "_covariance", sep="")
  cov <- model$sampler_adaptation[[adaptation_iter]][[name]][c(i,j),c(i,j)]
  ell <- ellipse(cov, centre = c(mean(clustx), mean(clusty)), level=0.6, draw=F)
  lines(ell, col=clusti)
}

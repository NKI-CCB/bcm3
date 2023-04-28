source(paste(Sys.getenv("BCM3_ROOT"), "/R/load.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/plots_functions.r", sep=""))

model <- bcm3.load.results(".", "output_t4_n5_e1_clust")
model <- bcm3.load.results(".", "output_t4_n5_e1_gmm")

sample_ix <- (dim(model$posterior$samples)[3]/2+1):(dim(model$posterior$samples)[3])
temperature_ix <- dim(model$posterior$samples)[2]
plot(t(model$posterior$samples[c(1,2),temperature_ix,sample_ix]), xlim=c(-10,10), ylim=c(-10,10), col=rgb(0,0,0,0.05), pch=19)

#plot_variable_distribution(model, 1)
plot_variable_distribution(model, 2, ylim=c(0, 0.3))

sum(model$posterior$samples[2,temperature_ix,sample_ix] > 0) / 5000

# 
library(mvtnorm)
x <- rmvnorm(5000, c(0,0), rbind(c(2,-1),c(-1,1)))
plot(x, xlim=c(-5,5), ylim=c(-5,5))
# 


library(ellipse)

i <- 1
j <- 2
adaptation_iter <- 1
x <- model$sampler_adaptation[[adaptation_iter]]$clustering_input_samples[i,] * model$sampler_adaptation[[adaptation_iter]]$clustering_input_sample_scaling[i]
y <- model$sampler_adaptation[[adaptation_iter]]$clustering_input_samples[j,] * model$sampler_adaptation[[adaptation_iter]]$clustering_input_sample_scaling[j]
plot(x, y, col=model$sampler_adaptation[[adaptation_iter]]$assignment+1, pch=20)

nclusters <- length(unique(model$sampler_adaptation[[adaptation_iter]]$assignment))

for (clusti in 1:nclusters) {
  clust_sample_ix <- which(model$sampler_adaptation[[adaptation_iter]]$assignment == clusti-1)
  clustx <- x[clust_sample_ix]
  clusty <- y[clust_sample_ix]
  name <- paste("cluster", clusti-1, "_covariance", sep="")
  cov <- model$sampler_adaptation[[adaptation_iter]][[name]][c(i,j),c(i,j)]
  ell <- ellipse(cov, centre = c(mean(clustx), mean(clusty)), level=0.6, draw=F)
  lines(ell, col=clusti)
}



i <- 1
j <- 2
plot(t(model$posterior$samples[c(i,j),temperature_ix,sample_ix]), pch='.')
for (clusti in 1:5) {
  mean <- model$sampler_adaptation$info0[[paste("cluster", clusti-1, "_mean", sep="")]][c(i,j)]
  cov <- model$sampler_adaptation$info0[[paste("cluster", clusti-1, "_covariance", sep="")]][c(i,j),c(i,j)]
  ell <- ellipse(cov, centre = mean, level=0.6, draw=F)
  lines(ell, col=clusti, lwd=2)
}


x <- rmvnorm(5000, c(0,0), rbind(c(2,-1),c(-1,1)))
plot(x, xlim=c(-5,5), ylim=c(-5,5))

Sigma <- cov(x)
eig <- eigen(Sigma)
shrunk_eigval <- eig$values
recov <- eig$vectors %*% diag(shrunk_eigval) %*% t(eig$vectors)

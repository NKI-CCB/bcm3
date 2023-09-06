source(paste(Sys.getenv("BCM3_ROOT"), "/R/load.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/plots_functions.r", sep=""))

model <- bcm3.load.results(".", "output_t6_n5_e1")

sample_ix <- (dim(model$posterior$samples)[3]/2+1):(dim(model$posterior$samples)[3])
temperature_ix <- dim(model$posterior$samples)[2]

plot(t(model$posterior$samples[,temperature_ix,sample_ix]))

plot_trace(model, var_ix=1)
plot_trace(model, var_ix=2)

library(ellipse)
library(pals)

clustcols <- brewer.set1(13)
adaptation_iter <- 1

adapt_sample_ix <- 1:2000
i <- 1
j <- 2
x <- model$posterior$samples[i,temperature_ix,adapt_sample_ix]
y <- model$posterior$samples[j,temperature_ix,adapt_sample_ix]

adapt_group <- model$sampler_adaptation[["adapt1"]][["block1"]]
nclusters <- (length(adapt_group)-2)/3

plot(x, y, pch='.')
for (clusti in 1:nclusters) {
  mean <- adapt_group[[paste("cluster", clusti-1, "_mean", sep="")]][c(i,j)]
  cov <- adapt_group[[paste("cluster", clusti-1, "_covariance", sep="")]][c(i,j),c(i,j)]
  ell <- ellipse(cov, centre = mean, level=0.6, draw=F)
  lines(ell, col=clustcols[clusti], lwd=2)
}


source(paste(Sys.getenv("BCM3_ROOT"), "/R/load.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/plots_functions.r", sep=""))

model1 <- bcm3.load.results(".", "output_t2_n5_e1_gmm")
model2 <- bcm3.load.results(".", "output_t2_n5_e1_globalcov")
model3 <- bcm3.load.results(".", "output_t2_n5_e1_gmm_noscaleadapt")

model <- model1
sample_ix <- (dim(model$posterior$samples)[3]/2+1):(dim(model$posterior$samples)[3])
temperature_ix <- dim(model$posterior$samples)[2]

plot(t(model1$posterior$samples[c(1,2),temperature_ix,sample_ix]), xlim=c(-10,10), ylim=c(-10,10), col=rgb(0,0,0,0.05), pch=19)
plot(t(model2$posterior$samples[c(1,2),temperature_ix,sample_ix]), xlim=c(-10,10), ylim=c(-10,10), col=rgb(0,0,0,0.05), pch=19)
plot(t(model3$posterior$samples[c(1,2),temperature_ix,sample_ix]), xlim=c(-10,10), ylim=c(-10,10), col=rgb(0,0,0,0.05), pch=19)

#plot_variable_distribution(model, 1)
#plot_variable_distribution(model, 2, ylim=c(0, 0.3))

sum(model1$posterior$samples[2,temperature_ix,sample_ix] > 0) / length(sample_ix)
sum(model2$posterior$samples[2,temperature_ix,sample_ix] > 0) / length(sample_ix)
sum(model3$posterior$samples[2,temperature_ix,sample_ix] > 0) / length(sample_ix)


# plot_variable_distribution(model1, 2)
# plot_variable_distribution(model2, 2)
# 
plot(model1$posterior$samples[1,temperature_ix,sample_ix[6601:7000]])
plot(model2$posterior$samples[1,temperature_ix,sample_ix[6601:7000]])
plot(model3$posterior$samples[1,temperature_ix,sample_ix[6601:7000]])


var_ix <- 1
lag_max <- 20

#acf(model1$posterior$samples[2,temperature_ix,sample_ix], lag.max = 50)
acf(t(model1$posterior$samples[2,,sample_ix]), lag.max = lag_max)
acf(t(model2$posterior$samples[2,,sample_ix]), lag.max = lag_max)

# acfs <- list()
# acfs[[1]] <- acf(model1$posterior$samples[var_ix,temperature_ix,sample_ix], lag.max = lag_max, plot=F)
# acfs[[2]] <- acf(model2$posterior$samples[var_ix,temperature_ix,sample_ix], lag.max = lag_max, plot=F)
# acfs[[3]] <- acf(model3$posterior$samples[var_ix,temperature_ix,sample_ix], lag.max = lag_max, plot=F)

acfs <- list()
acfs[[1]] <- acf(t(model1$posterior$samples[var_ix,,sample_ix]), lag.max = lag_max, plot=F)
acfs[[2]] <- acf(t(model2$posterior$samples[var_ix,,sample_ix]), lag.max = lag_max, plot=F)
acfs[[3]] <- acf(t(model3$posterior$samples[var_ix,,sample_ix]), lag.max = lag_max, plot=F)

boundary <- qnorm((1 + 0.95)/2)/sqrt(length(sample_ix))

library(pals)
cols <- brewer.set1(3)
plot(0, type='n', xlim=c(0,lag_max*1.1), ylim=c(-0.1,1))
lines(acfs[[1]]$lag[,1,temperature_ix], acfs[[1]]$acf[,temperature_ix,temperature_ix], col=cols[1], lwd=2)
lines(acfs[[2]]$lag[,1,temperature_ix], acfs[[2]]$acf[,temperature_ix,temperature_ix], col=cols[2], lwd=2)
lines(acfs[[3]]$lag[,1,temperature_ix], acfs[[3]]$acf[,temperature_ix,temperature_ix], col=cols[3], lwd=2)
legend("topright", legend=c("GMM", "Global covariance", "No adaptation"), col=cols, lwd=2)
abline(h=0, lty=1)
abline(h=c(boundary,-boundary), lty=2)


distances <- list()
diff <- model1$posterior$samples[,temperature_ix,sample_ix[-1]] - model1$posterior$samples[,temperature_ix,sample_ix[-length(sample_ix)]]
distances[[1]] <- sqrt(apply(diff*diff, 2, sum))
plot(distances[[1]])
diff <- model2$posterior$samples[,temperature_ix,sample_ix[-1]] - model2$posterior$samples[,temperature_ix,sample_ix[-length(sample_ix)]]
distances[[2]] <- sqrt(apply(diff*diff, 2, sum))
diff <- model3$posterior$samples[,temperature_ix,sample_ix[-1]] - model3$posterior$samples[,temperature_ix,sample_ix[-length(sample_ix)]]
distances[[3]] <- sqrt(apply(diff*diff, 2, sum))
points(distances[[2]], col="blue")
#mean(distances[distances < 10])

plot(density(distances[[1]]), ylim=c(0,1))
lines(density(distances[[2]]), col="blue")
lines(density(distances[[3]]), col="green")

h5 <- h5file("gmmtest.nc")
curpos <- h5[["proposal"]][["current_position"]][]
samples <- t(h5[["proposal"]][["samples"]][,])
mhratios <- h5[["proposal"]][["mh_ratios"]][]
h5close(h5)

points(samples, pch=20, col=rgb(0,0,1,0.2))
points(t(curpos), pch=19, col='red', xlim=c(-5,10), ylim=c(-5,10))
plot(mhratios)


library(mclust)
mc <- Mclust(t(model1$posterior$samples[,temperature_ix,sample_ix]), modelNames="VVV", G=1:3)

mc$parameters$mean
mc$parameters$variance$sigma[,,1]
mc$parameters$variance$sigma[,,2]

eig(mc$parameters$variance$sigma[,,1])
eig(mc$parameters$variance$sigma[,,2])

# 
library(mvtnorm)
x <- rmvnorm(5000, c(0,0), rbind(c(2,-0.5),c(-0.5,1)))
plot(x, xlim=c(-5,5), ylim=c(-5,5))
# 

library(ellipse)
model <- model1
if (!is.null(model$sampler_adaptation)) {
  clustcols <- brewer.set1(13)
  
  for (iter in length(model$sampler_adaptation)) {
    name <- names(model$sampler_adaptation)[iter]
    group <- model$sampler_adaptation[[name]]

    i <- 1
    j <- 2
    
        x <- model$posterior$samples[i,temperature_ix,sample_ix]
        y <- model$posterior$samples[j,temperature_ix,sample_ix]
        nclusters <- (length(model$sampler_adaptation[[name]])-2)/3
        
        
        plot(x, y, pch='.')
        for (clusti in 1:nclusters) {
          mean <- model$sampler_adaptation[[name]][[paste("cluster", clusti-1, "_mean", sep="")]][c(i,j)]
          cov <- model$sampler_adaptation[[name]][[paste("cluster", clusti-1, "_covariance", sep="")]][c(i,j),c(i,j)]
          ell <- ellipse(cov, centre = mean, level=0.6, draw=F)
          lines(ell, col=clustcols[clusti], lwd=2)
        }
  }
}

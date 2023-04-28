source(paste(Sys.getenv("BCM3_ROOT"), "/R/load.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/plots_functions.r", sep=""))

model1 <- bcm3.load.results(".", "output_t4_n5_e1")
model2 <- bcm3.load.results(".", "output_t4_n5_e1_globalcov")

sample_ix <- (dim(model$posterior$samples)[3]/2+1):(dim(model$posterior$samples)[3])
temperature_ix <- dim(model$posterior$samples)[2]
plot(t(model1$posterior$samples[c(1,2),temperature_ix,sample_ix]), xlim=c(-10,10), ylim=c(-10,10), col=rgb(0,0,0,0.05), pch=19)
plot(t(model2$posterior$samples[c(1,2),temperature_ix,sample_ix]), xlim=c(-10,10), ylim=c(-10,10), col=rgb(0,0,0,0.05), pch=19)

#plot_variable_distribution(model, 1)
#plot_variable_distribution(model, 2, ylim=c(0, 0.3))

sum(model$posterior$samples[2,temperature_ix,sample_ix] > 0) / 5000


# plot_variable_distribution(model1, 2)
# plot_variable_distribution(model2, 2)
# 
# plot(model1$posterior$samples[3,temperature_ix,sample_ix], pch='.')
# plot(model2$posterior$samples[3,temperature_ix,sample_ix], pch='.')
# 

#acf(model1$posterior$samples[2,temperature_ix,sample_ix], lag.max = 50)
acf(t(model1$posterior$samples[2,,sample_ix]), lag.max = 50)
acf(t(model2$posterior$samples[2,,sample_ix]), lag.max = 50)

var_ix <- 2
lag_max <- 50
acfs <- list()
acfs[[1]] <- acf(t(model1$posterior$samples[var_ix,,sample_ix]), lag.max = lag_max, plot=F)
acfs[[2]] <- acf(t(model2$posterior$samples[var_ix,,sample_ix]), lag.max = lag_max, plot=F)
#acfs[[3]] <- acf(t(model3$posterior$samples[3,,sample_ix]), lag.max = lag_max, plot=F)

boundary <- qnorm((1 + 0.95)/2)/sqrt(length(sample_ix))

library(pals)
cols <- brewer.set1(3)
plot(0, type='n', xlim=c(0,lag_max*1.1), ylim=c(-0.1,1))
lines(acfs[[1]]$lag[,temperature_ix,temperature_ix], acfs[[1]]$acf[,temperature_ix,temperature_ix], col=cols[1], lwd=2)
lines(acfs[[2]]$lag[,temperature_ix,temperature_ix], acfs[[2]]$acf[,temperature_ix,temperature_ix], col=cols[2], lwd=2)
#lines(acfs[[3]]$lag[,temperature_ix,temperature_ix], acfs[[3]]$acf[,temperature_ix,temperature_ix], col=cols[3], lwd=2)
legend("topright", legend=c("GMM", "Global covariance"), col=cols, lwd=2)
abline(h=0, lty=1)
abline(h=c(boundary,-boundary), lty=2)




# 
library(mvtnorm)
x <- rmvnorm(5000, c(0,0), rbind(c(2,-1),c(-1,1)))
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

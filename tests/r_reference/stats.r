library(mvtnorm)

means <- list()
means[[1]] <- c(-1,-1)
means[[2]] <- c(2, 3);

covariances <- list()
covariances[[1]] <- rbind(c(2, 1), c(1, 1))
covariances[[2]] <- rbind(c(1, -0.9), c(-0.9, 1))

weights <- c(0.25, 0.75)

n <- 5000

samples <- rbind(rmvnorm(n * weights[1], means[[1]], covariances[[1]]),
                 rmvnorm(n * weights[2], means[[2]], covariances[[2]]))

plot(samples)

#x <- c(2, 2)
x <- c(2, 0)
p <- weights[1] * dmvnorm(x, means[[1]], covariances[[1]]) + weights[2] * dmvnorm(x, means[[2]], covariances[[2]])
sprintf("%.16f",log(p))

p1 <- weights[1] * dmvnorm(x, means[[1]], covariances[[1]], log=F)
p2 <- weights[2] * dmvnorm(x, means[[2]], covariances[[2]], log=F)

sprintf("%.16f", p1 / (p1 + p2))
sprintf("%.16f", p2 / (p1 + p2))

exp(p1 - logsum)
exp(p2 - logsum)

library(crch)
library(sm)
library(extraDistr, warn.conflicts = F)
source(paste(Sys.getenv("BCM3_ROOT"), "/R/stats.r", sep=""))

.prior_color <- "#4053d3"
.posterior_color <- "#ddb310"
.posterior_predictive_color <- "#00b25d"
.posterior_predictive_color_alpha <- "#00b25d64"

# Plot prior and posterior distribution for a given variable
# You can specify either a variable name or a variable index
plot_variable_distribution <- function(model, var_ix=NULL, var_name=NULL, temperature_ix=NULL, sample_ix=NULL, xlab="", ylim=NULL, plot=T, adjust=1)
{
  if (!is.null(var_name)) {
    var_ix <- model_get_var_ix(model, var_name)
  }
  if (!is.null(var_ix)) {
    if (is.null(temperature_ix)) {
      temperature_ix <- dim(model$posterior$samples)[2]
    }
    if (is.null(sample_ix)) {
      sample_ix <- 1:dim(model$posterior$samples)[3]
    }
    
    varattrs <- model$prior$variable_attrs[[var_ix]]
    res <- plot_variable_distribution_impl(model$posterior$samples[var_ix, temperature_ix, sample_ix], 
                                           model$posterior$weights[temperature_ix, sample_ix],
                                           varattrs, xlab, ylim, plot=plot, adjust=adjust)
    if(!is.null(res)) {
      return(res)
    }
  } else {
    stop("Either a variable index or a variable name has to be specified")
  }
}

plot_all_densities <- function(model, sample_ix=NULL, imagetype="png")
{
  if (is.null(sample_ix)) {
    sample_ix <- (dim(model$posterior$samples)[3]/2+1):dim(model$posterior$samples)[3]
  }
  if (imagetype == "pdf") {
    pdf_tile(paste(model$output_folder, "/densities.pdf", sep=""), width=12, height=10, model$nvar)
  } else if (imagetype == "png") {
    png_tile(paste(model$output_folder, "/densities.png", sep=""), width=1600, height=1200, model$nvar)
  } else {
    stop("Unknown image type; should be either png or pdf")
  }
  par(mar=c(3,2,2,1))
  for (i in 1:model$nvar) {
    plot_variable_distribution(model, var_ix=i, sample_ix=sample_ix)
  }
  par(mfrow=c(1,1))
  res <- dev.off()
}

plot_all_traces <- function(model, burnin_cutoff=NULL, imagetype="png")
{
  if (is.null(burnin_cutoff)) {
    burnin_cutoff <- dim(model$posterior$samples)[3]/2+1
  }
  if (imagetype == "pdf") {
    pdf_tile(paste(model$output_folder, "/traces.pdf", sep=""), width=12, height=10, model$nvar)
  } else if (imagetype == "png") {
    png_tile(paste(model$output_folder, "/traces.png", sep=""), width=1600, height=1200, model$nvar)
  } else {
    stop("Unknown image type; should be either png or pdf")
  }
  par(mar=c(2,2,2,2))
  for (i in 1:model$nvar) {
    plot_trace(model, var_ix=i)
    abline(v=burnin_cutoff-0.5, col='grey', lty=2)
  }
  par(mfrow=c(1,1))
  res <- dev.off()
}

plot_variable_prior <- function(model, var_ix=NULL, var_name=NULL, xlab="", ylim=NULL, plot=T)
{
  if (!is.null(var_name)) {
    var_ix <- model_get_var_ix(model, var_name)
  }
  if (!is.null(var_ix)) {
    varattrs <- model$prior$variable_attrs[[var_ix]]
    res <- plot_variable_prior_impl(varattrs, xlab, ylim, plot=plot)
    return(res)
  } else {
    stop("Either a variable index or a variable name has to be specified")
  }
}

plot_trace <- function(model, var_ix=NULL, var_name=NULL, temperature_ix=NULL)
{
  if (!is.null(var_name)) {
    var_ix <- model_get_var_ix(model, var_name)
  }
  if (!is.null(var_ix)) {
    if (is.null(temperature_ix)) {
      temperature_ix <- dim(model$posterior$samples)[2]
    }
    varattrs <- model$prior$variable_attrs[[var_ix]]
    plot(model$posterior$samples[var_ix, temperature_ix, ], main=varattrs["name"], pch='.')
  } else {
    stop("Either a variable index or a variable name has to be specified")
  }
}

# Posterior predictive distribution as a barplot
ppd_barplot <- function(model, variable_samples, data, labels, sd_var_ix=NULL, sd_incr_var_ix=NULL, sd_samples=NULL, error_model="t", bounds=c(0.05, 0.95), ylim=NULL, ppdsamples = 20, pointsize = 1.5, ...) {
  bound_color <- rgb(43,131,186, max=255)
  bound_color_alpha <- rgb(43,131,186, 0.4*255, max=255)
  
  if (is.null(ncol(data))) {
    if (ncol(variable_samples) != length(data)) {
      stop("Number of columns of data and posterior samples should be the same")
    }
  } else {
    if (ncol(variable_samples) != ncol(data)) {
      stop("Number of columns of data and posterior samples should be the same")
    }
  }
  
  if (bounds[2] < bounds[1]) {
    bounds <- rev(bounds)
  }
  
  if (is.null(ylim)) {
    sample_min <- min(apply(variable_samples, 2, quantile, probs=bounds[1], na.rm=T))
    sample_max <- max(apply(variable_samples, 2, quantile, probs=bounds[2], na.rm=T))
    ylim <- c(min(data, sample_min, na.rm=T), max(data, sample_max, na.rm=T))
    ylim <- c(ylim[1] - 0.05 * (ylim[2]-ylim[1]), ylim[2] + 0.05 * (ylim[2]-ylim[1]))
  }
  
  xpos <- barplot(rep(0,ncol(variable_samples)), names.arg=labels, las=2, border=NA, ylim=ylim, ...)
  for (i in 1:ncol(variable_samples)) {
    barwidth <- 0.6
    lx <- xpos[i] - barwidth / 2
    ux <- xpos[i] + barwidth / 2
    
    # Posterior
    tboundscale <- 0
    if (!is.null(sd_samples)) {
      tboundscale <- sd_samples
    } else if (!is.null(sd_var_ix)) {
      tboundscale <- model$posterior$samples[sd_var_ix,dim(model$posterior$samples)[2],]
    }
    if (!is.null(sd_incr_var_ix)) {
      tboundscale <- model$posterior$samples[sd_var_ix,dim(model$posterior$samples)[2],] + variable_samples[,i] * model$posterior$samples[sd_incr_var_ix,sd_var_ix,dim(model$posterior$samples)[2],]
    }
    
    predy <- NULL
    if (error_model == "truncated_t") {
      predy <- rtt(ppdsamples * nrow(variable_samples), variable_samples[,i], tboundscale, 3, 0, 1)
    } else if (error_model == "t") {
      predy <- rtt(ppdsamples * nrow(variable_samples), variable_samples[,i], tboundscale, 3, -Inf, Inf)
    } else if (error_model == "normal") {
      predy <- rnorm(ppdsamples * nrow(variable_samples), variable_samples[,i], tboundscale)
    } else if (error_model == "truncated_normal") {
      predy <- rtnorm(ppdsamples * nrow(variable_samples), variable_samples[,i], tboundscale, 0, 1)
    } else {
      stop("Unknown error model")
    }
    ly <- quantile(predy, probs=bounds[1], na.rm=T)
    uy <- quantile(predy, probs=bounds[2], na.rm=T)
    iqrly <- quantile(predy, probs=0.25, na.rm=T)
    iqruy <- quantile(predy, probs=0.75, na.rm=T)
    polygon(c(lx,lx,ux,ux), c(max(ly, ylim[1]),min(uy,ylim[2]),min(uy, ylim[2]),max(ly, ylim[1])), col=.posterior_predictive_color_alpha, border=F)
    lines(c(lx,ux), c(ly,ly), col=.posterior_predictive_color)
    lines(c(lx,ux), c(uy,uy), col=.posterior_predictive_color)

    # Data
    if (is.null(nrow(data))) {
      points(xpos[i], data[i], pch=20, cex = pointsize)
    } else {
      points(rep(xpos[i],nrow(data)), data[,i], pch=20, cex = pointsize)
    }
  }
  return(xpos)
}

# Posterior predictive distribution as a lineplot
ppd_lineplot <- function(x.data, y.data, x.samples, y.samples, bounds=c(0.05, 0.95), xlim=NULL, ylim=NULL, data.color='black', data_plot_type = 'p', median_line = T, ...)
{

  if (bounds[2] < bounds[1]) {
    bounds <- rev(bounds)
  }
  
  ly <- apply(y.samples, 2, quantile, bounds[1], na.rm=T)
  my <- apply(y.samples, 2, quantile, 0.5, na.rm=T)
  uy <- apply(y.samples, 2, quantile, bounds[2], na.rm=T)
  if (all(is.na(y.data)) && all(is.na(my))) {
    # No data and no modeled values - nothing meaningful to plot
    plot(1, type="n")
    return()
  }
  
  if (is.null(ylim)) {
    sample_min <- min(ly, na.rm=T)
    sample_max <- max(uy, na.rm=T)
    ylim <- c(min(y.data, sample_min, na.rm=T), max(y.data, sample_max, na.rm=T))
    ylim <- c(ylim[1] - 0.05 * (ylim[2]-ylim[1]), ylim[2] + 0.05 * (ylim[2]-ylim[1]))
  }
  if (is.null(xlim)) {
    xlim <- range(c(x.data, x.samples), na.rm=T)
  }
  
  poly_point_ix <- which(!is.na(my))
  
  plot(1, type="n", xlim=xlim, ylim=ylim, col=.posterior_predictive_color, ...)
  if (median_line) {
    lines(x.samples[poly_point_ix], my[poly_point_ix], lwd=2)
  }
  polygon(c(x.samples[poly_point_ix], rev(x.samples[poly_point_ix])), c(ly[poly_point_ix],rev(uy[poly_point_ix])), col=.posterior_predictive_color_alpha, border=F)
  lines(x.samples, ly, col=.posterior_predictive_color)
  lines(x.samples, uy, col=.posterior_predictive_color)
  if (is.null(nrow(y.data))) {
    points(x.data, y.data, type=data_plot_type, pch=19, col=data.color)
  } else {
    for (i in 1:nrow(y.data)) {
      points(x.data, y.data[i,], type=data_plot_type, pch=19, col=data.color)
    }
  }
}

plot_bivariate_variable_distribution <- function(model, parix1, parix2, hscale=1, probrange=NULL, other_samples=NULL)
{
  library(ComplexHeatmap)
  library(ks)
  library(pals)
  library(MASS)
  library(circlize)
  library(mvtnorm)
  
  xrange <- prior_bounds(model, parix1, q=c(0,1))
  yrange <- prior_bounds(model, parix2, q=c(0,1))
  
  if (is.null(other_samples)) {
    samples1 <- model$posterior$samples[,parix1]
    samples2 <- model$posterior$samples[,parix2]
  } else {
    samples1 <- other_samples[,1]
    samples2 <- other_samples[,2]
  }
  
  x <- c(samples1, xrange[1]+(xrange[1]-samples1), xrange[2] + (xrange[2]-samples1),
         samples1, xrange[1]+(xrange[1]-samples1), xrange[2] + (xrange[2]-samples1),
         samples1, xrange[1]+(xrange[1]-samples1), xrange[2] + (xrange[2]-samples1))
  y <- c(samples2, samples2, samples2,
         yrange[1]+(yrange[1]-samples2), yrange[1]+(yrange[1]-samples2),yrange[1]+(yrange[1]-samples2),
         yrange[2]+(yrange[2]-samples2), yrange[2]+(yrange[2]-samples2),yrange[2]+(yrange[2]-samples2))
  
  bw <- Hpi(cbind(x, y))
  z <- matrix(NA, 20, 20)
  for (i in 1:20) {
    for (j in 1:20) {
      x1 <- xrange[1] + (xrange[2] - xrange[1]) * (i-1) / (19)
      x2 <- yrange[1] + (yrange[2] - yrange[1]) * (j-1) / (19)
      p <- mvtnorm::dmvnorm(cbind(x, y), t(c(x1, x2)), bw*hscale) / length(x)
      z[i,j] <- sum(p)
    }
  }
  
  if (is.null(probrange)) {
    probrange <- c(0, max(z))
  }
  
  Heatmap(t(z), row_order=20:1, cluster_rows = F, cluster_columns = F,
          col=colorRamp2(seq(probrange[1], probrange[2], len=20), viridis(20)),
          column_title=model$variables[parix1],
          column_title_side="bottom",
          row_title=model$variables[parix2],
          heatmap_legend_param=list(title="Probability\ndensity\n"))
}



#
# Internal functions
#

load_binary <- function(filename, numvar) {
  file_size <- file.info(filename)$size
  numsamples <- file_size / (4 * numvar)
  
  infile <- file(filename, "rb");
  data <- matrix(NA, numsamples, numvar);
  for (sample in 1:numsamples) {
    d <- readBin(infile, double(), n=numvar, size=4);
    data[sample,] <- d;
  }
  close(infile);
  return(data);
}

plot_variable_distribution_impl <- function(samples, weights, varattrs, xlab="", ylim=NULL, plot=T, adjust=1, lwd=2)
{
  name = varattrs["name"];
  
  minx <- min(samples);
  maxx <- max(samples);
  
  lbound <- NA
  ubound <- NA
  
  priortype <- varattrs["type"]
  if (is.na(priortype) || priortype == "regular") {
    distribution <- varattrs["distribution"]
    if (distribution == "normal") {
      mu = as.numeric(varattrs["mu"]);
      sigma = as.numeric(varattrs["sigma"]);
      minx <- min(minx, qnorm(0.01, mean=mu, sd=sigma));
      maxx <- max(maxx, qnorm(0.99, mean=mu, sd=sigma));
      priorx <- seq(minx, maxx, length.out=200);
      priory <- dnorm(priorx, mean=mu, sd=sigma);
    } else if (distribution == "gamma") {
      k = as.numeric(varattrs["k"]);
      theta = as.numeric(varattrs["theta"]);
      minx <- 0;#min(minx, qgamma(0.01, shape=k, scale=theta));
      maxx <- max(maxx, qgamma(0.99, shape=k, scale=theta));
      priorx <- seq(minx, maxx, length.out=200);
      priory <- dgamma(priorx, shape=k, scale=theta);
      lbound <- 0;
    } else if (distribution == "uniform") {
      a = as.numeric(varattrs["lower"]);
      b = as.numeric(varattrs["upper"]);
      minx <- a;
      maxx <- b;
      priorx <- seq(minx, maxx, length.out=200);
      priory <- dunif(priorx, min=a, max=b);
      lbound <- a;
      ubound <- b;
    } else if (distribution == "half_cauchy") {
      scale = as.numeric(varattrs["scale"]);
      minx <- 0;
      maxx <- max(maxx, qcauchy(0.95, scale=scale));
      priorx <- seq(minx, maxx, length.out=200);
      priory <- 2.0 * dcauchy(priorx, scale=scale);
      lbound <- 0;
    } else if (distribution == "beta") {
      a = as.numeric(varattrs["a"]);
      b = as.numeric(varattrs["b"]);
      minx <- 1e-6;
      maxx <- 1-1e-6;
      priorx <- seq(minx, maxx, length.out=200);
      priory <- dbeta(priorx, a, b);
      lbound <- 0;
      ubound <- 1;
    } else if (distribution == "exponential") {
      rate = as.numeric(varattrs["lambda"]);
      minx <- 0;
      maxx <- max(maxx, qexp(0.99, rate=rate));
      priorx <- seq(minx, maxx, length.out=200);
      priory <- dexp(priorx, rate=rate);
      lbound <- 0;
    } else if (distribution == "beta_prime") {
      a = as.numeric(varattrs["a"]);
      b = as.numeric(varattrs["b"]);
      scale = as.numeric(varattrs["scale"]);
      minx <- 0;
      maxx <- max(quantile(samples, probs=0.995), qbetapr(0.99, shape1=a, shape2=b, scale=scale));
      priorx <- seq(minx, maxx, length.out=200);
      priory <- dbetapr(priorx, shape1=a, shape2=b, scale=scale)
      lbound <- 0;
    } else if (distribution == "dirichlet") {
      # Todo
      minx <- 1e-6;
      maxx <- 1-1e-6;
      priorx <- seq(minx, maxx, length.out=200);
      priory <- dbeta(priorx, 1, 9);
      lbound <- 0;
      ubound <- 1;
    } else if (distribution == "exponential_mix") {
      rate1 = as.numeric(varattrs["lambda"]);
      rate2 = as.numeric(varattrs["lambda2"]);
      mix = as.numeric(varattrs["mix"]);
      minx <- 0
      maxx <- max(qexp(0.99, rate1), qexp(0.99, rate2))
      priorx <- seq(minx, maxx, length.out=200);
      priory <- mix * dexp(priorx, rate = rate1) + (1-mix) * dexp(priorx, rate = rate2)
      lbound <- 0;
    }
  }
  
  # Select bandwidth
  bw <- h.select(samples, weights=weights, method="cv")
  
  # If necessary, reflect samples around bounds and add the reflection to the sample, to get a more reasonable KDE near the bounds
  factor <- 1
  density_samples <- samples
  density_weights <- weights
  if (!is.na(lbound)) {
    reflected <- lbound - (samples - lbound)
    density_samples <- c(density_samples, reflected)
    density_weights <- c(density_weights, weights)
    factor <- factor + 1
  }
  if (!is.na(ubound)) {
    reflected <- ubound + (ubound - samples)
    density_samples <- c(density_samples, reflected)
    density_weights <- c(density_weights, weights)
    factor <- factor + 1
  }
  density_weights <- density_weights / sum(density_weights)
  
  # Calculate density
  d <- density(density_samples, bw=bw*adjust, weights=density_weights, from=minx, to=maxx);
  d$y <- d$y * factor
  
  maxy <- max(d$y, priory);
  if (is.null(ylim)) {
    ylim <- c(0, maxy)
  }

  if (plot) {
    plot(priorx, priory, col=.prior_color, main=name, xlab=xlab, ylab="Probability density", type="l", xlim=c(minx, maxx), ylim=ylim, lwd=lwd);
    
    if(!is.na(lbound)) {
      lines(c(lbound - 1e6, lbound, lbound), c(0, 0, priory[1]), col=.prior_color, lwd=lwd)
      lines(c(lbound - 1e6, lbound, lbound), c(0, 0, d$y[1]), col=.posterior_color, lwd=lwd)
    }
    if(!is.na(ubound)) {
      lines(c(ubound, ubound, ubound + 1e6), c(tail(priory, 1), 0, 0), col=.prior_color, lwd=lwd)
      lines(c(ubound, ubound, ubound + 1e6), c(tail(d$y, 1), 0, 0), col=.posterior_color, lwd=lwd)
    }
    
    lines(d$x, d$y, col=.posterior_color, lwd=lwd)
  } else {
    result <- list()
    result$xlim <- c(minx, maxx)
    result$ylim <- ylim
    result$priorx <- priorx
    result$priory <- priory
    result$dx <- d$x
    result$dy <- d$y
    result$rdens <- d
    return(result)
  }
}

plot_variable_prior_impl <- function(varattrs, xlab="", ylim=NULL, plot=T)
{
  name = varattrs["name"];
  
  priortype <- varattrs["type"]
  if (is.na(priortype) || priortype == "regular") {
    distribution <- varattrs["distribution"]
    if (distribution == "normal") {
      mu = as.numeric(varattrs["mu"]);
      sigma = as.numeric(varattrs["sigma"]);
      minx <- qnorm(0.01, mean=mu, sd=sigma);
      maxx <- qnorm(0.99, mean=mu, sd=sigma);
      priorx <- seq(minx, maxx, length.out=200);
      priory <- dnorm(priorx, mean=mu, sd=sigma);
    } else if (distribution == "gamma") {
      k = as.numeric(varattrs["k"]);
      theta = as.numeric(varattrs["theta"]);
      minx <- 0;#min(minx, qgamma(0.01, shape=k, scale=theta));
      maxx <- qgamma(0.99, shape=k, scale=theta);
      priorx <- seq(minx, maxx, length.out=200);
      priory <- dgamma(priorx, shape=k, scale=theta);
      lbound <- 0;
    } else if (distribution == "uniform") {
      a = as.numeric(varattrs["lower"]);
      b = as.numeric(varattrs["upper"]);
      minx <- a;
      maxx <- b;
      priorx <- seq(minx, maxx, length.out=200);
      priory <- dunif(priorx, min=a, max=b);
      lbound <- a;
      ubound <- b;
    } else if (distribution == "half_cauchy") {
      scale = as.numeric(varattrs["scale"]);
      minx <- 0;
      maxx <- qcauchy(0.95, scale=scale);
      priorx <- seq(minx, maxx, length.out=200);
      priory <- 2.0 * dcauchy(priorx, scale=scale);
      lbound <- 0;
    } else if (distribution == "beta") {
      a = as.numeric(varattrs["a"]);
      b = as.numeric(varattrs["b"]);
      minx <- 1e-6;
      maxx <- 1-1e-6;
      priorx <- seq(minx, maxx, length.out=200);
      priory <- dbeta(priorx, a, b);
      lbound <- 0;
      ubound <- 1;
    } else if (distribution == "exponential") {
      rate = as.numeric(varattrs["lambda"]);
      minx <- 0;
      maxx <- qexp(0.99, rate=rate);
      priorx <- seq(minx, maxx, length.out=200);
      priory <- dexp(priorx, rate=rate);
      lbound <- 0;
    } else if (distribution == "beta_prime") {
      a = as.numeric(varattrs["a"]);
      b = as.numeric(varattrs["b"]);
      scale = as.numeric(varattrs["scale"]);
      minx <- 0;
      maxx <- qbetapr(0.99, shape1=a, shape2=b, scale=scale);
      priorx <- seq(minx, maxx, length.out=200);
      priory <- dbetapr(priorx, shape1=a, shape2=b, scale=scale)
      lbound <- 0;
    } else if (distribution == "dirichlet") {
      # Todo
      minx <- 1e-6;
      maxx <- 1-1e-6;
      priorx <- seq(minx, maxx, length.out=200);
      priory <- dbeta(priorx, 1, 9);
      lbound <- 0;
      ubound <- 1;
    } else if (distribution == "exponential_mix") {
      rate1 = as.numeric(varattrs["lambda"]);
      rate2 = as.numeric(varattrs["lambda2"]);
      mix = as.numeric(varattrs["mix"]);
      minx <- 0
      maxx <- max(qexp(0.99, rate1), qexp(0.99, rate2))
      priorx <- seq(minx, maxx, length.out=200);
      priory <- mix * dexp(priorx, rate = rate1) + (1-mix) * dexp(priorx, rate = rate2)
      lbound <- 0;
    }
  }
  
  if (plot) {
    plot(priorx, priory, col=.prior_color, main=name, xlab=xlab, ylab="Probability density", type="l", xlim=c(minx, maxx), ylim=ylim, lwd=2);
  } else {
    result <- list()
    result$xlim <- c(minx, maxx)
    result$ylim <- ylim
    result$priorx <- priorx
    result$priory <- priory
    return(result)
  }
}

variable_statistic <- function(samples, xlab="", ylim=NULL, statistic, ...)
{
  if (statistic == "mean") {
    return(mean(samples))
  }
  if (statistic == "median") {
    return(median(samples))
  }
  if (statistic == "sd") {
    return(sd(samples))
  }
  if (statistic == "quantile") {
    arguments <- list(...)
    return(quantile(samples, probs=arguments$q))
  }
  if (statistic == "autocorrelation") {
    lag <- list(...)$lag
    ac <- acf(samples, plot=F, lag.max=lag)$acf[lag+1]
    return(ac)
  }
  if (statistic == "decorr_lag") {
    ac <- acf(samples, plot=F, lag.max=length(samples)/2)
    threshold <- 2.0 / sqrt(length(samples))
    sign <- ac$acf < threshold
    return(match(T, sign))
  }
  if (statistic == "ess") {
    ac <- acf(samples, plot=F, lag.max=length(samples)/2)
    first_neg <- Position(function(x) x < 0, ac$acf[,1,1])
    if (first_neg > 2) {
      ess <- length(samples) / (1 + 2*sum(ac$acf[2:(first_neg-1),1,1]))
    } else {
      ess <- length(samples)
    }
    return(ess)
  }
}

png_tile <- function(filename, width, height, nplots)
{
  png(filename, width=width, height=height)
  nx <- ceiling(sqrt(nplots))
  ny <- ceiling(nplots / nx)
  par(mfrow=c(ny, nx))
}

pdf_tile <- function(filename, width, height, nplots)
{
  pdf(filename, width=width, height=height, useDingbats = F)
  nx <- ceiling(sqrt(nplots))
  ny <- ceiling(nplots / nx)
  par(mfrow=c(ny, nx))
}

.write_variables_to_xml <- function(model, xml_root, var_ix) {
  for (i in var_ix) {
    var <- XML::xmlNode("variable")
    logspace <- model$prior$variable_attrs[[i]]["logspace"]
    if (is.na(logspace)) {
      XML::xmlAttrs(var) <- c(name = model$variables[i])
    } else {
      XML::xmlAttrs(var) <- c(name = model$variables[i], logspace = paste(logspace))
    }
    xml_root$children[[length(xml_root$children) + 1]] <- var
  }
  return(xml_root)
}

.write_bounds_to_xml <- function(model, xml_root, var_ix) {
  # BCM is not a CRAN package, so manually check whether the BCM R function have been loaded.
  stopifnot(exists('prior_bounds'))
  for (i in var_ix) {
    var <- XML::xmlNode("variable_bound")
    bounds <- prior_bounds(model, i, c(0,1))
    XML::xmlAttrs(var) <- c(name = model$variables[i], lower = bounds[1], upper = bounds[2])
    xml_root$children[[length(xml_root$children) + 1]] <- var
  }
  return(xml_root)
}

.write_transformations_to_xml <- function(model, xml_root, var_ix) {
  for (i in var_ix) {
    var <- XML::xmlNode("variable_transformation")
    
    distribution <- model$prior$variable_attrs[[i]]["distribution"]
    if (distribution == "normal") {
      transformation <- "none"
    } else if (distribution == "uniform") {
      a = as.numeric(model$prior$variable_attrs[[i]]["lower"]);
      b = as.numeric(model$prior$variable_attrs[[i]]["upper"]);
      if (a == 0 && b == 1) {
        transformation <- "logit"
      } else {
        transformation <- "logit_scale"
      }
    } else if (distribution == "beta") {
      transformation <- "logit"
    } else if (distribution == "exponential") {
      transformation <- "log"
    } else if (distribution == "gamma") {
      transformation <- "log"
    }
    
    XML::xmlAttrs(var) <- c(name = paste(model$variables[i]),
                            transform = transformation)
    
    if (distribution == "uniform" && transformation == "logit_scale") {
      a = as.numeric(model$prior$variable_attrs[[i]]["lower"]);
      b = as.numeric(model$prior$variable_attrs[[i]]["upper"]);
      XML::xmlAttrs(var) <- c(a = a, b = b)
    }
    
    xml_root$children[[length(xml_root$children) + 1]] <- var
  }
  return(xml_root)
}

#' Export a mvdens density object as xml file that can be used as prior in the BCM software.
#' 
#' See the mvdens R package.
#' @param bcm.model A BCM model results object obtained from one of BCM's load.model functions
#' @param fit The mvdens fit object to export
#' @param outfn Output filename
#' @export
#' @examples
#' # # This assumes that BCM is installed, an environment variable BCM_ROOT is specified, and there
#' # # is a model named "model" with output directory "output_dir" in the current working folder.
#' # source(paste(Sys.getenv("BCM_ROOT"), "/scripts/plots_functions.r", sep = ""))
#' # model <- load_sbmlpd_model("model", "output_dir")
#' # x <- model$posterior$samples[sample_ix,]
#' # p <- model$posterior$lposterior[sample_ix]
#' # bounds <- prior_bounds_all(model, q = c(0, 1))
#' # gmm <- gmm.BIC(x, K = 1:6, optimal.only = T)
#' # mvd.export.bcm(model, gmm, "posterior_gmm.xml")
mvd.export.bcm <- function(bcm.model, fit, outfn, var_ix = 1:bcm.model$nvar)
{
  # BCM is not a CRAN package, so manually check whether the BCM R function have been loaded.
  stopifnot(exists('prior_bounds_all'))
  
  xml_root <- XML::xmlNode("prior")
  xml_root <- .write_variables_to_xml(bcm.model, xml_root, var_ix)
  xml_root <- .write_bounds_to_xml(bcm.model, xml_root, var_ix)
  
  if (fit$type == "kde") {
    XML::xmlAttrs(xml_root) <- c(type = "KDE")
    comp <- XML::xmlNode("kde")
    XML::xmlAttrs(comp) <- c(samples = paste(apply(fit$x, 2, paste, collapse = ","), collapse = ";"),
                             H = paste(apply(fit$H, 2, paste, collapse = ","), collapse = ";"))
    xml_root$children[[length(xml_root$children) + 1]] <- comp
  } else if (fit$type == "kde.transformed") {
    XML::xmlAttrs(xml_root) <- c(type = "KDE")
    xml_root <- .write_transformations_to_xml(bcm.model, xml_root, var_ix)
    comp <- XML::xmlNode("kde")
    XML::xmlAttrs(comp) <- c(samples = paste(apply(fit$kde$x, 2, paste, collapse = ","), collapse = ";"),
                             H = paste(apply(fit$kde$H, 2, paste, collapse = ","), collapse = ";"))
    xml_root$children[[length(xml_root$children) + 1]] <- comp
  } else if (fit$type == "gmm") {
    XML::xmlAttrs(xml_root) <- c(type = "GMM")
    for (i in 1:fit$K) {
      comp <- XML::xmlNode("component")
      XML::xmlAttrs(comp) <- c(weight = fit$proportions[i],
                               mean = paste(fit$centers[i,], collapse = ";"),
                               covariance = paste(apply(fit$covariances[[i]], 2, paste, collapse = ","), collapse = ";"))
      xml_root$children[[length(xml_root$children) + 1]] <- comp
    }
  } else if (fit$type == "gmm.transformed") {
    XML::xmlAttrs(xml_root) <- c(type = "GMM")
    xml_root <- .write_transformations_to_xml(bcm.model, xml_root, var_ix)
    for (i in 1:fit$gmm$K) {
      comp <- XML::xmlNode("component")
      XML::xmlAttrs(comp) <- c(weight = fit$gmm$proportions[i],
                               mean = paste(fit$gmm$centers[i,], collapse = ";"),
                               covariance = paste(apply(fit$gmm$covariances[[i]], 2, paste, collapse = ","), collapse = ";"))
      xml_root$children[[length(xml_root$children) + 1]] <- comp
    }
  } else if (fit$type == "gmm.truncated") {
    XML::xmlAttrs(xml_root) <- c(type = "GMMtruncated")
    bounds <- prior_bounds_all(bcm.model, c(0, 1))
    for (i in 1:fit$K) {
      log_nc <- log(pmvnorm(lower=bounds[,1], upper=bounds[,2], mean=fit$centers[i,], sigma=fit$covariances[[i]]))
      comp <- XML::xmlNode("component")
      XML::xmlAttrs(comp) <- c(weight = fit$proportions[i],
                               mean = paste(fit$centers[i,], collapse = ";"),
                               covariance = paste(apply(fit$covariances[[i]], 2, paste, collapse = ","), collapse = ";"),
                               log_truncation_correction = log_nc)
      xml_root$children[[length(xml_root$children) + 1]] <- comp
    }
  } else if (fit$type == "vine.copula") {
    XML::xmlAttrs(xml_root) <- c(type = "VineCopula")
    for (j in 1:length(var_ix)) {
      i <- var_ix[j]
      margin_node <- XML::xmlNode("marginal")
      if (fit$marginal$type == "ecdf") {
        XML::xmlAttrs(margin_node) <- c(name = bcm.model$variables[i],
                                        type = "ecdf",
                                        bw = fit$marginal$bw[i],
                                        x = paste(sort(fit$marginal$x[,j]), collapse = ";"))
      } else if (fit$marginal$type == "ecdf.pareto") {
        XML::xmlAttrs(margin_node) <- c(name = bcm.model$variables[i],
                                        type = "ecdf-pareto",
                                        bw = fit$marginal$ecdf$bw[j],
                                        x = paste(sort(fit$marginal$ecdf$x[,j]), collapse = ";"))
        if (length(fit$marginal$lower.tails[[j]]) > 0) {
          XML::xmlAttrs(margin_node) <- c(XML::xmlAttrs(margin_node),
                                          lxi = fit$marginal$lower.tails[[j]]$xi,
                                          lbeta = fit$marginal$lower.tails[[j]]$beta,
                                          lu = fit$marginal$lower.tails[[j]]$u,
                                          lq = fit$marginal$lower.tails[[j]]$q,
                                          ld = fit$marginal$lower.tails[[j]]$d)
        }
        if (length(fit$marginal$upper.tails[[j]]) > 0) {
          XML::xmlAttrs(margin_node) <- c(XML::xmlAttrs(margin_node),
                                          uxi = fit$marginal$upper.tails[[j]]$xi,
                                          ubeta = fit$marginal$upper.tails[[j]]$beta,
                                          uu = fit$marginal$upper.tails[[j]]$u,
                                          uq = fit$marginal$upper.tails[[j]]$q,
                                          ud = fit$marginal$upper.tails[[j]]$d)
        }
      } else if (fit$marginal$type == "parametric") {
        dist <- fit$marginal$dists[[j]]
        XML::xmlAttrs(margin_node) <- c(name = bcm.model$variables[i], type = dist$type)
        if (dist$type == "normal") {
          XML::xmlAttrs(margin_node) <- c(mu = dist$mean, sigma = dist$sd)
        } else if (dist$type == "beta") {
          XML::xmlAttrs(margin_node) <- c(a = dist$shape1, b = dist$shape2, min = dist$min, max = dist$max)
        } else if (dist$type == "gamma") {
          XML::xmlAttrs(margin_node) <- c(shape = dist$shape, scale = dist$scale)
        }
      } else if (fit$marginal$type == "mixture") {
        dist <- fit$marginal$dists[[j]]
        XML::xmlAttrs(margin_node) <- c(name = bcm.model$variables[i], type = dist$type, p = paste(dist$p, collapse = ";"))
        if (dist$type == "normal") {
          XML::xmlAttrs(margin_node) <- c(mu = paste(dist$mu, collapse = ";"),
                                          sigma = paste(dist$sigma, collapse = ";"))
        } else if (dist$type == "beta") {
          XML::xmlAttrs(margin_node) <- c(a = paste(dist$a, collapse = ";"),
                                          b = paste(dist$b, collapse = ";"),
                                          min = dist$min,
                                          max = dist$max)
        } else if (dist$type == "gamma") {
          XML::xmlAttrs(margin_node) <- c(shape = paste(dist$shape, collapse = ";"),
                                          scale = paste(dist$scale, collapse = ";"))
        }
      }
      xml_root$children[[length(xml_root$children) + 1]] <- margin_node
    }
    xmlvar_rvm <- XML::xmlNode("RVineCopula")
    normalized_RVM <- VineCopula::RVineMatrixNormalize(fit$RVM)
    XML::xmlAttrs(xmlvar_rvm) <- c(names = paste(normalized_RVM$names, collapse = ";"),
                                   matrix = paste(normalized_RVM$Matrix, collapse = ";"),
                                   family = paste(normalized_RVM$family, collapse = ";"),
                                   par = paste(normalized_RVM$par, collapse = ";"),
                                   par2 = paste(normalized_RVM$par2, collapse = ";"),
                                   maxmat = paste(normalized_RVM$MaxMat, collapse = ";"),
                                   codirect = paste(as.numeric(normalized_RVM$CondDistr$direct), collapse = ";"),
                                   coindirect = paste(as.numeric(normalized_RVM$CondDistr$indirect), collapse = ";"))
    xml_root$children[[length(xml_root$children) + 1]] <- xmlvar_rvm
  } else if (fit$type == "gp") {
    XML::xmlAttrs(xml_root) <- c(type = "GaussianProcess")
    
    gp_params <- XML::xmlNode("GaussianProcess")
    XML::xmlAttrs(gp_params) <- c(kernel = fit$kernel.name,
                                  l = fit$l,
                                  s = fit$s,
                                  x = paste(apply(fit$x, 2, paste, collapse = ","), collapse = ";"),
                                  p = paste(fit$p, collapse=";"))
    xml_root$children[[length(xml_root$children) + 1]] <- gp_params
  } else if (fit$type == "resample") {
    XML::xmlAttrs(xml_root) <- c(type = "resample")
    comp <- XML::xmlNode("samples")
    XML::xmlAttrs(comp) <- c(samples = paste(apply(fit$x, 2, paste, collapse = ","), collapse = ";"),
                             p = paste(fit$p, collapse = ";"))
    xml_root$children[[length(xml_root$children) + 1]] <- comp
  } else if (fit$type == "mfa") {
    # Make a GMM out of the MFA
    XML::xmlAttrs(xml_root) <- c(type = "GMM")
    for (i in 1:fit$num_components) {
      comp <- XML::xmlNode("component")
      XML::xmlAttrs(comp) <- c(weight = fit$weights[i],
                               mean = paste(fit$factor_means[,i], collapse = ";"),
                               covariance = paste(apply(fit$BtBpD[[i]], 2, paste, collapse = ","), collapse = ";"))
      xml_root$children[[length(xml_root$children) + 1]] <- comp
    }
  } else if (fit$type == "mfa.transformed") {
    XML::xmlAttrs(xml_root) <- c(type = "GMM")
    xml_root <- .write_transformations_to_xml(bcm.model, xml_root, var_ix)
    for (i in 1:fit$mfa$num_components) {
      comp <- XML::xmlNode("component")
      XML::xmlAttrs(comp) <- c(weight = fit$mfa$weights[i],
                               mean = paste(fit$mfa$factor_means[,i], collapse = ";"),
                               covariance = paste(apply(fit$mfa$BtBpD[[i]], 2, paste, collapse = ","), collapse = ";"))
      xml_root$children[[length(xml_root$children) + 1]] <- comp
    }
  }
  
  if (is.null(outfn)) {
    return(xml_root)
  } else {
    XML::saveXML(xml_root, outfn)
  }
}

hill <- function(x, k, n) { return((x^n)/(k^n+x^n))}

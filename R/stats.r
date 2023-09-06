library(pracma)

# Retrieve the index of a given variable
model_get_var_ix <- function(model, var_name) {
  match(var_name, model$variables)
}

variable_mean <- function(model, var_ix=NULL, var_name=NULL, temperature_ix = dim(model$posterior$samples)[2], sample_ix = (dim(model$posterior$samples)[3]/2+1):dim(model$posterior$samples)[3])
{
  if (!is.null(var_name)) {
    var_ix <- model_get_var_ix(model, var_name)
  }
  if (!is.null(var_ix)) {
    varattrs <- model$prior$variable_attrs[[var_ix]]
    variable_statistic(model$posterior$samples[var_ix, temperature_ix, sample_ix], varattrs, statistic="mean")
  } else {
    stop("Either a variable index or a variable name has to be specified")
  }
}

variable_sd <- function(model, var_ix=NULL, var_name=NULL, temperature_ix = dim(model$posterior$samples)[2], sample_ix = (dim(model$posterior$samples)[3]/2+1):dim(model$posterior$samples)[3])
{
  if (!is.null(var_name)) {
    var_ix <- model_get_var_ix(model, var_name)
  }
  if (!is.null(var_ix)) {
    varattrs <- model$prior$variable_attrs[[var_ix]]
    variable_statistic(model$posterior$samples[var_ix, temperature_ix, sample_ix], varattrs, statistic="sd")
  } else {
    stop("Either a variable index or a variable name has to be specified")
  }
}

variable_median <- function(model, var_ix=NULL, var_name=NULL, temperature_ix = dim(model$posterior$samples)[2], sample_ix = (dim(model$posterior$samples)[3]/2+1):dim(model$posterior$samples)[3])
{
  if (!is.null(var_name)) {
    var_ix <- model_get_var_ix(model, var_name)
  }
  if (!is.null(var_ix)) {
    varattrs <- model$prior$variable_attrs[[var_ix]]
    variable_statistic(model$posterior$samples[var_ix, temperature_ix, sample_ix], varattrs, statistic="median")
  } else {
    stop("Either a variable index or a variable name has to be specified")
  }
}

variable_quantile <- function(model, var_ix=NULL, var_name=NULL, temperature_ix = dim(model$posterior$samples)[2], sample_ix = (dim(model$posterior$samples)[3]/2+1):dim(model$posterior$samples)[3], q)
{
  if (!is.null(var_name)) {
    var_ix <- model_get_var_ix(model, var_name)
  }
  if (!is.null(var_ix)) {
    varattrs <- model$prior$variable_attrs[[var_ix]]
    variable_statistic(model$posterior$samples[var_ix, temperature_ix, sample_ix], varattrs, statistic="quantile", q=q)
  } else {
    stop("Either a variable index or a variable name has to be specified")
  }
}

variable_autocorrelation <- function(model, var_ix=NULL, var_name=NULL, temperature_ix = dim(model$posterior$samples)[2], sample_ix = (dim(model$posterior$samples)[3]/2+1):dim(model$posterior$samples)[3], lag)
{
  if (!is.null(var_name)) {
    var_ix <- model_get_var_ix(model, var_name)
  }
  if (!is.null(var_ix)) {
    varattrs <- model$prior$variable_attrs[[var_ix]]
    variable_statistic(model$posterior$samples[var_ix, temperature_ix, sample_ix], varattrs, statistic="autocorrelation", lag=lag)
  } else {
    stop("Either a variable index or a variable name has to be specified")
  }
}

variable_decorrelation_lag <- function(model, var_ix=NULL, var_name=NULL, temperature_ix = dim(model$posterior$samples)[2], sample_ix = (dim(model$posterior$samples)[3]/2+1):dim(model$posterior$samples)[3])
{
  if (!is.null(var_name)) {
    var_ix <- model_get_var_ix(model, var_name)
  }
  if (!is.null(var_ix)) {
    varattrs <- model$prior$variable_attrs[[var_ix]]
    variable_statistic(model$posterior$samples[var_ix, temperature_ix, sample_ix], varattrs, statistic="decorr_lag")
  } else {
    stop("Either a variable index or a variable name has to be specified")
  }
}

variable_effective_sample_size <- function(model, var_ix=NULL, var_name=NULL, temperature_ix = dim(model$posterior$samples)[2], sample_ix = (dim(model$posterior$samples)[3]/2+1):dim(model$posterior$samples)[3])
{
  if (!is.null(var_name)) {
    var_ix <- model_get_var_ix(model, var_name)
  }
  if (!is.null(var_ix)) {
    varattrs <- model$prior$variable_attrs[[var_ix]]
    variable_statistic(model$posterior$samples[var_ix, temperature_ix, sample_ix], varattrs, statistic="ess")
  } else {
    stop("Either a variable index or a variable name has to be specified")
  }
}

variable_summary <- function(model, temperature_ix = dim(model$posterior$samples)[2], sample_ix = (dim(model$posterior$samples)[3]/2+1):dim(model$posterior$samples)[3])
{
  df <- data.frame(row.names=model$variables)
  for (j in 1:model$nvar) {
    df$mean[j] <- variable_mean(model, j, temperature_ix=temperature_ix, sample_ix=sample_ix)
    df$sd[j] <- variable_sd(model, j, temperature_ix=temperature_ix, sample_ix=sample_ix)
    df$median[j] <- variable_median(model, j, temperature_ix=temperature_ix, sample_ix=sample_ix)
    df$q025[j] <- variable_quantile(model, j, temperature_ix=temperature_ix, sample_ix=sample_ix, q=0.025)
    df$q975[j] <- variable_quantile(model, j, temperature_ix=temperature_ix, sample_ix=sample_ix, q=0.975)
    df$autocorrelation_lag1[j] <- variable_autocorrelation(model, j, temperature_ix=temperature_ix, sample_ix=sample_ix, lag=1)
    df$decorrelation_lag[j] <- variable_decorrelation_lag(model, j, temperature_ix=temperature_ix, sample_ix=sample_ix)
    df$ess[j] <- variable_effective_sample_size(model, j, temperature_ix=temperature_ix, sample_ix=sample_ix)
  }
  return(df)
}

prior_bounds <- function(model, var_ix, q=c(0.01,0.99))
{
  p <- rep(NA, length(q))
  varattrs <- model$prior$variable_attrs[[var_ix]]
  priortype <- varattrs["type"]
  if (is.na(priortype) || priortype == "regular") {
    distribution <- varattrs["distribution"]
    if (distribution == "normal") {
      mu = as.numeric(varattrs["mu"]);
      sigma = as.numeric(varattrs["sigma"]);
      p <- qnorm(q, mean=mu, sd=sigma);
    } else if (distribution == "gamma") {
      k = as.numeric(varattrs["k"]);
      theta = as.numeric(varattrs["theta"]);
      p <- qgamma(q, shape=k, scale=theta);
    } else if (distribution == "uniform") {
      a = as.numeric(varattrs["lower"]);
      b = as.numeric(varattrs["upper"]);
      p <- c(a,b)
      p <- qunif(q, min=a, max=b);
    } else if (distribution == "half_cauchy") {
      scale = as.numeric(varattrs["scale"]);
      p <- qcauchy(q, scale=scale);
    } else if (distribution == "beta") {
      a = as.numeric(varattrs["a"]);
      b = as.numeric(varattrs["b"]);
      p <- qbeta(q, a, b)
    } else if (distribution == "exponential") {
      rate = as.numeric(varattrs["lambda"]);
      p <- qexp(q, rate=rate)
    } else if (distribution == "beta_prime") {
      a = as.numeric(varattrs["a"]);
      b = as.numeric(varattrs["b"]);
      scale = as.numeric(varattrs["scale"]);
      p <- qbetapr(q, shape1=a, shape2=b, scale=scale)
    } else if (distribution == "dirichlet") {
      # Todo
      p <- qbeta(q, shape1=1, shape2=9)
    } else if (distribution == "exponential_mix") {
    }
  }
  return(p)
}

prior_probability_var <- function(model, var_ix, x, log=F)
{
  p <- rep(NA, length(x))
  varattrs <- model$prior$variable_attrs[[var_ix]]
  priortype <- varattrs["type"]
  if (is.na(priortype) || priortype == "regular") {
    distribution <- varattrs["distribution"]
    if (distribution == "normal") {
      mu = as.numeric(varattrs["mu"]);
      sigma = as.numeric(varattrs["sigma"]);
      p <- dnorm(x, mean=mu, sd=sigma, log=log);
    } else if (distribution == "gamma") {
      k = as.numeric(varattrs["k"]);
      theta = as.numeric(varattrs["theta"]);
      p <- dgamma(x, shape=k, scale=theta, log=log);
    } else if (distribution == "uniform") {
      a = as.numeric(varattrs["lower"]);
      b = as.numeric(varattrs["upper"]);
      p <- c(a,b)
      p <- dunif(x, min=a, max=b, log=log);
    } else if (distribution == "half_cauchy") {
      scale = as.numeric(varattrs["scale"]);
      p <- dcauchy(x, scale=scale, log=log);
    } else if (distribution == "beta") {
      a = as.numeric(varattrs["a"]);
      b = as.numeric(varattrs["b"]);
      p <- dbeta(x, a, b, log=log)
    } else if (distribution == "exponential") {
      rate = as.numeric(varattrs["lambda"]);
      p <- dexp(x, rate=rate, log=log)
    } else if (distribution == "beta_prime") {
      a = as.numeric(varattrs["a"]);
      b = as.numeric(varattrs["b"]);
      scale = as.numeric(varattrs["scale"]);
      p <- dbetapr(x, shape1=a, shape2=b, scale=scale)
    } else if (distribution == "dirichlet") {
      # Todo
      p <- dbeta(x, a=1, b=9)
    } else if (distribution == "exponential_mix") {
      rate1 = as.numeric(varattrs["lambda"]);
      rate2 = as.numeric(varattrs["lambda2"]);
      mix = as.numeric(varattrs["mix"]);
      p <- mix * dexp(x, rate = rate1) + (1-mix) * dexp(x, rate = rate2)
    }
  }
  return(p)
}

prior_probability <- function(model, x, log=F)
{
  stopifnot(ncol(x) == model$nvar)
  
  p <- matrix(NA, nrow(x), ncol(x))
  for (i in 1:model$nvar) {
    p[,i] <- prior_probability_var(model, i, x[,i], log=log)
  }
  
  if (log) {
    return(apply(p, 1, sum))
  } else {
    return(apply(p, 1, prod))
  }
}

prior_bounds_all <- function(model, q=c(0.01,0.99))
{
  bounds <- matrix(NA, nrow=model$nvar, ncol=2)
  for (i in 1:model$nvar) {
    bounds[i,] <- prior_bounds(model, i, q)
  }
  return(bounds)
}

marginal_likelihood <- function(model, sample_ix=(dim(model$posterior$samples)[3]/2+1):(dim(model$posterior$samples)[3]))
{
  mean_ll <- apply(model$posterior$llikelihood[,sample_ix], 1, mean)
  if (is.infinite(mean_ll[1])) {
    return(trapz(model$posterior$temperatures[-1], mean_ll[-1]))
  } else {
    return(trapz(model$posterior$temperatures, mean_ll))
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
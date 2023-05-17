source(paste(Sys.getenv("BCM3_ROOT"), "/R/load.r", sep=""))
source(paste(Sys.getenv("BCM3_ROOT"), "/R/plots_functions.r", sep=""))

model <- bcm3.load.results(".", "output_t1_n10_e2_gmm")

sample_ix <- (dim(model$posterior$samples)[3]/2+1):(dim(model$posterior$samples)[3])
temperature_ix <- dim(model$posterior$samples)[2]

plot(t(model$posterior$samples[,temperature_ix,sample_ix]))

plot_trace(model, var_ix=1)
plot_trace(model, var_ix=2)

model_t <- bcm3.load.results(".", "output_t1_n10_e2_gmm_t")

plot(t(model_t$posterior$samples[,temperature_ix,sample_ix]))


par(mfcol=c(2,2))
plot_variable_distribution(model, 1)
plot_variable_distribution(model_t, 1)
plot_variable_distribution(model, 2)
plot_variable_distribution(model_t, 2)

source("initializeValues.R")

library(rstan)
rstan_options (auto_write = F)
options (mc.cores = parallel::detectCores ())

niter <- 2000

model <- stan_model("stanModel.stan")

fit_stan <- sampling(model,
                     data = data_fit,
                     iter = niter,
                     chains = 4, 
                     seed = 0,
                     init_r = 0.8,
                     thin = 1, 
                     control=list(adapt_delta=0.9,
                                  max_treedepth = 10))

#saveRDS(fit_stan, "data/fit_stan_Main_bin.rds")

data_fit$muA = 1.0/3.1
data_fit$muB = 1.0/3.1
data_fit$sigmaLatentA = 1.0/3.1
data_fit$sigmaLatentB = 1.0/3.1

fit_stan <- sampling(model,
                     data = data_fit,
                     iter = niter,
                     chains = 4, 
                     seed = 0,
                     init_r = 0.8,
                     thin = 1, 
                     control=list(adapt_delta=0.9,
                                  max_treedepth = 10))

#saveRDS(fit_stan, "data/fit_stan_Gen1_bin.rds")

data_fit$muA = 1.0/2.1
data_fit$muB = 1.0/2.1
data_fit$sigmaLatentA = 1.0/2.1
data_fit$sigmaLatentB = 1.0/2.1

fit_stan <- sampling(model,
                     data = data_fit,
                     iter = niter,
                     chains = 4, 
                     seed = 0,
                     init_r = 0.8,
                     thin = 1, 
                     control=list(adapt_delta=0.9,
                                  max_treedepth = 10))

#saveRDS(fit_stan, "data/fit_stan_Gen2_bin.rds")

data_fit$muA = 1.0/2.75
data_fit$muB = 1.0/2.75
data_fit$sigmaLatentA = 1.0/2.75
data_fit$sigmaLatentB = 1.0/1.85

fit_stan <- sampling(model,
                     data = data_fit,
                     iter = niter,
                     chains = 4, 
                     seed = 0,
                     init_r = 0.5,
                     thin = 1,
                     control=list(adapt_delta=0.9,
                                  max_treedepth = 10))

#saveRDS(fit_stan, "data/fit_stan_Gen3_bin.rds")

data_fit$sIHRDelta = 2.83
fit_stan <- sampling(model,
                     data = data_fit,
                     iter = niter,
                     chains = 4, 
                     seed = 0,
                     init_r = 0.8,
                     thin = 1,
                     control=list(adapt_delta=0.9,
                                  max_treedepth = 10))

#saveRDS(fit_stan, "data/fit_stan_sev_bin.rds")

data_fit$f_init = c(0.115, 0.115, 0.183, 0.183, 0.178, 0.178, 0.19, 0.19)
fit_stan <- sampling(model,
                     data = data_fit,
                     iter = 2000,
                     chains = 4, 
                     seed = 0,
                     init_r = 1.0,
                     thin = 1,
                     control=list(adapt_delta=0.9,
                                  max_treedepth = 10))
#saveRDS(fit_stan, "data/fit_stan_imm_bin.rds")

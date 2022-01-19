library(data.table)
data.dt <- fread("data/sequencing_new.csv")

data.dt$DATA <- as.Date(data.dt$DATA, format = "%d/%m/%Y")
data.dt <- data.dt[data.dt$DATA > as.Date("2021-03-01"),]

data.dt <- data.dt[data.dt$VARIANT %in% c("Delta (B.1.617.2 i llinatges AY)", "Alfa (B.1.1.7)"),]

data.dt[data.dt$VARIANT == "Alfa (B.1.1.7)",]$VARIANT <- "Alpha"
data.dt[data.dt$VARIANT ==  "Delta (B.1.617.2 i llinatges AY)",]$VARIANT <- "Delta"

data.dt <- data.dt[data.dt$PROVA == "seq\xfcenciaci\xf3 completa",]
data.dt$PROVA <- NULL

lDates <- length(unique(data.dt$DATA))
inter.dt <- data.table(
  DATA = unique(data.dt$DATA),
  VARIANT = rep("Delta", lDates),
  RECOMPTE = rep(0, lDates)
)


data.dt <- rbind(data.dt,inter.dt)

inter.dt <- data.table(
  DATA = unique(data.dt$DATA),
  VARIANT = rep("Alpha", lDates),
  RECOMPTE = rep(0, lDates)
)

data.dt <- rbind(data.dt,inter.dt)

data.dt <- data.dt |>
  dplyr::group_by(VARIANT, DATA) |>
  dplyr::summarise(RECOMPTE = sum(RECOMPTE)) |>
  tidyr::pivot_wider(names_from = VARIANT, values_from = RECOMPTE) |>
  dplyr::mutate(total = Alpha + Delta) |>
  dplyr::mutate(fracDelta = Delta / total) |>
  dplyr::mutate(n_delta = Delta) |>
  dplyr::filter(DATA >= as.Date("2021-05-01")) |>
  dplyr::filter(DATA <= as.Date("2021-10-15"))
write.table(data.dt, "data/sequenceFraction.csv")

#------------------ get_uncertainty ----------------------
library(rstan)
rstan_options (auto_write = F)
options (mc.cores = parallel::detectCores ())
stanBinomModel <- stan_model("binomial_test.stan")
n_seq <- data.dt$total
delta_seq <- data.dt$n_delta
N_l <- length(n_seq)
data_bin <- list(N = N_l, n_seq = n_seq, delta_seq = delta_seq)
fit_binom <- sampling(stanBinomModel,
                      data = data_bin,
                      iter = 2000,
                      chains = 4)
smr_pred <- cbind(as.data.frame(summary(
  fit_binom, pars = "theta", probs = c(0.025, 0.5, 0.975))$summary))
fwrite(smr_pred, "data/conf_seq.csv")

library(outbreaks)
library(bayesplot)
library(ggplot2)
library(cmdstanr)
library(posterior)

bayesplot_theme_set(theme_bw(base_size = 16))
options(mc.cores = parallel::detectCores())
cmdstanr::check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
Ndata <- 763
data(influenza_england_1978_school, package = "outbreaks")
sol <- influenza_england_1978_school
sol$time <- as.numeric(sol$date - min(sol$date)) + 2
sol$I <- sol$in_bed
noisy_I <- sol$I / Ndata
iniTime <- 0
epi.data <- list(
  n_obs = length(noisy_I),
  t0 = iniTime,
  ts = sol$time,
  y = noisy_I,
  mu_beta = 0,
  sigma_beta = 1,
  mu_gamma = 0,
  sigma_gamma = 1,
  as = 99,
  bs = 1
)

file <- file.path("/home/marcio/Projects/R0_uncertainty/code", "stan", "ssli.stan")
SIR_code <- cmdstan_model(file)

pars.interest <- c("beta", "gamma", "S0", "R0", "sigma")
yrep_vars <- function(n_obs) paste0("y_rep[", seq_len(n_obs), "]")
run_sir <- function(SIR_code, epi.data, seed = 123, chains = 4) {
  mcmc_fit <- SIR_code$sample(
    data = epi.data,
    seed = seed,
    chains = chains,
    parallel_chains = min(chains, parallel::detectCores()),
    refresh = 0
  )
  list(mcmc = mcmc_fit)
}

time_one_run <- function(SIR_code, epi.data, seed, chains = 4) {
  t0 <- Sys.time()
  fit <- run_sir(SIR_code, epi.data, seed = seed, chains = chains)
  t1 <- Sys.time()
  
  elapsed_sec <- as.numeric(difftime(t1, t0, units = "secs"))
  list(elapsed_sec = elapsed_sec)
}
s_values <- c(1, 10, 100)
n_runs <- 100
chains <- 4
set.seed(1)
timings_all <- list()
for (s_fixed in s_values) {
  epi.data$sigma_beta  <- s_fixed
  epi.data$sigma_gamma <- s_fixed
  seeds <- sample.int(1e9, n_runs)
  timings <- numeric(n_runs)
  
  for (i in seq_len(n_runs)) {
    res <- time_one_run(SIR_code, epi.data, seed = seeds[i], chains = chains)
    timings[i] <- res$elapsed_sec
  }
  
  timings_all[[as.character(s_fixed)]] <- data.frame(
    run = seq_len(n_runs),
    s = factor(s_fixed, levels = s_values),   # important for legend order
    seed = seeds,
    elapsed_sec = timings
  )
}

timings_df <- do.call(rbind, timings_all)
print(aggregate(elapsed_sec ~ s, data = timings_df, summary))
timing_hist <- ggplot(timings_df, aes(x = elapsed_sec, fill = s, color = s)) +
  geom_histogram(
    bins = 12,
    position = "identity",
    alpha = 0.35
  ) +
  scale_x_continuous("Elapsed time per run (seconds)", expand = c(0, 0)) +
  scale_y_continuous("Count", expand = c(0, 0)) +
  ggtitle(sprintf("Runtime distribution: %d runs per s (chains=%d)", n_runs, chains))

timing_hist
ggsave(
  filename = "timing_hist.jpg",
  plot = timing_hist,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)

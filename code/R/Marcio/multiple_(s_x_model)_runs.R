library(outbreaks)
library(bayesplot)
library(ggplot2)
library(cmdstanr)
library(posterior)

bayesplot_theme_set(theme_bw(base_size = 16))
options(mc.cores = parallel::detectCores())
cmdstanr::check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)

# -----------------------------
# Data loading and prep
# -----------------------------
Ndata <- 763
data(influenza_england_1978_school, package = "outbreaks")

sol <- influenza_england_1978_school
sol$time <- as.numeric(sol$date - min(sol$date)) + 2
sol$I <- sol$in_bed
noisy_I <- sol$I / Ndata
iniTime <- 0
iniI <- 1 / Ndata
epi.data_base <- list(
  n_obs = length(noisy_I),
  t0 = iniTime,
  ts = sol$time,
  y_init = iniI,
  y = noisy_I,
  as = 99,
  bs = 1
)

gamma_from_mean_sd <- function(m, sd) {
  shape <- (m / sd)^2
  scale <- (sd^2) / m
  list(shape = shape, scale = scale)
}

lognormal_from_mean_sd <- function(m, sd) {
  sigma <- sqrt(log(1 + (sd^2 / m^2)))
  mu    <- log(m) - 0.5 * sigma^2
  list(mu = mu, sigma = sigma)
}

halfnormal_scale_from_mean <- function(m) {
  # mean = tau * sqrt(2/pi)  => tau = mean * sqrt(pi/2)
  tau <- m * sqrt(pi / 2)
  list(scale = tau)
}

stan_dir <- "/home/marcio/Projects/R0_uncertainty/code/stan/Marcio"
img_dir <- "/home/marcio/Projects/R0_uncertainty/figures/Marcio"
model_files <- list(
  gamma = file.path(stan_dir, "SIR_gamma.stan"),
  halfnormal = file.path(stan_dir, "SIR_halfnormal.stan"),
  lognormal = file.path(stan_dir, "SIR_lognormal.stan")
)

models <- lapply(model_files, cmdstan_model)
run_sir <- function(SIR_code, epi.data, seed = 123, chains = 4) {
  mcmc_fit <- SIR_code$sample(
    data = epi.data,
    seed = seed,
    chains = chains,
    parallel_chains = min(chains, parallel::detectCores()),
    refresh = 0
  )
  mcmc_fit
}

time_one_run <- function(SIR_code, epi.data, seed, chains = 4) {
  t0 <- Sys.time()
  invisible(run_sir(SIR_code, epi.data, seed = seed, chains = chains))
  t1 <- Sys.time()
  as.numeric(difftime(t1, t0, units = "secs"))
}

s_values <- c(1, 10, 100)
n_runs <- 100
chains <- 4
set.seed(1)

for (model_name in names(models)) {
  SIR_code <- models[[model_name]]
  timings_all <- list()
  for (s_fixed in s_values) {
    target_mean <- 1e-6
    target_sd <- s_fixed
    epi.data <- epi.data_base
    if (model_name == "gamma") {
      # beta ~ gamma(kb, thetab) with mean=kb*thetab and sd=sqrt(kb)*thetab
      p_beta <- gamma_from_mean_sd(target_mean, target_sd)
      p_gamma <- gamma_from_mean_sd(target_mean, target_sd)
      epi.data$kb <- p_beta$shape
      epi.data$thetab <- p_beta$scale
      epi.data$kg <- p_gamma$shape
      epi.data$thetag <- p_gamma$scale
      
    } else if (model_name == "lognormal") {
      # beta ~ lognormal(kb, thetab) where kb=mu, thetab=sigma on log scale
      p_beta <- lognormal_from_mean_sd(target_mean, target_sd)
      p_gamma <- lognormal_from_mean_sd(target_mean, target_sd)
      epi.data$kb <- p_beta$mu
      epi.data$thetab <- p_beta$sigma
      epi.data$kg <- p_gamma$mu
      epi.data$thetag <- p_gamma$sigma
      
    } else if (model_name == "halfnormal") {
      # beta ~ normal(0, thetab) T[0,]  (half-normal), match mean only
      p_beta <- halfnormal_scale_from_mean(target_mean)
      p_gamma <- halfnormal_scale_from_mean(target_mean)
      epi.data$thetab <- p_beta$scale
      epi.data$thetag <- p_gamma$scale
    }
    seeds <- sample.int(1e9, n_runs)
    timings <- numeric(n_runs)
    for (i in seq_len(n_runs)) {
      timings[i] <- time_one_run(SIR_code, epi.data, seed = seeds[i], chains = chains)
    }
    timings_all[[as.character(s_fixed)]] <- data.frame(
      run = seq_len(n_runs),
      s = factor(s_fixed, levels = s_values),
      seed = seeds,
      elapsed_sec = timings
    )
  }
  
  timings_df <- do.call(rbind, timings_all)
  print(paste("Model:", model_name))
  print(aggregate(elapsed_sec ~ s, data = timings_df, summary))

  timing_hist <- ggplot(timings_df, aes(x = elapsed_sec, fill = s, color = s)) +
    geom_histogram(
      bins = 12,
      position = "identity",
      alpha = 0.35
    ) +
    scale_x_continuous("Elapsed time per run (seconds)", expand = c(0, 0)) +
    scale_y_continuous("Count", expand = c(0, 0)) +
    ggtitle(sprintf(
      "Runtime distribution: %s.stan (%d runs per s, chains=%d)",
      model_name, n_runs, chains
    ))
  print(timing_hist)
  filename <- file.path(img_dir, sprintf("timing_hist_%s.jpg", model_name))
  ggsave(
    filename = filename,
    plot = timing_hist,
    width = 8,
    height = 6,
    units = "in",
    dpi = 300
  )
}

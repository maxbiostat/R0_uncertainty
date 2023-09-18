library(outbreaks)
library(bayesplot)
library(cmdstanr)
library(rstan)
stanfit <- function(fit) rstan::read_stan_csv(fit$output_files())

########## Data loading and prep

data(influenza_england_1978_school)
Ndata <- 763
sol <- influenza_england_1978_school
sol$time <- as.numeric(sol$date-min(sol$date)) + 2
sol$I <- sol$in_bed
forfit.sol <- sol
noisy_I <- forfit.sol$I/Ndata
iniTime <- 0
iniI <- 1/Ndata

epi.data <- list(
  n_obs = length(noisy_I),
  t0 = iniTime,
  ts = forfit.sol$time,
  y_init = iniI,
  y = noisy_I,
  mu_beta = 0,
  sigma_beta = 1,
  mu_gamma = 0,
  sigma_gamma = 1,
  as = 9, #254,
  bs = 1#350-254
)

########## Model fitting

options(mc.cores = parallel::detectCores())
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)

file <- file.path("../stan", "sir_simple_I.stan")
SIR_code <- cmdstan_model(file)

file_log <- file.path("../stan", "sir_simple_logI.stan")
SIR_log_code <- cmdstan_model(file_log)


pars.interest <- c("beta", "gamma", "S0", "R0", "sigma")

SIR.map <- SIR_code$optimize(data = epi.data, seed = 123)

fit_mcmc <- SIR_code$sample(data = epi.data,
                            seed = 666,
                            iter_warmup = 1500,
                            iter_sampling = 500,
                            chains = 100, parallel_chains = 4,
                            refresh = 0)

write.csv(fit_mcmc$time()$chains, 
          file = "../../data/saved_data/runtimes_real_SIR_bschool.csv",
          row.names = FALSE)

SIR.posterior <- fit_mcmc
SIR.posterior$cmdstan_diagnose()
SIR.posterior$cmdstan_summary()

draws.real <- lapply(pars.interest, function(par) {
  raw <- posterior::extract_variable_matrix(SIR.posterior, par)
  all.ess <- apply(raw, 2, posterior::ess_basic)
  out <- data.frame(ess = all.ess, parameter = par)
  out$chain <- paste0("chain_", 1:nrow(out))
  return(out)
})
draws.real.dt <- do.call(rbind, draws.real)

write.csv(draws.real.dt, 
          file = "../../data/saved_data/ESS_real_SIR_bschool.csv",
          row.names = FALSE)

print(stanfit(SIR.posterior), pars = pars.interest)

mcmc_pairs(SIR.posterior$draws(pars.interest))

p <- mcmc_trace(SIR.posterior$draws(pars.interest))
p + facet_text(size = 15)

simulated_trajectories <- SIR.posterior$draws('y_rep', format = "df")[, c(1:14)]
predicted <- data.frame(time = epi.data$ts,
                           lower = apply(simulated_trajectories, 
                           2, function(x) as.numeric(quantile(x, probs = .025)))[1:14],
                           post_mean = colMeans(simulated_trajectories),
                           post_median = apply(simulated_trajectories, 2, median),
                           upper = apply(simulated_trajectories, 
                           2, function(x) as.numeric(quantile(x, probs = .975)))[1:14],
                           parametrisation = "real")

### Log 'parametrisation'

SIR_log.map <- SIR_log_code$optimize(data = epi.data, seed = 123)

fit_mcmc_log <- SIR_log_code$sample(data = epi.data, seed = 123,
                                    
                            chains = 100, parallel_chains = 4,
                            refresh = 0)

write.csv(fit_mcmc_log$time()$chains, 
          file = "../../data/saved_data/runtimes_log_SIR_bschool.csv",
          row.names = FALSE)

SIR_log.posterior <- fit_mcmc_log
SIR_log.posterior$cmdstan_diagnose()
SIR_log.posterior$cmdstan_summary()

draws.log <- lapply(pars.interest, function(par) {
  raw <- posterior::extract_variable_matrix(SIR_log.posterior, par)
  all.ess <- apply(raw, 2, posterior::ess_basic)
  out <- data.frame(ess = all.ess, parameter = par)
  out$chain <- paste0("chain_", 1:nrow(out))
  return(out)
})
draws.log.dt <- do.call(rbind, draws.log)

write.csv(draws.log.dt, 
          file = "../../data/saved_data/ESS_log_SIR_bschool.csv",
          row.names = FALSE)

print(stanfit(SIR_log.posterior), pars = pars.interest)

mcmc_pairs(SIR_log.posterior$draws(pars.interest))

p <- mcmc_trace(SIR_log.posterior$draws(pars.interest))
p + facet_text(size = 15)

simulated_trajectories_log <- SIR_log.posterior$draws('y_rep',
                                                      format = "df")[, c(1:14)]

predicted_log <- data.frame(time = epi.data$ts,
                           lower = apply(simulated_trajectories_log, 
                                         2, function(x) as.numeric(quantile(x, probs = .025)))[1:14],
                           post_mean = colMeans(simulated_trajectories_log),
                           post_median = apply(simulated_trajectories_log, 2, median),
                           upper = apply(simulated_trajectories_log, 
                                         2, function(x) as.numeric(quantile(x, probs = .975)))[1:14],
                           parametrisation = "log")


#### Plotting and annotating

prediction.bands.SIR <- do.call(rbind,
                                list(predicted,
                                     predicted_log))

library(ggplot2)

predictions_SIR <- ggplot(data = prediction.bands.SIR,
                          aes(x = time, y = post_mean,
                              colour = parametrisation,
                              fill = parametrisation)) +
  geom_line() +
  geom_point(data = data.frame(time = epi.data$ts, I = epi.data$y),
             aes(x = time, y = I), inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .2) +
  scale_x_continuous("Time", expand = c(0, 0)) + 
  scale_y_continuous(expression(I(t)), expand = c(0, 0)) + 
  ggtitle("Mean prediction") +
  theme_bw(base_size = 16)

predictions_SIR

predictions_SIR_median <- ggplot(data = prediction.bands.SIR,
                          aes(x = time, y = post_median,
                              colour = parametrisation,
                              fill = parametrisation)) +
  geom_line() +
  geom_point(data = data.frame(time = epi.data$ts, I = epi.data$y),
             aes(x = time, y = I), inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .2) +
  scale_x_continuous("Time", expand = c(0, 0)) + 
  scale_y_continuous(expression(I(t)), expand = c(0, 0)) + 
  ggtitle("Median prediction") +
  theme_bw(base_size = 16)

predictions_SIR_median

R0 <- data.frame(R0 = SIR.posterior$draws('R0', format = "df"),
                 parametrisation = "real")
R0_log <- data.frame(R0 = SIR_log.posterior$draws('R0', format = "df"),
                     parametrisation = "log")

R0.posteriors <- do.call(rbind, list(R0, R0_log))
R0.posteriors <- R0.posteriors[, c(1:5)]

### Posterior quantities for R0 (mean, median, sd)
aggregate(R0.R0~parametrisation, mean, data = R0.posteriors)
aggregate(R0.R0~parametrisation, median, data = R0.posteriors)
aggregate(R0.R0~parametrisation, sd, data = R0.posteriors)

R0_posterior <- ggplot(data = R0.posteriors,
                       aes(x = R0.R0,
                           colour = parametrisation,
                           fill = parametrisation)) +
  geom_density(alpha = .4) +
  # geom_vline(xintercept = 1.5, linetype = "dotted", size = 1.01) + 
  stat_function(fun = function(x) dlnorm(x, meanlog = 0,
                                         sdlog = sqrt(epi.data$sigma_beta^2 + epi.data$sigma_gamma^2)),
                inherit.aes = FALSE, linetype = "longdash", size = 1.10) +
  stat_function(fun = function(x) dlnorm(x, meanlog = 0,
                                         sdlog = sqrt(epi.data$sigma_beta^2 + epi.data$sigma_gamma^2)),
                inherit.aes = FALSE, linetype = "twodash", size = 1.10) +
  stat_function(fun = function(x) dlnorm(x, meanlog = 0,
                                         sdlog = sqrt(epi.data$sigma_beta^2 + epi.data$sigma_gamma^2)),
                inherit.aes = FALSE, linetype = "F1", size = 1.10) +
  scale_x_continuous(expression(R[0]), expand = c(0, 0),
                     limits = c(0, 10)) + 
  scale_y_continuous("Density", expand = c(0, 0)) + 
  theme_bw(base_size = 16)

R0_posterior

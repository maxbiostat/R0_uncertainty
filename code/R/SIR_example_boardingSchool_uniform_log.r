library(outbreaks)
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
  as = 9, #254,
  bs = 1#350-254
)
plot(epi.data$ts, epi.data$y)

#### Model fitting

library(cmdstanr)
options(mc.cores = parallel::detectCores())
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
file <- file.path("../stan", "sir_simple_I(t)_uniform_log.stan")
SIR_code <- cmdstan_model(file)

SIR.map <- SIR_code$optimize(data = epi.data, seed = 123)
fit_mcmc <- SIR_code$sample(data = epi.data, seed = 123,
                            chains = 4, parallel_chains = 4,
                            max_treedepth = 13,
                            adapt_delta = 0.95)
SIR.posterior <- fit_mcmc
SIR.posterior$cmdstan_diagnose()
SIR.posterior$cmdstan_summary()

library("bayesplot")
print(SIR.posterior$draws(c("beta", "gamma", "S0", "R0", "sigma")))

mcmc_pairs(SIR.posterior$draws(c("beta", "gamma", "S0", "R0", "sigma")))

p <- mcmc_trace(SIR.posterior$draws(c("beta", "gamma", "S0", "R0", "sigma")))
p + facet_text(size = 15)
simulated_trajectories <- SIR.posterior$draws('y_rep', format = "df")[,c(1:14)]
predicted <- data.frame(time = epi.data$ts,
                           lower = apply(simulated_trajectories, 
                           2, function(x) as.numeric(quantile(x, probs = .025)))[1:14],
                           post_mean = colMeans(simulated_trajectories),
                           post_median = apply(simulated_trajectories, 2, median),
                           upper = apply(simulated_trajectories, 
                           2, function(x) as.numeric(quantile(x, probs = .975)))[1:14])

#### Plotting and annotating

prediction.bands.SIR <- predicted

library(ggplot2)

predictions_SIR <- ggplot(data = prediction.bands.SIR,
                          aes(x = time, y = post_mean)) +
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
                          aes(x = time, y = post_median)) +
  geom_line() +
  geom_point(data = data.frame(time = epi.data$ts, I = epi.data$y),
             aes(x = time, y = I), inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .2) +
  scale_x_continuous("Time", expand = c(0, 0)) + 
  scale_y_continuous(expression(I(t)), expand = c(0, 0)) + 
  ggtitle("Median prediction") +
  theme_bw(base_size = 16)

predictions_SIR_median

R0.posteriors <- SIR.posterior$draws('R0', format = "df")

hist(R0.posteriors$R0)


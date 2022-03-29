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
  mu_beta = 0,
  sigma_beta = 1,
  mu_gamma = 0,
  sigma_gamma = 1,
  as = 9, #254,
  bs = 1#350-254
)
plot(epi.data$ts, epi.data$y)

#### Model fitting

library(cmdstanr)
# rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
file <- file.path("stan", "sir_simple_I(t).stan")
SIR_code <- cmdstan_model(file)
# SIR_code <- stan_model(file = "stan/sir_simple_I(t).stan")

### Varying alpha

# SIR.map.s1 <- optimize(SIR_code, data = epi.data, hessian = TRUE, verbose = TRUE)
SIR.map.s1 <- SIR_code$optimize(data = epi.data, seed = 123)
# SIR.posterior.s1 <- sampling(SIR_code, data = epi.data, chains = 4, control = list(adapt_delta = .99))
fit_mcmc <- SIR_code$sample(data = epi.data, seed = 123, chains = 4, parallel_chains = 4)
# check_hmc_diagnostics(SIR.posterior.s1)
SIR.posterior.s1 <- fit_mcmc
SIR.posterior.s1$cmdstan_diagnose()
SIR.posterior.s1$cmdstan_summary()

library("bayesplot")
print(SIR.posterior.s1$draws(c("beta", "gamma", "S0", "R0", "sigma")))
# pairs(SIR.posterior.s1$draws(c("beta", "gamma", "S0", "R0", "sigma")))
mcmc_pairs(SIR.posterior.s1$draws(c("beta", "gamma", "S0", "R0", "sigma")))
# stan_trace(SIR.posterior.s1, pars = c("beta", "gamma", "S0", "R0", "sigma"))
p <- mcmc_trace(SIR.posterior.s1$draws(c("beta", "gamma", "S0", "R0", "sigma")))
p + facet_text(size = 15)
# simulated_trajectories.s1 <- extract(SIR.posterior.s1, 'y_rep')$y_rep
simulated_trajectories.s1 <- SIR.posterior.s1$draws('y_rep', format = "df")[,c(1:14)]
predicted.s1 <- data.frame(time = epi.data$ts,
                           lower = apply(simulated_trajectories.s1, 
                           2, function(x) as.numeric(quantile(x, probs = .025)))[1:14],
                           post_mean = colMeans(simulated_trajectories.s1),
                           upper = apply(simulated_trajectories.s1, 
                           2, function(x) as.numeric(quantile(x, probs = .975)))[1:14],
                           s = "1")

### sigma_beta = sigma_gamma = 10
epi.data$sigma_beta <- 10
epi.data$sigma_gamma <- 10
SIR.map.s10 <- SIR_code$optimize(data = epi.data, seed = 123)
fit_mcmc <- SIR_code$sample(data = epi.data, seed = 123, chains = 4, parallel_chains = 4)
SIR.posterior.s10 <- fit_mcmc
SIR.posterior.s10$cmdstan_diagnose()
SIR.posterior.s10$cmdstan_summary()
print(SIR.posterior.s10$draws(c("beta", "gamma", "S0", "R0", "sigma")))
mcmc_pairs(SIR.posterior.s10$draws(c("beta", "gamma", "S0", "R0", "sigma")))
p <- mcmc_trace(SIR.posterior.s10$draws(c("beta", "gamma", "S0", "R0", "sigma")))
p + facet_text(size = 15)
simulated_trajectories.s10 <- SIR.posterior.s10$draws('y_rep', format = "df")[,c(1:14)]
dim(simulated_trajectories.s10)
predicted.s10 <- data.frame(time = epi.data$ts,
                           lower = apply(simulated_trajectories.s10, 
                                         2, function(x) as.numeric(quantile(x, probs = .025)))[1:14],
                           post_mean = colMeans(simulated_trajectories.s10),
                           upper = apply(simulated_trajectories.s10, 
                                         2, function(x) as.numeric(quantile(x, probs = .975)))[1:14],
                           s = "10")

### sigma_beta = sigma_gamma = 100
epi.data$sigma_beta <- 100
epi.data$sigma_gamma <- 100

SIR.map.s100 <- SIR_code$optimize(data = epi.data, seed = 123)
fit_mcmc <- SIR_code$sample(data = epi.data, seed = 123, chains = 4, parallel_chains = 4)
SIR.posterior.s100 <- fit_mcmc
SIR.posterior.s100$cmdstan_diagnose()
SIR.posterior.s100$cmdstan_summary()
print(SIR.posterior.s100$draws(c("beta", "gamma", "S0", "R0", "sigma")))
mcmc_pairs(SIR.posterior.s100$draws(c("beta", "gamma", "S0", "R0", "sigma")))
p <- mcmc_trace(SIR.posterior.s100$draws(c("beta", "gamma", "S0", "R0", "sigma")))
p + facet_text(size = 15)
simulated_trajectories.s100 <- SIR.posterior.s100$draws('y_rep', format = "df")[,c(1:14)]
dim(simulated_trajectories.s100)
predicted.s100 <- data.frame(time = epi.data$ts,
                            lower = apply(simulated_trajectories.s100, 
                                          2, function(x) as.numeric(quantile(x, probs = .025)))[1:14],
                            post_mean = colMeans(simulated_trajectories.s100),
                            upper = apply(simulated_trajectories.s100, 
                                          2, function(x) as.numeric(quantile(x, probs = .975)))[1:14],
                            s = "100")

#### Plotting and annotating

prediction.bands.SIR <- do.call(rbind, list(predicted.s1, predicted.s10, predicted.s100))

library(ggplot2)

predictions_SIR <- ggplot(data = prediction.bands.SIR, aes(x = time, y = post_mean, colour = s, fill = s)) +
  geom_line() +
  geom_point(data = data.frame(time = epi.data$ts, I = epi.data$y), aes(x = time, y = I), inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .2) +
  scale_x_continuous("Time", expand = c(0, 0)) + 
  scale_y_continuous(expression(I(t)), expand = c(0, 0)) + 
  theme_bw(base_size = 16)

predictions_SIR

# R0.s1 <- data.frame(R0 = extract(SIR.posterior.s1, 'R0')$R0, s = "1")
R0.s1 <- data.frame(R0 = SIR.posterior.s1$draws('R0', format = "df"), s="1")
R0.s10 <- data.frame(R0 = SIR.posterior.s10$draws('R0', format = "df"), s="10")
R0.s100 <- data.frame(R0 = SIR.posterior.s100$draws('R0', format = "df"), s="100")

R0.posteriors <- do.call(rbind, list(R0.s1, R0.s10, R0.s100))
R0.posteriors <- R0.posteriors[,c(1:3)]

R0_posterior <- ggplot(data = subset(R0.posteriors, s != "100"), aes(x = R0, colour = s, fill = s)) +
  geom_density(alpha = .4) +
  geom_vline(xintercept = 1.5, linetype = "dotted", size = 1.01) + 
  stat_function(fun = function(x) dlnorm(x, meanlog = 0, sdlog = sqrt(epi.data$sigma_beta^2 + epi.data$sigma_gamma^2)),
                inherit.aes = FALSE, linetype = "longdash", size = 1.10) +
  stat_function(fun = function(x) dlnorm(x, meanlog = 0, sdlog = sqrt(epi.data$sigma_beta^2 + epi.data$sigma_gamma^2)),
                inherit.aes = FALSE, linetype = "twodash", size = 1.10) +
  stat_function(fun = function(x) dlnorm(x, meanlog = 0, sdlog = sqrt(epi.data$sigma_beta^2 + epi.data$sigma_gamma^2)),
                inherit.aes = FALSE, linetype = "F1", size = 1.10) +
  scale_x_continuous(expression(R[0]), expand = c(0, 0)) + 
  scale_y_continuous("Density", expand = c(0, 0)) + 
  theme_bw(base_size = 16)
R0_posterior

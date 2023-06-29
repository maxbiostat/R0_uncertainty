library(outbreaks)
library(bayesplot)
library(cmdstanr)
library(rstan)
stanfit <- function(fit) rstan::read_stan_csv(fit$output_files())
library(priorsense)

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
file <- file.path("../stan", "sir_simple_I_psa.stan")
SIR_code <- cmdstan_model(file)

pars.interest <- c("beta", "gamma", "S0", "R0", "sigma")

SIR.map.s1 <- SIR_code$optimize(data = epi.data, seed = 123)

fit_mcmc <- SIR_code$sample(data = epi.data, seed = 123,
                            chains = 4, parallel_chains = 4,
                            refresh = 0)

SIR.posterior.s1 <- fit_mcmc
SIR.posterior.s1$cmdstan_diagnose()
SIR.posterior.s1$cmdstan_summary()
fit.s1 <- stanfit(SIR.posterior.s1)


print(fit.s1, pars = pars.interest)

psa.s1 <- powerscale_sensitivity(fit.s1)
subset(psa.s1$sensitivity, variable %in% pars.interest)
pss.s1 <- powerscale_sequence(fit.s1)  


psa.ecdf <- priorsense::powerscale_plot_ecdf(pss.s1,
                     variables = pars.interest)

psa.ecdf

ggsave(
  plot = psa.ecdf,
  filename = "../../figures/PSA_simple_SIR_I.pdf",
  scale = 1,
  width = 297,
  height = 210,
  units = "mm",
  dpi = 300
)

priorsense::powerscale_plot_dens(pss.s1, variables = pars.interest)


mcmc_pairs(SIR.posterior.s1$draws(pars.interest))

p <- mcmc_trace(SIR.posterior.s1$draws(pars.interest))
p + facet_text(size = 15)
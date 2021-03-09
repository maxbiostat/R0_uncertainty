
true.gamma <- .2
true.R0 <- 1.15
true.beta <- true.gamma * true.R0
true.sigma <- .10

library(deSolve)
sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I
    dI <-  beta * S * I - gamma * I
    dR <-                 gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

true.s0 <- .92
init       <- c(S = true.s0, I = 1-true.s0, R = 0.0)
parameters <- c(beta = true.beta, gamma = true.gamma)
times      <- seq(0, 99, by = .5)
sol <- as.data.frame(ode(y = init, times = times, func = sir, parms = parameters)) 
sol <- sol[-1, ]

N <- nrow(sol)
noisy_R <- rlnorm(N, mean = log(sol$R), sd = true.sigma)
# hamper <- TRUE
# if(hamper){
#   noisy_R <- ifelse(noisy_R >= 1, 1, noisy_R)
# }

plot(noisy_R)

epi.data <- list(
  n_obs = length(noisy_R),
  ts = sol$time,
  y_init = noisy_R[1],
  y = noisy_R,
  ab = 1,
  bb = 1,
  ag = 1,
  bg = 1,
  ar = 1, 
  br = 1,
  as = 1,
  bs = 1
)
#####

library(rstan)
rstan_options(auto_write = TRUE)
SIR_code <- stan_model(file = "stan/simple_SIR.stan")

SIR.map <- optimizing(SIR_code, data = epi.data)
SIR.posterior <- sampling(SIR_code, data = epi.data, chains = 4)
print(SIR.posterior, pars = c("beta", "gamma", "s0", "R0"))
pairs(SIR.posterior, pars = c("beta", "gamma", "s0", "R0"))
stan_trace(SIR.posterior, pars = c("beta", "gamma", "s0", "R0"))
                          
simulated_trajectories <- extract(SIR.posterior, 'y_rep')$y_rep
measured_and_predicted <- data.frame(
  lower = apply(simulated_trajectories, 2, function(x) as.numeric(quantile(x, probs = .025))),
  R = noisy_R,
  time = sol$time,
  post_mean = colMeans(simulated_trajectories),
  upper = apply(simulated_trajectories, 2, function(x) as.numeric(quantile(x, probs = .975)))
)

ggplot(data = measured_and_predicted, aes(x = time, y = R)) +
  geom_point(alpha = .4) +
  geom_line(aes(x = time, y = post_mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .2) +
  scale_x_continuous("Time", expand = c(0, 0)) + 
  scale_y_continuous("R(t)", expand = c(0, 0)) + 
  theme_bw()

########Prelimniary stuff #########################################
LogMean <- function(realMean, realSD){
  ## takes REAL mean and REAL standard deviation; returns LOG mean
  realVar <- realSD^2
  mu <- log(realMean/ sqrt(1 + (realVar/realMean^2)))
  return(mu)
}
LogVar <- function(realMean, realSD){
  ## takes REAL mean and REAL standard deviation; returns LOG variance
  realVar <- realSD^2
  sigmasq <- log(1 + (realVar/realMean^2))
  return(sigmasq)
}
##
exact_I_max <- function(S0, I0, R_0){
  ans <- S0 + I0 - 1/R_0 * ( log(S0) + 1 + log(R_0))
  return(ans)
}
###################################################################
###################################################################

### Elicitation

mb <- 2 # E[beta]
vb <- 1^2 # Var(beta)

mg <- .4 # E[gamma]
vg <- .5^2 # Var(gamma)

# Gamma
alphab <- mb^2/vb
betab <- mb/vb

alphag <- mg^2/vg
betag <- mg/vg
qgamma(p = c(.025, .975), shape = alphag, rate = betag)

# Log-normal
mub <- LogMean(realMean = mb, realSD = sqrt(vb) )
sigmab <- sqrt( LogVar(realMean = mb, realSD = sqrt(vb) ))

mug <- LogMean(realMean = mg, realSD = sqrt(vg) )
sigmag <- sqrt( LogVar(realMean = mg, realSD = sqrt(vg) ))

qlnorm(p = c(.025, .975), meanlog = mub, sdlog = sigmab)
qlnorm(p = c(.025, .975), meanlog = mug, sdlog = sigmag)
### Comparing densities

minRb <- 0
maxRb <- 5

curve(dlnorm(x, meanlog = mub, sdlog = sigmab), lwd = 3,
      ylab = "Density", xlab = expression(beta), minRb, maxRb)
curve(dgamma(x, shape = alphab, rate = betab), lwd = 3, lty = 2, add = TRUE)
legend(x = "topright", legend = c("Gamma", "Log-normal"), lty = 2:1, bty = 'n')


minRg <- 0
maxRg <- 1
curve(dlnorm(x, meanlog = mug, sdlog = sigmag), lwd = 3,
      ylab = "Density", xlab = expression(gamma), minRg, maxRg)
curve(dgamma(x, shape = alphag, rate = betag), lwd = 3, lty = 2, add = TRUE)
legend(x = "topright", legend = c("Gamma", "Log-normal"), lty = 2:1, bty = 'n')


## Sampling 
M <- 1E6

beta.samples.G <- rgamma(n = M, shape = alphab, rate = betab)
gamma.samples.G <- rgamma(n = M, shape = alphag, rate = betag)

R0.samples.G <- beta.samples.G/gamma.samples.G

beta.samples.L <- rlnorm(n = M, meanlog = mub, sdlog = sigmab)
gamma.samples.L <- rlnorm(n = M, meanlog = mug, sdlog = sigmag)

R0.samples.L <- beta.samples.L/gamma.samples.L

beta.samples.H <- abs(rnorm(n = M, mean = mb, sd = sqrt(vb)))
gamma.samples.H <- abs(rnorm(n = M, mean = mg, sd = sqrt(vg)))

R0.samples.H <- beta.samples.H/gamma.samples.H

N <- 1
I0a <- 1E-7
I0b <- 1E-6
I0c <- 1E-5

Imax.samples.I0a.G <- exact_I_max(R_0 = R0.samples.G, I0 = I0a, S0 = N-I0a)
Imax.samples.I0a.L <-  exact_I_max(R_0 = R0.samples.L, I0 = I0a, S0 = N-I0a)
Imax.samples.I0a.H <-  exact_I_max(R_0 = R0.samples.H, I0 = I0a, S0 = N-I0a)

Imax.samples.I0b.G <- exact_I_max(R_0 = R0.samples.G, I0 = I0b, S0 = N-I0b)
Imax.samples.I0b.L <-  exact_I_max(R_0 = R0.samples.L, I0 = I0b, S0 = N-I0b)
Imax.samples.I0b.H <-  exact_I_max(R_0 = R0.samples.H, I0 = I0b, S0 = N-I0b)

Imax.samples.I0c.G <- exact_I_max(R_0 = R0.samples.G, I0 = I0c, S0 = N-I0c)
Imax.samples.I0c.L <-  exact_I_max(R_0 = R0.samples.L, I0 = I0c, S0 = N-I0c)
Imax.samples.I0c.H <-  exact_I_max(R_0 = R0.samples.H, I0 = I0c, S0 = N-I0c)

Imax.dt <- data.frame(
  Imax = c(Imax.samples.I0a.G, Imax.samples.I0b.G, Imax.samples.I0c.G,
           Imax.samples.I0a.L, Imax.samples.I0b.L, Imax.samples.I0c.L,
           Imax.samples.I0a.H, Imax.samples.I0b.H, Imax.samples.I0c.H),
  prior = rep(c("Gamma", "Lognormal", "Halfnormal"), each = 3*M),
  I0 = rep(c(I0a, I0b, I0c), each = M)
)
Imax.dt$I0 <- as.factor(Imax.dt$I0)

library(ggplot2)

# p0 <- ggplot(data = Imax.dt, aes(x = Imax, col = prior, fill = prior)) + 
#   geom_density(alpha = .2)+
#   scale_x_log10(expression(log(I[max])), expand = c(0, 0)) + 
#   scale_y_continuous("Density", expand = c(0, 0)) + 
#   geom_vline(xintercept = 1/2, linetype = "longdash") +
#   facet_wrap(I0~.) + 
#   theme_bw(base_size = 16)
# p0

p0 <- ggplot(data = Imax.dt,
             aes(x = Imax, col = prior, fill = prior)) + 
  geom_density(alpha = .2)+
  scale_x_continuous(expression(I[max]), expand = c(0, 0), limit = c(0, 1)) + 
  scale_y_continuous("Density", expand = c(0, 0)) + 
  geom_vline(xintercept = 1/2, linetype = "longdash") +
  facet_wrap(I0~.) + 
  theme_bw(base_size = 16)
p0

ggsave(plot = p0,
       filename = "../../figures/Imax_boarding_GammaLNHN.pdf",
       scale = 1,
       width = 297,
       height = 210,
       units = "mm",
       dpi = 300)

# p1 <- ggplot(data = Imax.dt, aes(x = Imax, col = I0, fill = I0)) + 
#   geom_density(alpha = .2)+
#   scale_x_continuous(expression(I[max]), expand = c(0, 0)) + 
#   scale_y_continuous("Density", expand = c(0, 0)) + 
#   facet_wrap(prior~.) + 
#   theme_bw(base_size = 16)
# p1

summarise <- function(x, alpha = 0.95){
  return(c(
    lwr =  quantile(x, probs = (1-alpha)/2),
    mean = mean(x),
    lwr =  quantile(x, probs = (1+alpha)/2)
  ))
}

print( aggregate(x = Imax.dt$Imax, 
                 by = list(I = Imax.dt$I0, P  = Imax.dt$prior), 
                 FUN = summarise), digits = 3)

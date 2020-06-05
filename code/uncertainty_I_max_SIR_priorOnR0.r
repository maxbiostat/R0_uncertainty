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

Imax.1.7 <- exact_I_max(R_0 = 1.7, I0 = 1.67E-7, S0 = 1-1.67E-7)
Imax.2.1 <- exact_I_max(R_0 = 2.1, I0 = 1.67E-7, S0 = 1-1.67E-7)
Imax.2.5 <- exact_I_max(R_0 = 2.5, I0 = 1.67E-7, S0 = 1-1.67E-7)

curve(exact_I_max(x, I0 = 1.67E-7, S0 = 1-1.67E-7), 1.05, 5, lwd = 3, 
      xlab = expression(R[0]), ylab = expression(I[max]))

segments(x0 = 2.5, y0 = Imax.2.5, x1 = 2.5, y1 = 0, lty = 2)
segments(x0 = 0, y0 = Imax.2.5, x1 = 2.5, y1 = Imax.2.5, lty = 2)

segments(x0 = 1.7, y0 = Imax.1.7, x1 = 1.7, y1 = 0, lty = 2, col = 2)
segments(x0 = 0, y0 = Imax.1.7, x1 = 1.7, y1 = Imax.1.7, lty = 2, col = 2)

segments(x0 = 2.1, y0 = Imax.2.1, x1 = 2.1, y1 = 0, lty = 2, col = 3)
segments(x0 = 0, y0 = Imax.2.1, x1 = 2.1, y1 = Imax.2.1, lty = 2, col = 3)

legend(x = "bottomright", legend = c("R_0 = 1.7", "R_0 = 2.1", "R_0 = 2.5"), col = c(2, 3, 1), bty = 'n', lwd = 2, lty = 2)

### Elicitation

m <- 2.5 # E[R0]
v <- .65^2 # Var(R0)

# Gamma
alpha <- m^2/v
beta <- m/v
qgamma(p = c(.025, .975), shape = alpha, rate = beta)

# Log-normal
mu <- LogMean(realMean = m, realSD = sqrt(v) )
sigma <- sqrt( LogVar(realMean = m, realSD = sqrt(v) )    ) 

qlnorm(p = c(.025, .975), meanlog = mu, sdlog = sigma)
### Comparing densities

minR <- 1.01
maxR <- 5

curve(dlnorm(x, meanlog = mu, sdlog = sigma), lwd = 3,
      ylab = "Density", xlab = expression(R[0]), minR, maxR)
curve(dgamma(x, shape = alpha, rate = beta), lwd = 3, lty = 2, add = TRUE)
legend(x = "topright", legend = c("Gamma", "Log-normal"), lty = 2:1, bty = 'n')

crazyR <- 6

pgamma(q = crazyR, shape = alpha, rate = beta)
plnorm(q = crazyR, meanlog = mu, sdlog = sigma)

## Sampling 
M <- 1E6

R0.samples.G <- rgamma(n = M, shape = alpha, rate = beta)
R0.samples.L <- rlnorm(n = M, meanlog = mu, sdlog = sigma)

N <- 1
I0a <- 1E-7
I0b <- 1E-6
I0c <- 1E-5

Imax.samples.I0a.G <- exact_I_max(R_0 = R0.samples.G, I0 = I0a, S0 = N-I0a)
Imax.samples.I0a.L <-  exact_I_max(R_0 = R0.samples.L, I0 = I0a, S0 = N-I0a)

Imax.samples.I0b.G <- exact_I_max(R_0 = R0.samples.G, I0 = I0b, S0 = N-I0b)
Imax.samples.I0b.L <-  exact_I_max(R_0 = R0.samples.L, I0 = I0b, S0 = N-I0b)

Imax.samples.I0c.G <- exact_I_max(R_0 = R0.samples.G, I0 = I0c, S0 = N-I0c)
Imax.samples.I0c.L <-  exact_I_max(R_0 = R0.samples.L, I0 = I0c, S0 = N-I0c)

Imax.dt <- data.frame(
  Imax = c(Imax.samples.I0a.G, Imax.samples.I0b.G, Imax.samples.I0c.G,
           Imax.samples.I0a.L, Imax.samples.I0b.L, Imax.samples.I0c.L),
  prior = rep(c("Gamma", "Lognormal"), each = 3*M),
  I0 = rep(c(I0a, I0b, I0c), each = M)
)
Imax.dt$I0 <- as.factor(Imax.dt$I0)

library(ggplot2)

p0 <- ggplot(data = Imax.dt, aes(x = Imax, col = prior, fill = prior)) + 
  geom_density(alpha = .2)+
  scale_x_continuous(expression(I[max]), expand = c(0, 0)) + 
  scale_y_continuous("Density", expand = c(0, 0)) + 
  geom_vline(xintercept = 1/2, linetype = "longdash") +
  facet_wrap(I0~.) + 
  theme_bw(base_size = 16)
p0

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
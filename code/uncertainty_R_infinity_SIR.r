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
get_R_infinity <- function(R_0, N, S0){
  a <- R_0/N
  b <- R_0 - log(S0)
  ## Returns R(\infty)
  ans <- -VGAM::lambertW(-a * exp(-b))/a
  return(N-ans)
}
###################################################################
###################################################################

RInf.1.7 <- get_R_infinity(R_0 = 1.7, N = 1, S0 = 1-1E-5)
RInf.2.1 <- get_R_infinity(R_0 = 2.1, N = 1, S0 = 1-1E-5)
RInf.2.5 <- get_R_infinity(R_0 = 2.5, N = 1, S0 = 1-1E-5)

curve(get_R_infinity(x, N = 1, S0 = 1-1E-5), 1.05, 5, lwd = 3, 
      xlab = expression(R[0]), ylab = expression(R(infinity)))

segments(x0 = 2.5, y0 = RInf.2.5, x1 = 2.5, y1 = 0, lty = 2)
segments(x0 = 0, y0 = RInf.2.5, x1 = 2.5, y1 = RInf.2.5, lty = 2)

segments(x0 = 1.7, y0 = RInf.1.7, x1 = 1.7, y1 = 0, lty = 2, col = 2)
segments(x0 = 0, y0 = RInf.1.7, x1 = 1.7, y1 = RInf.1.7, lty = 2, col = 2)

segments(x0 = 2.1, y0 = RInf.2.1, x1 = 2.1, y1 = 0, lty = 2, col = 3)
segments(x0 = 0, y0 = RInf.2.1, x1 = 2.1, y1 = RInf.2.1, lty = 2, col = 3)

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
I0c <- 1E-1

Rinf.samples.I0a.G <- get_R_infinity(R_0 = R0.samples.G, N = N, S0 = N-I0a)
Rinf.samples.I0a.L <-  get_R_infinity(R_0 = R0.samples.L, N = N, S0 = N-I0a)

Rinf.samples.I0b.G <- get_R_infinity(R_0 = R0.samples.G, N = N, S0 = N-I0b)
Rinf.samples.I0b.L <-  get_R_infinity(R_0 = R0.samples.L, N = N, S0 = N-I0b)

Rinf.samples.I0c.G <- get_R_infinity(R_0 = R0.samples.G, N = N, S0 = N-I0c)
Rinf.samples.I0c.L <-  get_R_infinity(R_0 = R0.samples.L, N = N, S0 = N-I0c)

Rinf.dt <- data.frame(
  Rinf = c(Rinf.samples.I0a.G, Rinf.samples.I0b.G, Rinf.samples.I0c.G,
           Rinf.samples.I0a.L, Rinf.samples.I0b.L, Rinf.samples.I0c.L),
  prior = rep(c("Gamma", "Lognormal"), each = 3*M),
  I0 = rep(c(I0a, I0b, I0c), each = M)
)
Rinf.dt$I0 <- as.factor(Rinf.dt$I0)

library(ggplot2)

p0 <- ggplot(data = Rinf.dt, aes(x = Rinf, col = prior, fill = prior)) + 
  geom_density(alpha = .2)+
  scale_x_continuous(expression(R(infinity)), expand = c(0, 0)) + 
  scale_y_continuous("Density", expand = c(0, 0)) + 
  facet_wrap(I0~.) + 
  theme_bw(base_size = 16)
p0

summarise <- function(x, alpha = 0.95){
  return(c(
    lwr =  quantile(x, probs = (1-alpha)/2),
    mean = mean(x),
    lwr =  quantile(x, probs = (1+alpha)/2)
  ))
}

print( aggregate(x = Rinf.dt$Rinf, 
                 by = list(I = Rinf.dt$I0, P  = Rinf.dt$prior), 
                 FUN = summarise), digits = 5)


## This script will show the induced densities on R0 from three choices of priors on (beta, gamma), under a constant population SIR model
###############
## Preliminaries
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
kappa <- function(m, v){
  r2 <- sqrt(2)
  d <- m^2/(2*v)
  ans <-  (r2*gsl::gamma_inc(1, d)*v + (2*sqrt(pi)-gsl::gamma_inc(1/2, d))*m*sqrt(v))/r2
  return(ans)
}
f_hn_ratio <- function(z, m1, v1, m2, v2){ ## exact density
  s1 <- sqrt(v1)
  s2 <- sqrt(v2)
  const <- 2*pi * s1 * s2* exp(pnorm(0, mean = m1, sd = s1, log.p = TRUE, lower.tail = FALSE) +
                                 pnorm(0, mean = m2, sd = s2, log.p = TRUE, lower.tail = FALSE))
  if(z==0) z <- 1e-30
  k <- exp(-(m1/z-m2)^2/(2*(v1/z^2 + v2)))
  mu <- (m1*v2*z + m2*v1)/(v2*z^2 + v1)
  v <- (v1*v2)/(v2*z^2 + v1)
  dens <- kappa(mu, v)
  return(k * dens/const)
}
f_hn_ratio <- Vectorize(f_hn_ratio)

devtools::source_url("https://raw.githubusercontent.com/maxbiostat/R0_uncertainty/master/code/gamma_ratio.R")
###############
mb <- 2
vb <- 1^2
mg <- 0.4
vg <- 0.5^2

# Gammas
alphab <- mb^2/vb
betab <- mb/vb

alphag <- mg^2/vg
betag <- mg/vg

# Log-normals
mub <- LogMean(realMean = mb, realSD = sqrt(vb) )
sigmab <- sqrt( LogVar(realMean = mb, realSD = sqrt(vb) ))

mug <- LogMean(realMean = mg, realSD = sqrt(vg) )
sigmag <- sqrt( LogVar(realMean = mg, realSD = sqrt(vg) ))

####
f1 <- function(x) dgamma.ratio(x = x, k1 = alphab, k2 = alphag, t1 = betab, t2 = betag, N = 1); f1 <- Vectorize(f1)
f2 <- function(x) f_hn_ratio(z = x, m1 = mb, v1 = vb, m2 = mg, v2 = vg); f2 <- Vectorize(f2)
f3 <- function(x) dlnorm(x = x, meanlog = mub-mug, sd = sqrt(sigmab^2 + sigmag^2)); f3 <- Vectorize(f3)


PrLessThan1.gamma <- integrate(function(x) f1(x), 0, 1)$value
PrLessThan1.halfnormal <- integrate(function(x) f2(x), 0, 1)$value
PrLessThan1.lognormal <- integrate(function(x) f3(x), 0, 1)$value

PrLessThan10.gamma <- integrate(function(x) f1(x), 0, 10)$value
PrLessThan10.halfnormal <- integrate(function(x) f2(x), 0, 10)$value
PrLessThan10.lognormal <- integrate(function(x) f3(x), 0, 10)$value

PrLessThan100.gamma <- integrate(function(x) f1(x), 0, 100)$value
PrLessThan100.halfnormal <- integrate(function(x) f2(x), 0, 100)$value
PrLessThan100.lognormal <- integrate(function(x) f3(x), 0, 100)$value

mR <- 0
MR <- 10
curve(f1, mR, MR, lwd = 4, col = 1, xlab = expression(R[0]), ylab = "Density", ylim = c(0, .2))
curve(f2, mR, MR, lwd = 4, col = 2, lty = 2, add = TRUE)
curve(f3, mR, MR, lwd = 4, col = 3, lty = 3, add = TRUE)
legend(x = "topright",
       legend = c("Gammas", "Half-normals", "Log-normals"),
       col = 1:3,
       lty = 1:3,
       lwd = 2,
       bty = 'n')

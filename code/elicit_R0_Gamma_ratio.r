#### Script to elicit the (four) parameter of a Gamma ratio given the mean and the variance (standard deviation)
GR.mean <- function(k1, k2, t1, t2){
  (k1*t1)/(t2*(k2 - 1))
}
GR.sd <- function(k1, k2, t1, t2){
  sqrt(
    (t1/t2)^2 *((k1 + k2 -1)*k1)/((k2-2)*(k2-1)^2)
  )
}
lossR0 <- function(pars, N, Mean, Var, alpha = .5){
  a1 <- pars[1]
  a2 <- pars[2]
  b1 <- pars[3]
  b2 <- pars[4]
  m <- GR.mean(k1 = a1, k2 = a2, t1 = b1*N, t2 = b2)
  M.err <- sqrt((m-Mean)^2)
  s <- GR.sd(k1 = a1, k2 = a2, t1 = b1*N, t2 = b2)
  S.err <- sqrt((s-sqrt(Var))^2)
  return(alpha*M.err + (1-alpha)*S.err)
}
Nn <- 1
Opt <- optim(par = c(1, 2.5, 1, 1), lossR0, lower = c(0, 2, 0, 0),
      method = "L-BFGS-B", Mean = 1, Var = 2^2, N = Nn, alpha = .9) 
Opt.par <- Opt$par
GR.mean(k1 = Opt.par[1], k2 = Opt.par[2], t1 = Opt.par[3] * Nn, t2 = Opt.par[4])
GR.sd(k1 = Opt.par[1], k2 = Opt.par[2], t1 = Opt.par[3] * Nn, t2 = Opt.par[4])
Opt.par
curve(dgamma(x, shape = Opt.par[1], scale = Opt.par[3]), lwd = 3, 0, 5,
      xlab = expression(beta), ylab = "Density")
curve(dgamma(x, shape = Opt.par[2], scale = Opt.par[4]), lwd = 3, 0, 5,
      xlab = expression(gamma), ylab = "Density" , col = 2)
devtools::source_url("https://raw.githubusercontent.com/maxbiostat/CODE/master/R/DISTRIBUTIONS/gamma_ratio.R")
curve(dgamma.ratio(x, k1 = Opt.par[1], k2 = Opt.par[2], t1 = Opt.par[3], t2 = Opt.par[4], N = Nn), lwd = 3, 0, 10,
      xlab = expression(gamma), ylab = "Density", col = 3)
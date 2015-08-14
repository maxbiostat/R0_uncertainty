# opt.quantile.gamma <- function(b, mean, CI, alpha){ # optimizing in terms of the obtained coverage
#   astar <- mean*b   
#   covstar <- pgamma(CI[2], astar, b) - pgamma(CI[1], astar, b) # obtained coverage
#   return(abs(alpha - covstar))
# } # DOESN'T WORK PROPERLY AS IS
opt.quantile.gamma <- function(b, mean, CI, alpha){ # minimizing the sum of absolute errors for CI bounds 
  astar <- mean*b   
  t.CI <- qgamma(c((1-alpha)/2, (1+alpha)/2), astar, b) # tentative CI
  dev.lwr <- abs(CI[1]-t.CI[1])
  dev.upr <- abs(CI[2]-t.CI[2])
  return(dev.lwr + dev.upr)
}
elicit.gamma <- function(mean, CI, alpha, M = 1E+7){
  # mean is the expected value
  # interval is a vector with upper and lower bounds 
  # alpha is the nominal coverage probability of interval
  # M is the upper bound of the search for b.opt
  b.opt <- optimise(opt.quantile.gamma, interval = c(1E-3, M), mean = mean, CI = CI, alpha = alpha)$minimum
  parms <- c(a.opt = mean*b.opt, b.opt = b.opt)
    return(parms)
}

# Testing
# a.true <- 4
# b.true <- 6
# alpha.true <- .95
# CI.true <- qgamma(c((1-alpha.true)/2, (1+alpha.true)/2), a.true, b.true)
# elicit.gamma(mean = a.true/b.true, CI = CI.true, alpha = alpha.true)

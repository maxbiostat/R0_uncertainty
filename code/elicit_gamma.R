# opt.quantile.gamma <- function(b, Mean, CI, alpha){ # optimizing in terms of the obtained coverage
#   astar <- Mean*b   
#   covstar <- pgamma(CI[2], astar, b) - pgamma(CI[1], astar, b) # obtained coverage
#   return(abs(alpha - covstar))
# } # DOESN'T WORK PROPERLY AS IS
opt.quantile.gamma <- function(b, Mean, CI, alpha){ # minimizing the sum of absolute errors for CI bounds 
  astar <- Mean*b   
  t.CI <- qgamma(c((1-alpha)/2, (1+alpha)/2), astar, b) # tentative CI
  dev.lwr <- abs(CI[1]-t.CI[1])
  dev.upr <- abs(CI[2]-t.CI[2])
  return(dev.lwr + dev.upr)
}
opt.gini.gamma <- function(b, Mean, gini){ # minimizing the difference between calculated and target Gini indexes 
  astar <- Mean*b   
  calc.gini <- gamma((2*astar + 1)/2)/(astar*gamma(astar)*sqrt(pi))
  return(abs(calc.gini-gini))
}
elicit.gamma <- function(Mean, CI, alpha, gini = NULL, M = 1E+7){
  ## Function to elicit parameters of a Gamma distribution using some prior information on mean, CI or Gini coefficient
  # mean is the expected value
  # interval is a vector with upper and lower bounds 
  # alpha is the nominal coverage probability of interval
  # M is the upper bound of the search for b.opt
  # gini is the Gini coefficient
  if(!is.null(gini)){
    b.opt <- optimise(opt.gini.gamma, interval = c(1E-3, 100), Mean = Mean, gini = gini)$minimum
  }else{
    b.opt <- optimise(opt.quantile.gamma, interval = c(1E-3, M), Mean = Mean, CI = CI, alpha = alpha)$minimum
  }
  parms <- c(a.opt = Mean*b.opt, b.opt = b.opt)
    return(parms)
}

# Testing
# a.true <- 4
# b.true <- 6
# alpha.true <- .95
# CI.true <- qgamma(c((1-alpha.true)/2, (1+alpha.true)/2), a.true, b.true)
# elicit.gamma(Mean = a.true/b.true, CI = CI.true, alpha = alpha.true)
# elicit.gamma(Mean = a.true/b.true, gini = .53) ## Brazil has a Gini index of 0.53 http://data.worldbank.org/indicator/SI.POV.GINI

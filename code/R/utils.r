pars_mm_lognormal <- function(k, theta){ ## moment-matching log-normal approximation of Gamma
  ## returns mu and sigma (standard deviation NOT variance)
  return(
    list(
      mu = log(k) + log(theta) - 0.5 * log(1 + 1/k),
      sigma = sqrt( log(1 + 1/k) )
    )
  )
}

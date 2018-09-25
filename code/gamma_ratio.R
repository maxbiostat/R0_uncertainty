## Generates density values for a gamma ratio distribution which arises in
## a certain problem of calculating the distribution of R0 in a specific SIR epidemic model. 
## In fact, when N = 1, this is distribution is simply the distribution of the ratio of two gamma variates
# Let R0 be given by \frac{\beta N}{\gamma}.
# If \beta and \gamma are Gamma-distributed with
# parameters {k1,t1} and {k2, t2}, respectively and N
# is a big positive constant, then R0  is distributed according to a gamma ratio distribution (or
# beta distribution of the second kind) with parameters k1,t1,k2,t2 and N.
# In this toy example, if we can represent our uncertainty about \beta and \gamma by means of Gamma random variables,
# we can calculate the distribution and moments analytically.
## Copyleft(or the one to blame): Carvalho, LMF (2018)
dgamma.ratio <- function(x, k1, k2, t1, t2, N, numstab = TRUE, log = FALSE){
  if(numstab){ ## "numerically stable" version
    lc <- (k1+ k2)*(log(N) + log(t1) + log(t2))-{lbeta(k1, k2) + k1*(log(N) + log(t1)) + k2*log(t2)}
    ld <- (k1-1)*log(x) - (k1+k2) * log(t2*x + N*t1)
    if(log){
      return(lc + ld)
    }else{
      return(exp(lc + ld))
    }
  }else{
    c <- {(N*t1*t2)^(k1+k2)}/(beta(k1, k2)*((N*t1)^k1)*(t2^k2)) # normalising constant
    d <- x^(k1-1) * ((t2*x + N*t1))^(-(k1+k2)) # density 
    return(c * d)
  }
}
dgamma.ratio <- Vectorize(dgamma.ratio)
##
rgamma.ratio <- function(n, k1, k2, t1, t2, N){
  A <- rgamma(n = n, shape = k1, scale = t1)*N
  B <- rgamma(n = n, shape = k2, scale = t2)
  return(A/B)
}
##
pgamma.ratio <- function(q, k1, k2, t1, t2, N, log.p = FALSE){
  integrate(dgamma.ratio, 0, q, k1 = k1, t1 = t1, k2 = k2, t2 = t2, N = N)$value
}
pgamma.ratio <- Vectorize(pgamma.ratio)
#
pgamma.ratio.exact <-  function(q, k1, k2, t1, t2, N, log.p = FALSE){
  # q <- 3
  # k1 <- 2
  # k2 <- 40
  # t1 <- 1/4000
  # t2 <- 1/200
  # N <- 1000
  lHPG <- BAS::hypergeometric2F1(k1 + k2, k1 ,  k1 + 1, -t2/(t1*N) * q)
  if(is.finite(lHPG)){
    lc <- (k1+ k2)*(log(N) + log(t1) + log(t2))-{lbeta(k1, k2) + k1*(log(N) + log(t1)) + k2*log(t2)}
    lF <- ( k1*log(q) - (k1 + k2)*log(t2*q + t1*N) + (k1 + k2)*log((t2/(t1*N))*q + 1) + lHPG) - log(k1)
    ans <- lc + lF    
  }else{
    ans <- 0
  }
  if(!log.p) ans <- exp(ans)
  return(ans)
}
pgamma.ratio.exact <- Vectorize(pgamma.ratio.exact)
##
dgamma.ratio.Ncorrected <- function(x, k1, k2, t1, t2, N, log = FALSE){
  S <- pgamma.ratio(N, k1, k2, t1, t2, N, log.p = TRUE)
  lcnew <- -S 
  d <- dgamma.ratio(x, k1, k2, t1, t2, N, log = log)
  ans <- lcnew + d
  if(!log) ans <- exp(ans)
  return(ans)
}
dgamma.ratio.Ncorrected <- Vectorize(dgamma.ratio.Ncorrected)
##
qgamma.ratio <- function(p, k1, k2, t1, t2, N, log = FALSE){# TODO: derive a proper upper limit (N doesn't work)
   optQuant <- function(x){
    (pgamma.ratio(x, k1 = k1, t1 = t1, k2 = k2, t2 = t2, N = N)-p)^2
  }
  u <- (k1/(k2-1))*(t1*N/t2)
  v <- ((N*t1/t2)^2)* (k1*(k1+k2-1))/((k2-2)*(k2-1)^2)  
  # Upr.old <- u + p*10 * sqrt(v)
  # m <- - (sqrt(v^2 + 4*u^2*v) -v - 2*u^2)/(2* u)
  # a <- (2*u^2)/{sqrt(v) * sqrt(v + 4*u^2 -v)}
  # Upr <- m/((1-p)^(1/a)) ## Pareto-based limit
  # cat(Upr.old, " ", Upr, "\n")
  Upr <- u + p*10 * sqrt(v)
  Ans <- optimise(optQuant, lower = 0, upper = Upr)# this upper limit should work
  final <- Ans$minimum
  if(log) final <- log(final)
  return(final)
}
qgamma.ratio <- Vectorize(qgamma.ratio)
##
entropy.gamma.ratio <- function(k1, t1, k2, t2, N){
  expectlog <- function(x){ -log(x) * dgamma.ratio(x, k1 = k1, t1 = t1,
                                                   k2 = k2, t2 = t2, N = N)}
  return(integrate(expectlog, 0, Inf)$value)
}
entropy.gamma.ratio <- Vectorize(entropy.gamma.ratio)
##
get_GR_mean <- function(k1, k2, t1, t2, N){
  return((k1/(k2-1))*(t1*N/t2))
}
#
get_GR_mode <- function(k1, k2, t1, t2, N){
  {(k1-1)*t1*N}/{t2*(k2+1)}
}
#
get_GR_var <- function(k1, k2, t1, t2, N){
  return (
    ((N*t1/t2)^2)* (k1*(k1+k2-1))/((k2-2)*(k2-1)^2)
  ) 
}
# 
get_GR_concavity_bound <- function(k1, k2, t1, t2, N){
  ## density is only concave for R_0 < concavity_bound
  a <- k1
  b <- t2
  c <- N*t1
  d <- k2
  return(
    (c*sqrt((a-1)*d^2+(a^2+a-2)*d+2*a^2-2*a)+(a-1)*c*d+(2*a-2)*c)/(b*d^2+3*b*d+2*b)
  )
}
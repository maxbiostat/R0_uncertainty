## Generates density values for a gamma ratio distribution which arises in
## a certain problem of calculating the distribution of R0 in a specific SIR epidemic model. 
## In fact, when N = 1, this is distribution is simply the distribution of the ratio of two gamma variates
# Let R0 be given by \frac{\beta N}{\gamma}.
# If \beta and \gamma are Gamma-distributed with
# parameters {k1,t1} and {k2,t2}, respectively and N
# is a big positive constant, then R0  is distributed according to a gamma ratio distribution (or
# beta distribution of the second kind) with parameters k1,t1,k2,t2 and N.
# In this toy example, if we can represent our uncertainty about \beta and \gamma by means of Gamma random variables,
# we can calculate the distribution and moments analytically.
## Copyleft(or the one to blame): Carvalho, LMF (2016)
dgamma.ratio <- function(x, k1, k2, t1, t2, N, numstab = TRUE){
  if(numstab){ ## "numerically stable" version
    lc <- (k1+ k2)*(log(N) + log(t1) + log(t2))-{lbeta(k1, k2) + k1*(log(N) + log(t1)) + k2*log(t2)}
    ld <- (k1-1)*log(x) - (k1+k2) * log(t2*x + N*t1)
    return(exp(lc + ld))
  }else{
    c <- {(N*t1*t2)^(k1+k2)}/(beta(k1,k2)*((N*t1)^k1)*(t2^k2)) # normalising constant
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
pgamma.ratio <- function(q, k1, k2, t1, t2, N){
  integrate(dgamma.ratio, 0, q, k1 = k1, t1 = t1, k2 = k2, t2 = t2, N = N)$value
}
pgamma.ratio <- Vectorize(pgamma.ratio)
##
qgamma.ratio <- function(p, k1, k2, t1, t2, N){# TODO: derive a proper upper limit (N doesn't work)
  optQuant <- function(x){
    (pgamma.ratio(x, k1 = k1, t1 = t1, k2 = k2, t2 = t2, N = N)-p)^2
  }
  Ans <- nlminb(start = 1, objective = optQuant, lower = 0)
  return(Ans$par)
}
qgamma.ratio <- Vectorize(qgamma.ratio)
##
dgamma.ratio.Ncorrected <- function(x, k1, k2, t1, t2, N){
  c <- {(N*t1*t2)^(k1+k2)}/(beta(k1, k2)*((N*t1)^k1)*(t2^k2)) 
  S <- integrate(dgamma.ratio, k1 = k1, k2 = k2, t1 = t1, t2 = t2, N = N, 0, N)$value
  cnew <- (1/S) * c # regularized normalising constant
  d <- x^(k1-1) * ((t2*x + N*t1))^(-(k1+k2)) # density 
  return(cnew * d)
}
##
entropy.gamma.ratio <- function(k1, t1, k2, t2, N){
  expectlog <- function(x){ -log(x) * dgamma.ratio(x, k1 = k1, t1 = t1,
                                                   k2 = k2, t2 = t2, N = N)}
  return(integrate(expectlog, 0, Inf)$value)
}
entropy.gamma.ratio <- Vectorize(entropy.gamma.ratio)
# mean =  (k1/(k2-1))*(t1*N/t2)
# mode = {(k1-1)*t1*N}/{t2*(k2+1)}
# var = ((N*t1/t2)^2)* (k1*(k1+k2-1))/((k2-2)*(k2-1)^2)
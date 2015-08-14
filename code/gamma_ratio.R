## Generates density values for a gamma ratio distribution which arises in
## a certain problem of calculating the distribution of R0 in a specific SIR epidemic model 
## In fact, when N = 1, this is distribution is simply the distribution of the ratio of two gamma variates
# Let R0 be given by \frac{\beta N}{\gamma}.
# If \beta and \gamma are Gamma-distributed with # parameters {k1,t1} and {k2,t2}, respectively  and N
# is a big positive constant, then R0  is distributed according to a gamma ratio distribution beta distribution of the second kind)
# with parameters k1,t1,k2,t2 and N.
# In this toy example, if we can represent our uncertainty about \beta and \gamma by means of Gamma random variables,
# we can calculate the distribution and moments analytically.
## Copyleft(or the one to blame): Carvalho, LMF (2014)

dgamma.ratio <- function(x, k1, k2, t1, t2, N){
  c <- {(N*t1*t2)^(k1+k2)}/(beta(k1,k2)*((N*t1)^k1)*(t2^k2)) # normalizing constant
  d <- x^(k1-1) * ((t2*x + N*t1))^(-(k1+k2)) # density 
  return(c * d)
}

dgamma.ratio.Ncorrected<- function(x, k1, k2, t1, t2, N){
  c <- {(N*t1*t2)^(k1+k2)}/(beta(k1, k2)*((N*t1)^k1)*(t2^k2)) 
  S <- integrate(dgamma.ratio, k1 = k1, k2 = k2, t1 = t1, t2 = t2, N = N, 0, N)$value
  cnew <- (1/S) * c # regularized normalizing constant
  d <- x^(k1-1) * ((t2*x + N*t1))^(-(k1+k2)) # density 
  return(cnew * d)
}

entropy.gamma.ratio <- function(k1, t1, k2, t2, N){
  expectlog <- function(x){ -log(x) * dgamma.ratio(x, k1 = k1, t1 = t1,
                                                   k2 = k2, t2 = t2, N = N)}
  return(integrate(expectlog, 0, Inf)$value)
}
# mean =  (k1/(k2-1))*(t1*N/t2)
# mode = {(k1-1)*t1*N}/{t2*(k2+1)}
# var = ((N*t1/t2)^2)* (k1*(k1+k2-1))/((k2-2)*(k2-1)^2)
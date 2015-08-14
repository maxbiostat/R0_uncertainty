entropy.gamma.ratio <- function(k1, t1, k2, t2, N){
  expectlog <- function(x){ -log(x) * dgamma.ratio(x, k1 = k1, t1 = t1,
                                                   k2 = k2, t2 = t2, N = N)}
  return(integrate(expectlog, 0, Inf)$value)
}
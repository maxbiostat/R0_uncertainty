### This piece of code implements the KL divergence minimization for transforms
# The idea is as follows: 
# Let y = M(x), where both x and y can be multivariate and M() is a continuous function.
# Imagine each expert gives a distribution for x, f_i(x), without taking into account the induced distribution on y, g_i(y)
# If one pools the distributions for x using L(F(x); alpha) one obtains \pi(x), which in turn induces a distribution q_1(y)
# At this point we might want to choose alpha such that the divergence between each induced distribution g_i(y) and q_1(y)
# is minimized.
# The problem is then pick an \alpha such that L = sum(KL(g_i(y)||q_1(y))) is minimum
source("elicit_gamma.R")
source("gamma_ratio.R")
source("../maxent_aux.R")
## Four distributions for x = (\beta, \gamma)
# N's
N0 <- 1000
N1 <- 2000
N2 <- 1500
N3 <- 2000
Nv <- c(N0, N1, N2, N3)
# Parameters for the distributions of \beta
k10 <- 2 ; t10 <- 1/4000 
k11 <- 1.5; t11 <- 1/4000 
k12 <- 2; t12 <- 1/4500
k13 <- 4; t13 <- 1/8000 
k1v <- c(k10, k11, k12, k13)
t1v <- c(t10, t11, t12, t13)
# Same for \gamma
k20 <- 40; t20 <- 1/200
k21 <- 50 ; t21 <- 1/200
k22 <- 15; t22 <- 1/50
k23 <- 30; t23 <- 1/300
k2v <- c(k20, k21, k22, k23)
t2v <- c(t20, t21, t22, t23)
########################################
########################################
optkltransform <- function(alpha, a1p, b1p, a2p, b2p, Np){
  # Let's first compute q_1(y) = M(\pi(x))
  K <- length(alpha)
  ds <- rep(NA, K) # the distances from each f_i to \pi
  for (i in 1:K){
    ds[i] <-  kl.transform(alpha = alpha, a1v = a1p, b1v = b1p,
                           a2v = a2p, b2v = b2p, Nv = Np, 
                           ga1= a1p[i], gb1 = b1p[i],
                           ga2 = a2p[i], gb2 = b2p[i], gN = Np[i])
  } 
  return(ds)
}

optkltransform.inv <- function(alpha.inv, a1p, b1p, a2p, b2p, Np){
  alpha <- alpha.01(alpha.inv)
  sum(optkltransform(alpha, a1p, b1p, a2p, b2p, Np))
}

a <- optim(c(0, 0, 0), optkltransform.inv, a1p = k1v, b1p = t1v, a2p = k2v, b2p = t2v, Np = Nv)
#            method = "SANN", control=list(maxit = 100000))
(round(alpha.opt <- alpha.01(a$par), 2))
########################################
########################################
(k1star <- sum(alpha.opt*k1v)) 
(t1star <- 1/sum(alpha.opt*t1v))
(k2star <- sum(alpha.opt*k2v))
(t2star <- 1/sum(alpha.opt*t2v))
(Nstar <-  sum(alpha.opt*Nv))
## The pooled distribution for x[1] = \beta
curve(dgamma(x, k1star, t1star), 0, 3E-3, main = expression("Pooled distribution for the transmission rate", beta),
      ylab = expression(pi(beta)), xlab = expression(beta), col = 2, lwd = 3)
for (i in 1:length(Nv)){
  curve(dgamma(x, k1v[i], 1/(t1v[i]) ), 0, 3E-3, lty = i+1, lwd = 2, add = TRUE)  
}

## Now \pi(x[2]) = \pi(\gamma)
curve(dgamma(x, k2star, t2star), 0, .5, main = expression("Pooled distribution for the recovery/removal rate", gamma),
      ylab = expression(pi(gamma)), xlab = expression(gamma), col = 3, lwd = 2)
for (i in 1:length(Nv)){
  curve(dgamma(x, k2v[i], 1/(t2v[i]) ), 0, .5, lty = i+1, lwd = 2, add = TRUE)  
}  
## And finally what we got for R_0

curve(dgamma.ratio(x, k1 = k1star, t1 = 1/t1star, k2 = k2star, t2 = 1/t2star, N = Nstar), 0, 15,
      main = expression("Pooled distribution for the reproductive number", R[0]),
      ylab = expression(pi(R[0])), xlab = expression(R[0]), col = 4, lwd = 2)
for (i in 1:length(Nv)){
  curve(dgamma.ratio(x, k1 = k1v[i], t1 = t1v[i], k2 = k2v[i], t2 = t2v[i], N = Nv[i]), 0, 15,
        lty = i+1, lwd = 2, add = TRUE)  
}


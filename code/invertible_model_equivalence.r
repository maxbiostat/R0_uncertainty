## In this script we will investigate the result that says that if the model M(x) = y is invertible,
## then the order of pooling and inducing does not matter.
## We will (re-)use some of the information of the example from avchuk & Martz, 1994 (IEEE)
## Suppose f_i(x) = dbeta(x, a_i, b_i) and M(x) = x/(x+1)
source("maxent_aux.R")
#### Parameters for the Beta distributions
a0 <- 18.1 ; b0 <- .995
a1 <- 3.44 ; b1 <- .860 
a2 <- 8.32 ; b2 <- .924
a3 <- 1.98 ; b3 <- .848

av <- c(a0, a1, a2, a3)
bv <- c(b0, b1, b2, b3)
K <- length(av)

#### Other auxiliary functions
M <- function(x) x/(x + 1)
Minv <- function(y) -y/(y-1) ## M_inverse
J <- function(y) 1/(y-1)^2 ## Jacobian
### Pool-then-induce
PtI.dens <- function(y, va, vb, alphas){
  parms <- pool.par(alpha = alphas, a = va, b = vb)
  return(dbeta(Minv(y), shape1 = parms[1], shape2 = parms[2])*J(y))
} 
## Induce-then-pool
ItP.dens <- function(y, va, vb, alphas){
  Ds <- list(
    f0 = function(y){dbeta(Minv(y), shape1 = va[1], shape2 = vb[1])*J(y)},
    f1 = function(y){dbeta(Minv(y), shape1 = va[2], shape2 = vb[2])*J(y)},
    f2 = function(y){dbeta(Minv(y), shape1 = va[3], shape2 = vb[3])*J(y)},
    f3 = function(y){dbeta(Minv(y), shape1 = va[4], shape2 = vb[4])*J(y)}
  ) # list with the densities
  return(dpoolnorm.unit(y, D = Ds, alphas = alphas))
}
####
# Generate from PtI.dens (mind you that if the remark is true this doesn't matter)
# 
smp <- runif(n = K, min = 1, max = 100)
(alph <- smp/sum(smp) )# rep(1/K, K)
parameters <-  pool.par(alpha = alph, a = av, b = bv)
X <- rbeta(100000, shape1 = parameters[1], shape2 = parameters[2])
hist(X)
Y <- M(X)
dPtI <- function(y) PtI.dens(y, va = av, vb = bv, alphas = alph)
dItP <- function(y) ItP.dens(y, va = av, vb = bv, alphas = alph)

hist(Y, probability = TRUE)
curve(dPtI, min(Y), max(Y), lwd = 2, col = "green", add = TRUE)
curve(dItP, min(Y), max(Y), lwd = 2, lty = 3, col = "red", add = TRUE)
legend(x = "topleft", legend = c("Pool-then-induce", "Induce-then-pool"),
       lwd = 2, col = c(3, 2), lty = c(1,3), bty = "n")
dItP(.3)
dItP(.3)

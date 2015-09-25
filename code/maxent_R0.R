source("elicit_gamma.R")
source("gamma_ratio.R")
source("../../CODE/maxent_aux.R")
source("R0Example_four_gammas_parameters.r")
#################
K <- length(k2v)

###############
# First approach: combine each beta and each gamma; construct pooled R0
# This is the Pool-then-Induce (PI) approach
###############

optentgammaRPI <- function(alpha, a1p, b1p, a2p, b2p, Np){
  entropy.gamma.ratio(k1 = sum(alpha*a1p), t1 = sum(alpha*b1p), k2 = sum(alpha*a2p),
                       t2 = sum(alpha*b2p), N = sum(alpha*Np))
}
optentgammaRPI.inv <- function(alpha.inv, a1p, b1p, a2p, b2p, Np){
  alpha <- alpha.01(alpha.inv)
  -optentgammaRPI(alpha = alpha, a1p, b1p, a2p, b2p, Np)
}
aPI <- optim(rep(0, K-1), optentgammaRPI.inv, a1p = k1v, b1p = t1v, a2p = k2v, b2p = t2v, Np = Nv) 
round(alphasPI <- alpha.01(aPI$par), 3)

# BetaN
(k1star <- as.numeric(crossprod(k1v, alphasPI))) 
(t1Nstar <- as.numeric(crossprod(t1v*Nv, alphasPI))) # theta_1i*N_i

# Gamma
(k2star <- as.numeric(crossprod(k2v, alphasPI)))
(t2star <- as.numeric(crossprod(t2v, alphasPI)))

entropy.gamma.ratio(t1 = t1Nstar, t2 = t2star, k1 = k1star, k2 = k2star, N = 1)

###############
# Second approach: combine each induced distribution; pool these
# This is the Induce-then-Pool (IP) approach
###############

Ds <- list(
  f0 = function(x){dgamma.ratio(x, t1 = t10 , t2 = t20, k1 = k10, k2 = k20 , N = N0)},
  f1 = function(x){dgamma.ratio(x, t1 = t11 , t2 = t21, k1 = k11, k2 = k21 , N = N1)},
  f2 = function(x){dgamma.ratio(x, t1 = t12 , t2 = t22, k1 = k12, k2 = k22 , N = N2)},
  f3 = function(x){dgamma.ratio(x, t1 = t13 , t2 = t23, k1 = k13, k2 = k23 , N = N3)}
) # list with the densities


optentgammaRIP <- function(alpha, D){
  dpool.entropy(alphas = alpha, D = D)
}

optentgammaRIP.inv <- function(alpha.inv, D){
  alpha <- alpha.01(alpha.inv)
  -optentgammaRIP(alpha = alpha, D = D)
}

aIP <- optim(rep(0, K-1), optentgammaRIP.inv, D = Ds) 
round(alphasIP <- alpha.01(aIP$par), 3)

###############
# Plotting
###############
pdf("../figures/ItP_vs_PtI_maxEnt.pdf")
curve(dgamma.ratio(x, t1 = t1Nstar, t2 = t2star, k1 = k1star, k2 = k2star, N = 1 ),
      0, 10, main = "Pooled distributions", xlab = expression(R[0]),
      ylab = "Density", lwd = 2)
curve(dpoolnorm.positive(x, D = Ds, alpha = alphasIP), 0, 10, lwd = 2,
      lty = 2, add = TRUE, col = 2)
legend(x = "topright", bty = "n", legend = c("Pool-then-induce", "Induce-then-pool"),
       col = 1:2, lty = 1:2, lwd = 2)
dev.off()
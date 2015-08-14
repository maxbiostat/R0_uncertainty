source("elicit_gamma.R")
source("gamma_ratio.R")
source("../maxent_aux.R")
# Four distributions
# N's
N0 <- 1000
N1 <- 2000
N2 <- 1500
N3 <- 2000
Nv <- c(N0, N1, N2, N3)
# Betas
k10 <- 2 ; t10 <- 1/4000 
k11 <- 1.5; t11 <- 1/4000 
k12 <- 2; t12 <- 1/4500
k13 <- 4; t13 <- 1/8000 
k1v <- c(k10, k11, k12, k13)
t1v <- c(t10, t11, t12, t13)

# Gammas
k20 <- 40; t20 <- 1/200
k21 <- 50 ; t21 <- 1/200
k22 <- 15; t22 <- 1/50
k23 <- 30; t23 <- 1/300
k2v <- c(k20, k21, k22, k23)
t2v <- c(t20, t21, t22, t23)

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
round(alphasPI <- alpha.01(aPI$par), 2)

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
round(alphasIP <- alpha.01(aIP$par), 2)

###############
# Plotting
###############
curve(dgamma.ratio(x, t1 = t1Nstar, t2 = t2star, k1 = k1star, k2 = k2star, N = 1 ),
      0, 10, main = "Pooled distributions", xlab = expression(R[0]),
      ylab = "Density", lwd = 2)
curve(dpoolnorm.positive(x, D = Ds, alpha = alphasIP), 0, 10, lwd = 2,
      lty = 2, add = TRUE, col = 2)

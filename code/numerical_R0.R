source("elicit_gamma.R")
source("gamma_ratio.R")
source("../../CODE/maxent_aux.R")
source("R0Example_four_gammas_parameters.r")
#################
K <- length(k2v)
alphas <- rep(1, K)/K # equal weights
###############
# First approach: combine each beta and each gamma; construct pooled R0
###############

# BetaN
(k1star <- as.numeric(crossprod(k1v, alphas))) 
(t1Nstar <- as.numeric(crossprod(t1v*Nv, alphas))) # theta_1i*N_i

# Gamma
(k2star <- as.numeric(crossprod(k2v, alphas)))
(t2star <- as.numeric(crossprod(t2v, alphas)))

###############
# Second approach: combine each induced distribution; pool these
###############
Ds <- list(
  f0 = function(x){dgamma.ratio(x, t1 = t10 , t2 = t20, k1 = k10, k2 = k20 , N = N0)},
  f1 = function(x){dgamma.ratio(x, t1 = t11 , t2 = t21, k1 = k11, k2 = k21 , N = N1)},
  f2 = function(x){dgamma.ratio(x, t1 = t12 , t2 = t22, k1 = k12, k2 = k22 , N = N2)},
  f3 = function(x){dgamma.ratio(x, t1 = t13 , t2 = t23, k1 = k13, k2 = k23 , N = N3)}
) # list with the densities

###############
# Plotting
###############
pdf("../figures/ItP_vs_PtI_equalWeights.pdf")
curve(dgamma.ratio(x, t1 = t1Nstar, t2 = t2star, k1 = k1star, k2 = k2star, N = 1 ),
      0, 15, main = "Pooled distributions", xlab = expression(R[0]),
      ylab = "Density", lwd = 2)
curve(dpoolnorm.positive(x, D = Ds, alpha = alphas), 0, 15, lwd = 2,
      lty = 2, add = TRUE, col = 2)
legend(x = "topright", bty = "n", legend = c("Pool-then-induce", "Induce-then-pool"),
       col = 1:2, lty = 1:2, lwd = 2)
dev.off()
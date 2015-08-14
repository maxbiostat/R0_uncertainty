## This script implements the analysis of the data in WHO Response team, 2014 (NEJM) on the 
## basic reproductive number (R0) of Ebola in three West African countries
source("../maxent_aux.R")
source("elicit_gamma.R")

# Data taken from the http://www.nejm.org/doi/full/10.1056/NEJMoa1411100
Ebola.data.R0 <- data.frame( matrix(NA, ncol = 3 , nrow = 3))
rownames(Ebola.data.R0) <- c("mean", "lwr", "upr")
colnames(Ebola.data.R0) <- c("Guinea", "Sierra Leone", "Liberia")   

Ebola.data.R0[, 1] <- c(1.71, 1.44, 2.01)
Ebola.data.R0[, 2] <- c(2.02, 1.79, 2.26)
Ebola.data.R0[, 3] <- c(1.83, 1.72, 1.94)

(pars <- as.matrix(apply(Ebola.data.R0, 2, function(c) elicit.gamma(mean = c[1], CI = c[2:3], alpha = .95))) )


curve(fgamma(x, par = pars[,1]), 1, 3, ylab = "Density", main = "Approximate Gamma distributions",
      xlab = expression(beta), lwd = 2, col = 2)
curve(fgamma(x, par = pars[,2]), 1, 3, lwd = 2, col = 3, add = TRUE)
curve(fgamma(x, par = pars[,3]), 1, 3, lwd = 2, col = 4, add = TRUE)
abline( v = Ebola.data.R0[,1], lty = 2, lwd = 1, col = 2)
abline( v = Ebola.data.R0[,2], lty = 2, lwd = 1, col = 3)
abline( v = Ebola.data.R0[,3], lty = 2, lwd = 1, col = 4)
legend(x = "topright", bty = "n", col = 2:4,
       legend = colnames(Ebola.data.R0), pch = 16)


apply(pars, 2, function(x) qgamma(c(.025, .975), x[1], x[2]) )


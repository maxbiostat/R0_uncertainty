source("../code/gamma_ratio.R")
par.N <- 1000
par.k1 <- 2
par.t1 <- 1/6000
par.k2 <- 40
par.t2 <- 1/200
NN <- 1E5
X <- rgamma.ratio(n = NN, k1 = par.k1, t1 = par.t1, 
                  t2 = par.t2, k2 = par.k2, N = par.N)
hist(X, probability = TRUE)
curve(dgamma.ratio(x, k1 = par.k1, t1 = par.t1, 
                   t2 = par.t2, k2 = par.k2, N = par.N), 0, max(X), add = TRUE, col = "blue", lwd = 3)
curve(dnorm(x, mean = (par.k1/(par.k2-1))*(par.t1*par.N/par.t2),
            sd = sqrt (((par.N*par.t1/par.t2)^2)* (par.k1*(par.k1+par.k2-1))/((par.k2-2)*(par.k2-1)^2))),
      0, max(X), add = TRUE,
      col = "black", lwd = 3)
legend(x = "topright", legend = c("GR", "Normal approx"),
       col = c("blue", "black"), lwd = 2, bty = "n")
median(X)
summary(X)
pgamma.ratio(median(X), k1 = par.k1, t1 = par.t1, 
             t2 = par.t2, k2 = par.k2, N = par.N)
quantile(X, .025)
qgamma.ratio(.025, k1 = par.k1, t1 = par.t1, 
             t2 = par.t2, k2 = par.k2, N = par.N)
qgamma.ratio(.5, k1 = par.k1, t1 = par.t1, 
       t2 = par.t2, k2 = par.k2, N = par.N)
quantile(X, .95)
qgamma.ratio(.95, k1 = par.k1, t1 = par.t1, t2 = par.t2, k2 = par.k2, N = par.N)

quantile(X, .999)
qgamma.ratio(.999, k1 = par.k1, t1 = par.t1, t2 = par.t2, k2 = par.k2, N = par.N)
## 

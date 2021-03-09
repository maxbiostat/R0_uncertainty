source("gamma_ratio.R")
N <- 1000
## Informative priors
# Gamma
inf.k.beta <- 2 ; inf.t.beta <- 1/4000 
inf.k.gamma <- 40; inf.t.gamma <- 1/200
# Log-normal
inf.mu.beta <- log( (inf.k.beta * inf.t.beta)/sqrt(1 + 1/inf.k.beta)) ; inf.sigma.beta <- log(1 + 1/inf.k.beta)
inf.mu.gamma <- log( (inf.k.gamma * inf.t.gamma)/sqrt(1 + 1/inf.k.gamma)) ; inf.sigma.gamma <- log(1 + 1/inf.k.gamma)
## "Non" informative priors
# Gamma
non.k.beta <- 1 ; non.t.beta <- 100 
non.k.gamma <- 1; non.t.gamma <- 100
# Log-normal
non.mu.beta <- log( (non.k.beta * non.t.beta)/sqrt(1 + 1/non.k.beta)) ; non.sigma.beta <- log(1 + 1/non.k.beta)
non.mu.gamma <- log( (non.k.gamma * non.t.gamma)/sqrt(1 + 1/non.k.gamma)) ; non.sigma.gamma <- log(1 + 1/non.k.gamma)

M <- 1e6
##
g.beta.inf <- rgamma(n = M, shape = inf.k.beta, scale = inf.t.beta)
g.gamma.inf <- rgamma(n = M, shape = inf.k.gamma, scale = inf.t.gamma)
l.beta.inf <- rlnorm(n = M, meanlog = inf.mu.beta, sdlog = sqrt(inf.sigma.beta))
l.gamma.inf <- rlnorm(n = M, meanlog = inf.mu.gamma, sdlog = sqrt(inf.sigma.gamma))

# g.beta.non <- rgamma(n = M, shape = non.k.beta, scale = non.t.beta)
# g.gamma.non <- rgamma(n = M, shape = non.k.gamma, scale = non.t.gamma)
# l.beta.non <- rlnorm(n = M, meanlog = non.mu.beta, sdlog = sqrt(non.sigma.beta))
# l.gamma.non <- rlnorm(n = M, meanlog = non.mu.gamma, sdlog = sqrt(non.sigma.gamma))

## beta (transmission rate)
plot_beta <- function(){
    curve(dgamma(x, shape = inf.k.beta, scale = inf.t.beta), 0, max(g.beta.inf),
          ylim = c(0, 2000), cex.axis = 1.0, cex.names = 1.5, cex.lab = 1.5,
          xlab = expression(beta), ylab = expression(f[beta](beta)), lwd = 2)
    curve(dgamma(x, shape = non.k.beta, scale = non.t.beta),
          lty = 2, lwd = 2, add = TRUE)
    curve(dlnorm(x, meanlog = inf.mu.beta, sdlog = sqrt(inf.sigma.beta)),
          lwd = 2, col = "grey50", add = TRUE)
    curve(dlnorm(x, meanlog = non.mu.beta, sdlog = sqrt(non.sigma.beta)),
          lwd = 2, lty = 2, col = "grey50", add = TRUE)
    legend(x = "topright", legend = c("Gamma", "Log-normal"),
           col = c("black", "grey50"), pch = 20, bty = "n", cex = 1.2)
    legend(x = "bottomright", legend = c("Informative", "Non-informative"),
           lty = c(1, 2), lwd = 2, bty = "n", cex = 1.2)
}
plot_beta()


## gamma (recovery rate)
plot_gamma <- function(){
  curve(dgamma(x, shape = inf.k.gamma, scale = inf.t.gamma), 0, max(g.gamma.inf),
        cex.axis = 1.0, cex.names = 1.5, cex.lab = 1.5,
        xlab = expression(gamma), ylab = expression(f[gamma](gamma)), lwd = 2)
  curve(dgamma(x, shape = non.k.gamma, scale = non.t.gamma),
        lty = 2, lwd = 2, add = TRUE)
  curve(dlnorm(x, meanlog = inf.mu.gamma, sdlog = sqrt(inf.sigma.gamma)),
        lwd = 2, col = "grey50", add = TRUE)
  curve(dlnorm(x, meanlog = non.mu.gamma, sdlog = sqrt(non.sigma.gamma)),
        lwd = 2, lty = 2, col = "grey50", add = TRUE)
}
plot_gamma()

## R0 = beta* N / gamma
plot_R0 <- function(){
  curve(dgamma.ratio(x, k1 = inf.k.beta, k2 = inf.k.gamma,
                     t1 = inf.t.beta, t2 = inf.t.gamma, N = N),
        0, 15, ylim = c(0, .35), cex.axis = 1.0, cex.names = 1.5, cex.lab = 1.5, 
        xlab = expression(R[0]), ylab = expression(f[R[0]](R[0])), lwd = 2)
  # abline(v = get_GR_concavity_bound(k1 = inf.k.beta, k2 = inf.k.gamma,
  #                                   t1 = inf.t.beta, t2 = inf.t.gamma, N = N), lty = 3, lwd = 2)
  curve(dlnorm(x, meanlog = log(N) + inf.mu.beta-inf.mu.gamma,
               sdlog = sqrt(inf.sigma.beta + inf.sigma.gamma)),
        lwd = 2, col = "grey50", add = TRUE)
  curve(dgamma.ratio(x, k1 = non.k.beta, k2 = non.k.gamma,
                     t1 = non.t.beta, t2 = non.t.gamma, N = N), lty = 2, lwd = 2, add = TRUE)
  curve(dlnorm(x, meanlog = log(N) + non.mu.beta-non.mu.gamma,
               sdlog = sqrt(non.sigma.beta + non.sigma.gamma)),
        lwd = 2, col = "grey50", lty = 2, add = TRUE)
}
plot_R0()


# pdf("../figures/informative_example.pdf")
par(mfrow = c(2, 2), mgp = c(2.5, 1, 0), mar = c(5, 4, 4 , 2) + 0.1)
  plot_beta()
  plot_gamma()
  plot_R0()
# dev.off()

get_LN_concavity_bound <- function(mu, sigma){
  v <- sigma^2
  exp( c(-1, 1)*  sqrt(v^2 +4*v)/2-(3*v)/2 + mu)
}
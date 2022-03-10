functions {
           vector SI(real t,
                     vector y,
                     vector params) { 
                                      vector[3] dydt;  
                                      dydt[1] = - params[1] * y[1] * y[2];
                                      dydt[2] = params[1] * y[1] * y[2] - params[2] * y[2];
                                      dydt[3] = params[2] * y[2];
                                      return dydt;
                                    }
}
data {
      int<lower = 0> n_obs; // Number of days sampled
      array[n_obs] real ts; // Time points that were sampled
      real t0; // initial time
      real y_init; // initial measured population
      vector[n_obs] y; // measured population at measurement times
      real<lower=0> kb;
      real<lower=0> thetab;
      real<lower=0> kg;
      real<lower=0> thetag;
      real<lower=0> as;
      real<lower=0> bs;
}
transformed data{
                 real mu_b = log(kb / thetab) - (0.5 * log(1 + 1 / kb));
                 real sigma_b = sqrt(log(1 + 1 / kb));
                 real mu_g = log(kg / thetag) - (0.5 * log(1 + 1 / kg));
                 real sigma_g = sqrt(log(1 + 1 / kg));
}

parameters{
           array[2] real<lower=0, upper=1> r_init; // initial population
           real<lower=0> beta;
           real<lower=0> gamma;
           real<lower=0, upper=1> S0;
           real<lower=0, upper=1> sigma;   // error scale
}
transformed parameters {
                        real<lower=0> R0 = beta/gamma;
                        vector[3] y0; // Initial conditions for both S and I
                        y0[1] = S0;
                        y0[2] = 1 - S0;
                        y0[3] = 0;
                        vector[2] theta;
                        theta[1] = beta;
                        theta[2] = gamma;
                        array[n_obs] vector[3] y_hat = ode_rk45(SI, y0, t0, ts, theta);
  
}
model {
       target += lognormal_lpdf(beta | mu_b, sigma_b);
       target += lognormal_lpdf(gamma | mu_g, sigma_g);
       target += beta_lpdf(S0 | as, bs);
       target += normal_lpdf(sigma | 0.5, 0.5);
       target += lognormal_lpdf(y  | log(y_hat[, 2]), sigma);
}

generated quantities{
                     array[n_obs] real<lower=0> y_rep;
                     for (i in 1:n_obs) y_rep[i] = lognormal_rng(log(y_hat[i, 2]), sigma);
}

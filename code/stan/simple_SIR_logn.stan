functions {
  real[] dr_dt(real t,
               real[] r,
               real[] theta,
               real[] x_r,
               int[] x_i) {    
    real beta = theta[1];
    real gamma = theta[2];
    real s0 = theta[3];
    real R0 = beta/gamma;
    real deriv [2];
    deriv[1] = gamma*(1-r[1]-s0*exp(-R0*r[1]));
    deriv[2] = 1;
    return deriv;
  }
}
data {
  int<lower = 0> n_obs; // Number of days sampled
  real ts[n_obs]; // Time points that were sampled
  real y_init; // initial measured population
  real y[n_obs]; // measured population at measurement times
  real<lower=0> ab;
  real<lower=0> bb;
  real<lower=0> ag;
  real<lower=0> bg;
  real<lower=0> as;
  real<lower=0> bs;
}
transformed data{
  real xr[2] = rep_array(0.0, 2);
  int xi[2] = rep_array(0, 2);
}
parameters {
  real<lower=0> beta;
  real<lower=0> gamma;
  real<lower=0, upper=1> s0;
  real<lower=0> sigma;   // error scale
}
transformed parameters {
    real theta[3];   // theta = { beta, gamma, s0 }
    real  r[n_obs, 2];
    real mu_b = (3 / 2) * (log(ag / bg)) - log((1 / bg) + 1);
    real sigma_b = sqrt(log(((1 / bb) + 1) * bb / ab));
    real mu_g = (3 / 2) * (log(ab / bb)) - log((1 / bb) + 1);
    real sigma_g = sqrt(log(((1 / bg) + 1) * bg / ag));
    theta[1] = beta;
    theta[2] = gamma;
    theta[3] = s0;
    r = integrate_ode_rk45(dr_dt, {y_init, 0.0}, 0, ts, theta,
                                         xr, xi,
                                         1e-6, 1e-5, 1e3);
}
model {
  beta ~ lognormal(mu_b, sigma_b);
  gamma ~ lognormal(mu_g, sigma_g);
  s0 ~ beta(as, bs);
  sigma ~ normal(0, 1);
  y ~ lognormal(log(r[, 1]), sigma);
}
generated quantities{
  real<lower=0> R0 = beta/gamma;
  real<lower=0> y_rep [n_obs];
  for (i in 1:n_obs) y_rep[i] = lognormal_rng(log(r[i, 1]), sigma);
}


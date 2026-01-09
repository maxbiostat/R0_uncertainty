functions {
  vector SI(real t, vector y, vector params) {
    vector[3] dydt;
    dydt[1] = - params[1] * y[1] * y[2];                  // dS/dt
    dydt[2] = params[1] * y[1] * y[2] - params[2] * y[2]; // dI/dt
    dydt[3] = params[2] * y[2];                           // dR/dt
    return dydt;
  }
}

data {
  int<lower = 0> n_obs;          // number of observations
  array[n_obs] real ts;          // observation times
  real t0;                       // initial time
  real<lower=0, upper=1> y_init; // initial infected fraction I(0)
  vector[n_obs] y;               // observed infected fraction at ts
  real<lower=0> thetab; // scale for beta
  real<lower=0> thetag; // scale for gamma
  real<lower=0> as;
  real<lower=0> bs;
}

parameters {
  real<lower=0> beta;
  real<lower=0> gamma;
  real<lower=0, upper=1> S0_raw;
  real<lower=0, upper=1> sigma;  // observation noise (lognormal scale)
}

transformed parameters {
  real<lower=0> R0 = beta / gamma;
  real<lower=0, upper=(1 - y_init)> S0 = (1 - y_init) * S0_raw;

  vector[3] y0;    // initial state: [S(0), I(0), R(0)]
  y0[1] = S0;
  y0[2] = y_init;
  y0[3] = 1 - S0 - y_init;
  vector[2] theta;
  theta[1] = beta;
  theta[2] = gamma;
  array[n_obs] vector[3] y_hat = ode_rk45(SI, y0, t0, ts, theta);
}

model {
  // half-normal (priors truncationg lower bounds)
  target += normal_lpdf(beta  | 0, thetab);
  target += normal_lpdf(gamma | 0, thetag);
  target += beta_lpdf(S0_raw | as, bs);
  target += normal_lpdf(sigma | 0.5, 0.5);
  target += lognormal_lpdf(y | log(y_hat[, 2]), sigma);
}

generated quantities {
  array[n_obs] real<lower=0> y_rep;
  for (i in 1:n_obs)
    y_rep[i] = lognormal_rng(log(y_hat[i, 2]), sigma);
}

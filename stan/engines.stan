data {
  // number of observations
  int<lower = 0> N;
  // observed data
  int<lower = 0, upper = N> N_obs;
  vector[N_obs] y_obs;
  // censored data
  int<lower = 0, upper = N> N_cens;
  vector[N_cens] y_cens;
  // covariates
  int<lower = 0> K;
  matrix[N_obs, K] X_obs;
  matrix[N_cens, K] X_cens;
  // priors
  real alpha_loc;
  real<lower = 0.> alpha_scale;
  vector[K] beta_loc;
  vector<lower = 0.>[K] beta_scale;
  real<lower = 0.> sigma_scale;
}
parameters {
  real alpha;
  vector[K] beta;
  real<lower = 0.> sigma;
}
transformed parameters {
  vector[N_obs] mu_obs;
  vector[N_cens] mu_cens;
  mu_obs = alpha + X_obs * beta;
  mu_cens = alpha + X_cens * beta;
}
model {
  sigma ~ cauchy(0, sigma_scale);
  alpha ~ normal(alpha_loc, alpha_scale);
  beta ~ normal(beta_loc, beta_scale);
  y_obs ~ normal(mu_obs, sigma);
  target += normal_lccdf(y_cens | mu_cens, sigma);
}

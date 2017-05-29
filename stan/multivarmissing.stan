data {
  int N;
  // X
  int<lower = 0, upper = N> N_x_obs;
  int<lower = 0, upper = N> N_x_miss;
  int<lower = 1, upper = N> x_obs_idx[N_x_obs];
  vector[N_x_obs] x_obs;
  int<lower = 1, upper = N> x_miss_idx[N_x_miss];
  // Y
  int<lower = 0, upper = N> N_y_obs;
  int<lower = 0, upper = N> N_y_miss;
  int<lower = 1, upper = N> y_obs_idx[N_y_obs];
  vector[N_y_obs] y_obs;
  int<lower = 1, upper = N> y_miss_idx[N_y_miss];
  // Z
  int<lower = 0, upper = N> N_z_obs;
  int<lower = 0, upper = N> N_z_miss;
  int<lower = 1, upper = N> z_obs_idx[N_z_obs];
  vector[N_z_obs] z_obs;
  int<lower = 1, upper = N> z_miss_idx[N_z_miss];
  // priors
  real mu_z;
  real<lower = 0.> sigma_z;
  real<lower = 0.> sigma_x_scale;
  real<lower = 0.> sigma_y_scale;
  vector[2] gamma_loc;
  vector<lower = 0.>[2] gamma_scale;
  vector[3] beta_loc;
  vector<lower = 0.>[3] beta_scale;
}
parameters {
  vector[2] gamma;
  vector[3] beta;
  vector[N_x_miss] x_miss;
  vector[N_y_miss] y_miss;
  vector[N_z_miss] z_miss;
  real<lower = 0.> sigma_x;
  real<lower = 0.> sigma_y;
}
transformed parameters {
  vector[N] x;
  vector[N] y;
  vector[N] z;
  vector[N] mu_x;
  vector[N] mu_y;
  x[x_miss_idx] = x_miss;
  x[x_obs_idx] = x_obs;
  y[y_miss_idx] = y_miss;
  y[y_obs_idx] = y_obs;
  z[z_miss_idx] = z_miss;
  z[z_obs_idx] = z_obs;
  mu_x = gamma[1] + gamma[2] * z;
  mu_y = beta[1] + beta[2] * x + beta[3] * z;
}
model {
  x ~ normal(mu_x, sigma_x);
  y ~ normal(mu_y, sigma_y);
  z_miss ~ normal(mu_z, sigma_z);
  gamma ~ normal(gamma_loc, gamma_scale);
  beta ~ normal(beta_loc, beta_scale);
	sigma_x ~ cauchy(0., sigma_x_scale);
	sigma_y ~ cauchy(0., sigma_y_scale);
}

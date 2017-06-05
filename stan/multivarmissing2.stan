data {
  // number of obs
  int<lower = 1> N;
  // number of variables
  int<lower = 1> K;
  // X
  int<lower = 0, upper = N * K> N_obs;
  vector[N_obs] X_obs;
  int<lower = 1, upper = N> X_obs_row[N_obs];
  int<lower = 1, upper = K> X_obs_col[N_obs];
  int<lower = 0, upper = N * K> N_miss;
  int<lower = 1, upper = N> X_miss_row[N_miss];
  int<lower = 1, upper = N> X_miss_col[N_miss];
  // priors
  vector<lower = 0.>[K] Sigma_scale_scale;
  real<lower = 0.> Sigma_corr_L_eta;
  vector[K] mu_loc;
  vector<lower = 0.>[K] mu_scale;
}
parameters {
  vector[K] mu;
  vector<lower = 0.>[K] Sigma_scale;
  cholesky_factor_corr[K] Sigma_corr_L;
  vector[N_miss] X_miss;
}
transformed parameters {
  vector[K] X[N];
  for (i in 1:N_obs) {
   X[X_obs_row[i], X_obs_col[i]] = X_obs[i];
  }
  for (i in 1:N_miss) {
    X[X_miss_row[i], X_miss_col[i]] = X_miss[i];
  }
}
model {
  Sigma_corr_L ~ lkj_corr_cholesky(Sigma_corr_L_eta);
  Sigma_scale ~ cauchy(0., Sigma_scale_scale);
  for (i in 1:N) {
    X[i] ~ multi_normal_cholesky(mu, Sigma_corr_L);
  }
}

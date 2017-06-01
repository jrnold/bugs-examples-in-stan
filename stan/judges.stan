data {
  // number of items
  int<lower = 1> K;
  // number of individuals
  int<lower = 1> N;
  // observed votes
  int<lower = 0, upper = N * K> Y_obs;
  int y_idx_leg[Y_obs];
  int y_idx_vote[Y_obs];
  int y[Y_obs];
  // ideal points
  vector[N] xi_loc;
  vector<lower = 0.>[N] xi_scale;
  vector[N] xi_skew;
  // priors
  vector[K] alpha_loc;
  vector<lower = 0.>[K] alpha_scale;
  vector[K] beta_loc;
  vector<lower = 0.>[K] beta_scale;
}
parameters {
  // item difficulties
  vector[K] alpha;
  // item cutpoints
  vector[K] beta;
  // unknown ideal points
  vector[N] xi;
}
transformed parameters {
  // create xi from observed and parameter ideal points
  vector[Y_obs] mu;
  for (i in 1:Y_obs) {
    mu[i] = beta[y_idx_vote[i]] * xi[y_idx_leg[i]] - alpha[y_idx_vote[i]];
  }
}
model {
  alpha ~ normal(alpha_loc, alpha_scale);
  beta ~ normal(beta_loc, beta_scale);
  xi ~ skew_normal(xi_loc, xi_scale, xi_skew);
  y ~ binomial_logit(1, mu);
}
generated quantities {
  vector[Y_obs] log_lik;
  for (i in 1:Y_obs) {
    log_lik[i] = binomial_logit_lpmf(y[i] | 1, mu[i]);
  }
}

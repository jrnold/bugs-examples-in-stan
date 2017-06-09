// ideal point model
//
// identification:
// - ideal points ~ normal(0, 1)
// - signs of item discrimination using skew normal
data {
  // number of individuals
  int N;
  // number of items
  int K;
  // observed votes
  int<lower = 0, upper = N * K> Y_obs;
  int y_idx_leg[Y_obs];
  int y_idx_vote[Y_obs];
  int y[Y_obs];
  // priors
  // on items
  real alpha_loc;
  real<lower = 0.> alpha_scale;
  vector[K] beta_loc;
  vector<lower = 0.>[K] beta_scale;
  vector[K] beta_skew;
}
parameters {
  // item difficulties
  vector[K] alpha;
  // item discrimination
  vector[K] beta;
  // unknown ideal points
  vector[N] xi_raw;
}
transformed parameters {
  // create xi from observed and parameter ideal points
  vector[Y_obs] mu;
  vector[N] xi;
  xi = (xi_raw - mean(xi_raw)) ./ sd(xi_raw);
  for (i in 1:Y_obs) {
    mu[i] = alpha[y_idx_vote[i]] + beta[y_idx_vote[i]] * xi[y_idx_leg[i]];
  }
}
model {
  alpha ~ normal(alpha_loc, alpha_scale);
  beta ~ skew_normal(beta_loc, beta_scale, beta_skew);
  // soft center ideal points
  // in transformed block enforce hard-centering
  xi_raw ~ normal(0., 1.);
  y ~ bernoulli_logit(mu);
}
generated quantities {
  vector[Y_obs] log_lik;
  for (i in 1:Y_obs) {
    log_lik[i] = bernoulli_logit_lpmf(y[i] | mu[i]);
  }
}

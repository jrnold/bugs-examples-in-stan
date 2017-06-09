// ideal point model
// identification:
// - ideal points ~ normal(0, 1)
// - signs of item discrimination using bounds
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
  int K_beta_pos;
  int<lower = 1, upper = K> beta_idx_pos[K_beta_pos];
  int K_beta_neg;
  int<lower = 1, upper = K> beta_idx_neg[K_beta_neg];
  int K_beta_unc;
  int<lower = 1, upper = K> beta_idx_unc[K_beta_unc];
}
parameters {
  // item difficulties
  vector[K] alpha;
  // item discrimination
  vector<lower = 0.>[K_beta_pos] beta_pos;
  vector<upper = 0.>[K_beta_neg] beta_neg;
  vector[K_beta_unc] beta_unc;
  // unknown ideal points
  vector[N] xi_raw;
}
transformed parameters {
  // create xi from observed and parameter ideal points
  vector[Y_obs] mu;
  vector[N] xi;
  vector[K] beta;

  beta[beta_idx_neg] = beta_neg;
  beta[beta_idx_pos] = beta_pos;
  beta[beta_idx_unc] = beta_unc;
  xi = (xi_raw - mean(xi)) / sd(xi);
  for (i in 1:Y_obs) {
    mu[i] = alpha[y_idx_vote[i]] + beta[y_idx_vote[i]] * xi[y_idx_leg[i]];
  }

}
model {
  alpha ~ normal(alpha_loc, alpha_scale);
  beta_neg ~ normal(beta_loc, beta_scale);
  beta_pos ~ normal(beta_loc, beta_scale);
  beta_unc ~ normal(beta_loc, beta_scale);
  xi_raw ~ normal(0., 1.);
  y ~ bernoulli_logit(mu);
}
generated quantities {
  vector[Y_obs] log_lik;

  for (i in 1:Y_obs) {
    log_lik[i] = bernoulli_logit_lpmf(y[i] | mu[i]);
  }
}

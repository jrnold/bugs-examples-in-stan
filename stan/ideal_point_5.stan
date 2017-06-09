// ideal point model
// identification:
// - xi ~ normal(0, 1)
// - signs of xi
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
  int N_xi_pos;
  int<lower = 1, upper = N> xi_idx_pos[N_xi_pos];
  int N_xi_neg;
  int<lower = 1, upper = N> xi_idx_neg[N_xi_neg];
  int N_xi_unc;
  int<lower = 1, upper = N> xi_idx_unc[N_xi_unc];
}
parameters {
  // item difficulties
  vector[K] alpha;
  // item discrimination
  vector[K] beta;
  // unknown ideal points
  vector<lower = 0.>[N_xi_pos] xi_pos;
  vector<upper = 0.>[N_xi_neg] xi_neg;
  vector[N_xi_unc] xi_unc;
}
transformed parameters {
  // create xi from observed and parameter ideal points
  vector[Y_obs] mu;
  vector[N] xi;
  xi[xi_idx_neg] = xi_neg;
  xi[xi_idx_pos] = xi_pos;
  xi[xi_idx_unc] = xi_unc;
  for (i in 1:Y_obs) {
    mu[i] = alpha[y_idx_vote[i]] + beta[y_idx_vote[i]] * xi[y_idx_leg[i]];
  }
}
model {
  alpha ~ normal(alpha_loc, alpha_scale);
  beta ~ normal(beta_loc, beta_scale);
  xi_neg ~ normal(0., 1.);
  xi_pos ~ normal(0., 1.);
  xi_unc ~ normal(0., 1.);
  y ~ bernoulli_logit(mu);
}
generated quantities {
  vector[Y_obs] log_lik;
  for (i in 1:Y_obs) {
    log_lik[i] = bernoulli_logit_lpmf(y[i] | mu[i]);
  }
}

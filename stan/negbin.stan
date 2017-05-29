data {
  int N;
  int y[N];
  int K;
  matrix[N, K] X;
  // priors
  real alpha_mean;
  real<lower = 0.> alpha_scale;
  vector[K] beta_mean;
  vector<lower = 0.>[K] beta_scale;
  real<lower = 0.> reciprocal_phi_scale;
}
parameters {
  real alpha;
  vector[K] beta;
  real<lower = 0.> reciprocal_phi;
}
transformed parameters {
  vector[N] eta;
  real<lower = 0.> phi;
  eta = alpha + X * beta;
  phi = 1. / reciprocal_phi;
}
model {
  reciprocal_phi ~ cauchy(0., reciprocal_phi_scale);
  alpha ~ normal(alpha_mean, alpha_scale);
  beta ~ normal(beta_mean, beta_scale);
  y ~ neg_binomial_2_log(eta, phi);
}
generated quantities {
  vector[N] mu;
  vector[N] log_lik;
  vector[N] y_rep;
  mu = exp(eta);
  for (i in 1:N) {
    log_lik[i] = neg_binomial_2_log_lpmf(y[i] | eta[i], phi);
    y_rep[i] = neg_binomial_2_rng(mu[i], phi);
  }
}

data {
  // response
  int N;
  int y[N];
  // covariates
  int K;
  matrix[N, K] X;
  // priors
  real alpha_mean;
  real<lower = 0.> alpha_scale;
  vector[K] beta_mean;
  vector<lower = 0.>[K] beta_scale;
}
parameters {
  real alpha;
  vector[K] beta;
}
transformed parameters {
  // linear predictor
  vector[N] eta;
  vector[N] mu;
  eta = alpha + X * beta;
  // mu = Phi(eta);
  // Phi_approx is faster
  mu = Phi_approx(eta);
}
model {
  alpha ~ normal(alpha_mean, alpha_scale);
  beta ~ normal(beta_mean, beta_scale);
  y ~ bernoulli(mu);
}
generated quantities {
  // log-likelihood of each obs
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = bernoulli_lpmf(y[i] | mu[i]);
  }
}

functions {
  real pct_correct_pred(int[] y, vector mu) {
    real out;
    int N;
    N = num_elements(mu);
    out = 0.;
    for (i in 1:N) {
      if (y[i]) {
        out = out + int_step(mu[i] >= 0.5);
      } else {
        out = out + int_step(mu[i] < 0.5);
      }
    }
    out = out / N;
    return out;
  }
  real expected_pct_correct_pred(int[] y, vector mu) {
    real out;
    int N;
    N = num_elements(mu);
    out = 0.;
    for (i in 1:N) {
      if (y[i]) {
        out = out + mu[i];
      } else {
        out = out + (1. - mu[i]);
      }
    }
    out = out / N;
    return out;
  }
}
data {
  // response
  int N;
  int y[N];
  // covariates
  int K;
  matrix[N, K] X;
  // priors
  real alpha_loc;
  real<lower = 0.> alpha_scale;
  vector[K] beta_loc;
  vector<lower = 0.>[K] beta_scale;
}
parameters {
  real alpha;
  vector[K] beta;
}
transformed parameters {
  // linear predictor
  vector[N] eta;
  eta = alpha + X * beta;
}
model {
  alpha ~ normal(alpha_loc, alpha_scale);
  beta ~ normal(beta_loc, beta_scale);
  // y ~ bernoulli(inv_logit(eta));
  // this is faster and more numerically stable
  y ~ bernoulli_logit(eta);
}
generated quantities {
  // log-likelihood of each obs
  vector[N] log_lik;
  // probability
  vector[N] mu;
  // simulated outcomes
  int y_rep[N];
  // percent correctly predicted
  real PCP;
  real ePCP;
  for (i in 1:N) {
    mu[i] = inv_logit(eta[i]);
    log_lik[i] = bernoulli_logit_lpmf(y[i] | eta[i]);
    y_rep[i] = bernoulli_rng(mu[i]);
  }
  PCP = pct_correct_pred(y, mu);
  ePCP = expected_pct_correct_pred(y, mu);
}

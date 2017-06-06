data {
  // multivariate outcome
  int<lower = 1> N;
  int<lower = 2> K;
  vector[K] y[N];
  // covariates
  int<lower = 0> P;
  vector[P] X[N];
  // prior
  vector[K] alpha_loc;
  vector<lower = 0.>[K] alpha_scale;
  vector[P] beta_loc[K];
  vector<lower = 0.>[P] beta_scale[K];
  real<lower = 0.> Sigma_corr_shape;
  real<lower = 0.> Sigma_scale_scale;
}
parameters {
  // regression intercept
  vector[K] alpha;
  // regression coefficients
  vector[P] beta[K];
  // Cholesky factor of the correlation matrix
  cholesky_factor_corr[K] Sigma_corr_L;
  vector<lower = 0.>[K] Sigma_scale;
  // student-T degrees of freedom
  real<lower = 2.> nu;
}
transformed parameters {
  vector[K] mu[N];
  matrix[K, K] Sigma;
  // covariance matrix
  Sigma = crossprod(diag_pre_multiply(Sigma_scale, Sigma_corr_L));
  for (i in 1:N) {
    for (k in 1:K) {
      mu[i, k] = alpha[k] + dot_product(X[i], beta[k]);
    }
  }
}
model {
  for (k in 1:K) {
    alpha[k] ~ normal(alpha_loc[k], alpha_scale[k]);
    beta[k] ~ normal(beta_loc[k], beta_scale[k]);
  }
  nu ~ gamma(2, 0.1);
  Sigma_scale ~ cauchy(0., Sigma_scale_scale);
  Sigma_corr_L ~ lkj_corr_cholesky(Sigma_corr_shape);
  y ~ multi_student_t(nu, mu, Sigma);
}

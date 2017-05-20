data {
  // multivariate outcome
  int<lower = 1> N;
  int<lower = 2> K;
  vector[K] y[N];
  // covariates
  int<lower = 0> P;
  vector[P] X[N];
  // prior
  vector[K] a_location;
  vector<lower = 0.>[K] a_scale;
  vector[P] b_location[K];
  vector<lower = 0.>[P] b_scale[K];
  real<lower = 0.> Sigma_shape;
}
parameters {
  vector[K] a;
  vector[P] b[K];
  corr_matrix[K] Sigma;
  real<lower = 2.> nu;
}
transformed parameters {
  vector[K] mu[N];
  for (i in 1:N) {
    for (k in 1:K) {
      mu[i, k] = a[k] + dot_product(X[i], b[k]);
    }
  }
}
model {
  for (k in 1:K) {
    a[k] ~ normal(a_location[k], a_scale[k]);
    b[k] ~ normal(b_location[k], b_scale[k]);
  }
  nu ~ gamma(2, 0.1);
  Sigma ~ lkj_corr(Sigma_shape);
  for (i in 1:N) {
    y[i] ~ multi_student_t(nu, mu[i], Sigma);
  }
}

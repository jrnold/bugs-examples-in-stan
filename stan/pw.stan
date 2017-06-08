data {
  // number of observations
  // need at least two to estimates
  int<lower = 2> N;
  // response
  vector[N] y;
  // regression design matrix
  int<lower = 1> K;
  matrix[N, K] X;
  // priors
  real alpha_loc;
  real<lower = 0.> alpha_scale;
  vector[K] beta_loc;
  vector<lower = 0.>[K] beta_scale;
  real<lower = 0.> sigma_scale;
  real<lower = 0.> theta_a;
  real<lower = 0.> theta_b;
}
parameters {
  // regression coefficients
  real alpha;
  vector[K] beta;
  // error scale
  real<lower=0> sigma;
  // lag coefficients
  real<lower = 0, upper = 1> theta_raw;
}
transformed parameters {
  // observation means
  vector[N] mu;
  // lag coefficient;
  real<lower = -1, upper = 1> theta;
  // convert range of theta from (0, 1) to (-1, 1)
  theta = (2. * theta_raw - 1.);
  // regression
  mu[1] = alpha + dot_product(beta, X[1, ]);
  mu[2:N] = alpha * (1 - theta) + (X[2:N, ] - theta * X[1:(N - 1), ]) * beta;
}
model {
  alpha ~ cauchy(alpha_loc, alpha_scale);
  beta ~ cauchy(beta_loc, beta_scale);
  theta_raw ~ beta(theta_a, theta_b);
  sigma ~ cauchy(0, sigma_scale);
  y[1] ~ normal(mu[1], sigma / sqrt(1 + theta ^ 2));
  y[2:N] ~ normal(mu[2:N], sigma);
}

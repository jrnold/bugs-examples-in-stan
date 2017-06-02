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
  // error terms
  vector[N] epsilon;
  // lag coefficient;
  real<lower = -1, upper = 1> theta;
  // convert range of theta from (0, 1) to (-1, 1)
  theta = (2. * theta_raw - 1.);
  // regression
  mu = alpha + X * beta;
  // construct errors
  epsilon[1] = y[1] - mu[1];
  for (i in 2:N) {
    epsilon[i] = y[i] - mu[i] - theta * epsilon[i - 1];
  }
}
model {
  alpha ~ cauchy(alpha_loc, alpha_scale);
  beta ~ cauchy(beta_loc, beta_scale);
  theta_raw ~ beta(theta_a, theta_b);
  sigma ~ cauchy(0, sigma_scale);
  for (i in 2:N) {
    y[i] ~ normal(mu[i] + theta * epsilon[i - 1], sigma);
  }
}

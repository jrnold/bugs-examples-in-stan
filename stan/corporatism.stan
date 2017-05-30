data {
  // number of observations
  int N;
  // response variable
  vector[N] y;
  // number of predictors in the regression
  int K;
  // design matrix of country-year obs
  matrix[N, K] X;
  // number of countries
  int n_country;
  // countries for each observation
  int<lower = 1, upper = n_country> country[N];
  // design matrix of country-variables
  int J;
  matrix[n_country, J] U;
  // priors
  // mean and scale of normal prior on beta
  vector[K] beta_mean;
  vector<lower = 0.>[K] beta_scale;
  // mean and scale of normal prior on gamma
  real gamma_mean;
  real<lower = 0.> gamma_scale;
  // scale for half-Cauchy prior on tau
  real<lower = 0.> tau_scale;
}
parameters {
  // obs. errors.
  real<lower = 0.> sigma;
  // country-specific terms
  vector[n_country] gamma;
  vector[J] delta;
  // regression coefficients
  vector[K] beta[n_country];
  // scale on country priors
  real<lower = 0.> tau;
}
transformed parameters {
  vector[N] mu;
  vector[n_country] alpha;
  alpha = gamma + U * delta;
  for (i in 1:N) {
    mu[i] = alpha[country[i]] + X[i] * beta[country[i]];
  }
}
model {
  gamma ~ normal(gamma_mean, gamma_scale);
  tau ~ cauchy(0., tau_scale);
  for (k in 1:K) {
    beta[k] ~ normal(beta_mean, beta_scale);
  }
  alpha ~ normal(gamma, tau);
  y ~ normal(mu, sigma);
}
generated quantities {
}

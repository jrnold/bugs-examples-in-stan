data {
  int N;
  vector[N] y;
  int K;
  matrix[N, K] X;
  int Y;
  int<lower = 1, upper = Y> year[N];
  // priors
  real sigma_scale;
  vector[K] beta_loc;
  vector<lower = 0.>[K] beta_scale;
  real alpha_loc;
  real<lower = 0.> alpha_scale;
}
parameters {
  vector[Y] alpha;
  vector[K] beta;
  real<lower = 2.> nu;
  real<lower = 0.> sigma;
  real<lower = 0.> tau;
}
transformed parameters {
  vector[N] mu;
  for (i in 1:N) {
    mu[i] = alpha[year[i]] + X[i] * beta;
  }
}
model{
  // priors for error variance
  sigma ~ cauchy(0., sigma_scale);
  // priors for year intercepts
  alpha ~ normal(alpha_loc, alpha_scale);
	// priors for the regression coefficients
	beta ~ normal(beta_loc, beta_scale);
	// degrees of freedom
	nu ~ gamma(2, 0.1);
	// likelihood
	y ~ student_t(nu, mu, sigma);
}
generated quantities {
  real delta;
  delta = beta[3] + beta[4];
}

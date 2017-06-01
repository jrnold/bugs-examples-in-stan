data {
  int N;
  vector[N] y;
  real U;
  real mu_mean;
  real mu_scale;
  real sigma_scale;
}
parameters {
  real mu;
  real<lower = 0.> sigma;
}
model {
  mu ~ normal(mu_mean, mu_scale);
  sigma ~ cauchy(0., sigma_scale);
  for (i in 1:N) {
    y[i] ~ normal(mu, sigma) T[, U];
  }
}

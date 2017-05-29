data {
  real<lower = 0., upper = 100.> y;
  real<lower = 0.> y_sd;
  real<lower = 0., upper = 100.> mu_mean;
  real<lower = 0.> mu_sd;
}
parameters {
  real mu;
}
model {
  mu ~ normal(mu_mean, mu_sd);
  y ~ normal(mu, y_sd);
}

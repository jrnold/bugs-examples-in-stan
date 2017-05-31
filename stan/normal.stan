data {
  int N;
  vector[N] y;
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
  y ~ normal(mu, sigma);
}

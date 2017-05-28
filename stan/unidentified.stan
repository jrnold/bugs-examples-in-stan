data {
  real y;
  vector[2] theta_mean;
  vector<lower = 0.>[2] theta_scale;
}
parameters {
  vector[2] theta;
}
transformed parameters {
  real mu;
  mu = sum(theta);
}
model {
  y ~ normal(mu, 1.);
  theta ~ normal(theta_mean, theta_scale);
}

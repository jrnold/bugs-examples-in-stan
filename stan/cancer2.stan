data {
  int<lower = 0> r[2];
  int<lower = 1> n[2];
  // param for beta prior on p
  real a_loc;
  real<lower = 0.> a_scale;
  real b_loc;
  real<lower = 0.> b_scale;
}
parameters {
  real a;
  real b;
}
transformed parameters {
  vector<lower = 0., upper = 1.>[2] p;
  p[1] = inv_logit(a + b);
  p[2] = inv_logit(a);
}
model {
  a ~ normal(a_loc, a_scale);
  b ~ normal(a_loc, b_scale);
  r ~ binomial(n, p);
}
generated quantities {
  real delta;
  int delta_up;
  real lambda;
  int lambda_up;

  delta = p[1] - p[2];
  delta_up = delta > 0;
  lambda = logit(p[1]) - logit(p[2]);
  lambda_up = lambda > 0;

}

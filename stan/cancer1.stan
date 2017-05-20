data {
  int<lower = 0> r[2];
  int<lower = 1> n[2];
  // param for beta prior on p
  vector<lower = 0.>[2] p_a;
  vector<lower = 0.>[2] p_b;
}
parameters {
  vector<lower = 0., upper = 1.>[2] p;
}
model {
  p ~ beta(p_a, p_b);
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

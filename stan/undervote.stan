data {
  int n[4];
  int y[4];
  vector[4] pi_a;
  vector[4] pi_b;
}
parameters {
  vector<lower = 0., upper = 1.>[4] pi;
}
model {
  y ~ binomial(n, pi);
  pi ~ beta(pi_a, pi_b);
}
generated quantities {
  vector[2] delta;
  int good[2];
  delta[1] = pi[2] - pi[1];
  delta[2] = pi[4] - pi[3];
  good[1] = int_step(delta[1]);
  good[2] = int_step(delta[2]);
}

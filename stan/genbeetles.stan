data {
  int N;
  int r[N];
  int n[N];
  vector[N] x;
}
parameters {
  real alpha;
  real beta;
  real<lower = 0.> nu;
}
transformed parameters {
  vector[N] mu;
  for (i in 1:N) {
    mu[i] = pow(inv_logit(alpha + beta * x[i]), nu) ;
  }
}
model {
  alpha ~ normal(0., 10.);
  beta ~ normal(0., 2.5);
  nu ~ gamma(0.25, 0.25);
  r ~ binomial(n, mu);
}
generated quantities {
  // probability where the maximum marginal effect
  real pdot;
  pdot = pow(inv_logit(nu), nu);
}

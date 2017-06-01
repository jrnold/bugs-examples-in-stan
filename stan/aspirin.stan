data {
  int N;
  vector[N] y;
  vector[N] s;
  real mu_loc;
  real<lower = 0.> mu_scale;
  real<lower = 0.> tau_scale;
  real<lower = 0.> tau_df;
}
parameters {
  vector[N] theta;
  real mu;
  real<lower = 0.> tau;
}
model {
  mu ~ normal(mu_loc, mu_scale);
  tau ~ student_t(tau_df, 0., tau_scale);
  theta ~ normal(mu, tau);
  y ~ normal(theta, s);
}
generated quantities {
  vector[N] shrinkage;
  {
    real tau2;
    tau2 = pow(tau, 2.);
    for (i in 1:N) {
      real v;
      v = pow(s[i], 2);
      shrinkage[i] = v / (v + tau2);
    }
  }
}

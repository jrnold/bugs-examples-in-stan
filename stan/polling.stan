data {
  int<lower = 1> N;
  int<lower = 1> T;
  vector[N] y;
  vector<lower = 0.>[N] s;
  int<lower = 1, upper = T> time[N];
  int H;
  int<lower = 1, upper = H> house[N];
  // initial and final values
  real xi_init;
  real xi_final;
  vector[H] delta_mean;
  vector<lower = 0.>[H] delta_scale;
}
parameters {
  vector[T] omega;
  real<lower = 0.> tau;
  vector[H] delta;
}
transformed parameters {
  vector[T - 1] xi;
  vector[N] mu;
  xi[1] = xi_init;
  for (i in 2:(T - 1)) {
    xi[i] = xi[i - 1] + tau * omega[i];
  }
  for (i in 1:N) {
    mu[i] = xi[time[i]] + delta[house[i]];
  }
}
model {
  delta ~ normal(delta_mean, delta_scale);
  xi_final ~ normal(xi[T - 1], tau);
  y ~ normal(mu, s);
}

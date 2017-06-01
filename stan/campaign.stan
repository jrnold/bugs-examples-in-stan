// polling model
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
  real delta_loc;
  real<lower = 0.> zeta_scale;
  real<lower = 0.> tau_scale;
}
parameters {
  vector[T - 1] omega;
  real<lower = 0.> tau;
  vector[H] delta_raw;
  real<lower = 0.> zeta;
}
transformed parameters {
  vector[H] delta;
  vector[T - 1] xi;
  vector[N] mu;
  // this is necessary. If not centered the model is unidentified
  delta = (delta_raw - mean(delta_raw)) / sd(delta_raw) * zeta;
  xi[1] = xi_init;
  for (i in 2:(T - 1)) {
    xi[i] = xi[i - 1] + tau * omega[i - 1];
  }
  for (i in 1:N) {
    mu[i] = xi[time[i]] + delta[house[i]];
  }
}
model {
  // house effects
  delta_raw ~ normal(0., 1.);
  zeta ~ normal(0., zeta_scale);
  // latent state innovations
  omega ~ normal(0., 1.);
  // scale of innovations
  tau ~ cauchy(0, tau_scale);
  // final known effect
  xi_final ~ normal(xi[T - 1], tau);
  // daily polls
  y ~ normal(mu, s);
}

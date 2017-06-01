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
  vector[2] xi_init;
  vector[2] xi_final;
  real delta_loc;
  real<lower = 0.> zeta_scale;
  vector<lower = 0.>[2] tau_scale;
}
parameters {
  vector[2] omega[T - 1];
  vector<lower = 0.>[2] tau;
  vector[H] delta_raw;
  real<lower = 0.> zeta;
}
transformed parameters {
  vector[H] delta;
  vector[2] xi[T - 1];
  vector[N] mu;
  // this is necessary. If not centered the model is unidentified
  delta = (delta_raw - mean(delta_raw)) / sd(delta_raw) * zeta;
  xi[1] = xi_init;
  for (i in 2:(T - 1)) {
    // slope needs to be defined before the original data
    xi[i, 2] = xi[i - 1, 2] + tau[2] * omega[i - 1, 2];
    xi[i, 1] = xi[i - 1, 1] + xi[i, 2] + tau[1] * omega[i - 1, 1];
  }
  for (i in 1:N) {
    mu[i] = xi[time[i], 1] + delta[house[i]];
  }
}
model {
  // house effects
  delta_raw ~ normal(0., 1.);
  zeta ~ normal(0., zeta_scale);
  // latent state innovations
  for (i in 1:size(omega)) {
    omega[i] ~ normal(0., 1.);
  }
  // scale of innovations
  tau ~ normal(0, tau_scale);
  // final known effect
  xi_final ~ normal(xi[T - 1], tau);
  // daily polls
  y ~ normal(mu, s);
}

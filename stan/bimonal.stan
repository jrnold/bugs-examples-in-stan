data {
  int N;
  int K;
  int<lower = 0, upper=N> N_obs;
  int<lower = 0, upper=N> N_miss;
  vector[N_obs] x_obs;
  int<lower = 1, upper = K> x_obs_idx[2, N_obs];
  int<lower = 1, upper = K> x_miss_idx[2, N_obs];
}
parameters {
  vector[K] mu;
  cov_matrix[K] Sigma;
  vector[N_miss] x_miss;
}
transformed parameters {
  vector[K] X[N];
  for (i in 1:N_obs) {
    X[x_obs_idx[i, 1], x_obs_idx[i, 2]] = x_obs[i];
  }
  for (i in 1:N_miss) {
    X[x_miss_idx[i, 1], x_miss_idx[i, 2]] = x_miss[i];
  }
}
model{
  for (i in 1:N) {
    X[i] ~ multi_normal(mu, Sigma);
  }
}

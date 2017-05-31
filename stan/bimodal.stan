data {
  int N;
  int<lower = 0, upper=N> N_obs;
  int<lower = 0, upper=N> N_miss;
  vector[N_obs] x_obs;
  int<lower = 1, upper = N> x_obs_idx[N_obs, 2];
  int<lower = 1, upper = N> x_miss_idx[N_miss, 2];
  vector[2] mu;
}
parameters {
  cov_matrix[2] Sigma;
  vector[N_miss] x_miss;
}
transformed parameters {
  // using an array of vectors is more convenient when sampling
  // multi_normal than using an matrix
  vector[2] X[N];
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

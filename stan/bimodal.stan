data {
  // number of obs
  int N;
  int<lower = 0, upper=N * 2> N_obs;
  int<lower = 0, upper=N * 2> N_miss;
  vector[N_obs] x_obs;
  int<lower = 1, upper = N> x_obs_row[N_obs];
  int<lower = 1, upper = 2> x_obs_col[N_obs];
  int<lower = 1, upper = N> x_miss_row[N_miss];
  int<lower = 1, upper = 2> x_miss_col[N_miss];
  real<lower = 1.> df;
}
parameters {
  vector[2] mu;
  cov_matrix[2] Sigma;
  vector[N_miss] x_miss;
}
transformed parameters {
  // using an array of vectors is more convenient when sampling
  // multi_normal than using an matrix
  vector[2] X[N];
  for (i in 1:N_obs) {
    X[x_obs_row[i], x_obs_col[i]] = x_obs[i];
  }
  for (i in 1:N_miss) {
    X[x_miss_row[i], x_miss_col[i]] = x_miss[i];
  }
}
model{
  for (i in 1:N) {
    X[i] ~ multi_student_t(df, mu, Sigma);
  }
}

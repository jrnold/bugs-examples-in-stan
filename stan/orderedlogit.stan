data {
  // number of observations
  int N;
  // number of response categories
  int K;
  // response
  int<lower = 1, upper = K> y[N];
  // regression design matrix
  int D;
  matrix[N, D] X;
}
parameters {
  // ordered logistic distribution cutpoints
  vector[K - 1] gamma;
  // intercept and coefficients in regression
  real alpha;
  vector[P] beta;
}
transformed parameters {
  // linear predictor in logit scale;
  vector[N] eta;
  eta = alpha + X * beta;
}
model {
  for (i in 1:N) {
    y[i] ~ ordered_logistic(eta[i], gamma);
  }
}

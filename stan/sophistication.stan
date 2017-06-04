data {
  // number of respondents
  int N;
  // number of items
  int K;
  // binary responses
  int<lower = 0, upper = 1> y_bern[K, N];
  // interviewer overall rating
  vector[N] y_norm;
  // interviewers
  int J;
  int<lower = 1, upper = J> interviewer[N];
  // priors
  real alpha_loc;
  real<lower = 0.> alpha_scale;
  real beta_loc;
  real<lower = 0.> beta_scale;
  real<lower = 0.> gamma_scale;
  real<lower = 0.> sigma_scale;
  real<lower = 0.> tau_scale;
  real delta_loc;
  real delta_scale;
}
parameters {
  // respondent latent score
  vector[N] xi_raw;
  // item discrimination
  vector[K] beta;
  // item difficulty
  vector[K] alpha;
  // coefficient in interviewer rating
  real<lower = 0.> gamma;
  // error in interviewer rating
  real<lower = 0.> sigma;
  // interviewer random effects
  vector[J] nu;
  // location of interviewer random effects
  real delta;
  // scale of interviewer random effects
  real<lower = 0.> tau;
}
transformed parameters {
  // interviewer rating
  vector[N] theta;
  // abilities
  vector[N] xi;
  xi = (xi_raw - mean(xi_raw));
  // respondent latent score
  for (i in 1:N) {
    theta[i] = gamma * xi[i] + nu[interviewer[i]];
  }
}
model {
  // priors
  xi_raw ~ normal(0., 1.);
  beta ~ normal(beta_loc, beta_scale);
  alpha ~ normal(alpha_loc, alpha_scale);
  gamma ~ normal(0., gamma_scale);
  sigma ~ cauchy(0., sigma_scale);
  tau ~ cauchy(0., tau_scale);
  delta ~ normal(delta_loc, delta_scale);
  // binary responses
  for (k in 1:K) {
    y_bern[k] ~ bernoulli_logit(beta[k] * xi - alpha[k]);
  }
  // interviewer random effects
  nu ~ normal(delta, tau);
  // interviewer score
  y_norm ~ normal(theta, sigma);
}

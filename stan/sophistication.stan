data {
  // number of respondents
  int N;
  // number of items
  int K;
  // binary responses
  int<lower = 0, upper = 1> y[N, K];
  // interviewer overall rating
  vector[N] z;
  // interviewers
  int J;
  int<lower = 1, upper = J> interviewer[N];
  // priors
  vector[N] xi_loc;
  vector<lower = 0.>[N] xi_scale;
  vector[K] beta_loc;
  vector<lower = 0.>[K] beta_scale;
  vector[K] alpha_loc;
  vector<lower = 0.>[K] alpha_scale;
  real<lower = 0.> gamma_scale;
  real<lower = 0.> sigma_scale;
  real<lower = 0.> tau_scale;
  real delta_loc;
  real delta_scale;
}
parameters {
  // respondent latent score
  vector[N] xi;
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
  // respondent latent score
  for (i in 1:N) {
    theta[i] = gamma * xi[i] + nu[interviewer[i]];
  }
}
model {
  // priors
  xi ~ normal(xi_loc, xi_scale);
  beta ~ normal(beta_loc, beta_scale);
  alpha ~ normal(alpha_loc, alpha_scale);
  gamma ~ normal(0., gamma_scale);
  sigma ~ cauchy(0., sigma_scale);
  tau ~ cauchy(0., tau_scale);
  delta ~ normal(delta_loc, delta_scale);
  // binary responses
  for (k in 1:K) {
    y[k] ~ bernoulli_logit(beta[k] * xi - alpha[k]);
  }
  // interviewer random effects
  nu ~ normal(delta, tau);
  // interviewer score
  z ~ normal(theta, sigma);
}

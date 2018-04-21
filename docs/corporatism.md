
# Corporatism: Hierarchical model for economic growth {#corporatism}


```r
library("rstan")
library("tidyverse")
```

The following program implements a regression model of economic growth among 16 OECD countries, 1971-1984 [@Western1998a, @AlvarezGarrettLange1991a].[^corporatism-src]
The model is hierarchical in that it specifies country-specific coefficients for the following predictors: lagged growth, demand, import price movements, export price movements, leftist government and an intercept.
The magnitudes of the country-specific coefficients are conditional on (time-invariant) extent of labor organization within each country; these regression relationships constitute the second level of the model.

The data come from N=16 countries, and $T=14$ years (1971:1984) with $K=6$ covariates at the lowest ("micro") level of the hierarchy, and $J=2$ covariates (an intercept and the labor organization variable) at the second level.


```r
data("corporatism", package = "bayesjackman")
```


```r
corporatism_country <- corporatism %>%
  dplyr::select(country, labor.org) %>%
  distinct()
```


```r
corporatism_mod <- stan_model("stan/corporatism.stan")
```

<pre>
  <code class="stan">data {
  // number of observations
  int N;
  // response variable
  vector[N] y;
  // number of predictors in the regression
  int K;
  // design matrix of country-year obs
  matrix[N, K] X;
  // number of countries
  int n_country;
  // countries for each observation
  int<lower = 1, upper = n_country> country[N];
  // design matrix of country-variables
  int J;
  matrix[n_country, J] U;
  // priors
  // mean and scale of normal prior on beta
  vector[K] beta_mean;
  vector<lower = 0.>[K] beta_scale;
  // mean and scale of normal prior on gamma
  real gamma_mean;
  real<lower = 0.> gamma_scale;
  // scale for half-Cauchy prior on tau
  real<lower = 0.> tau_scale;
}
parameters {
  // obs. errors.
  real<lower = 0.> sigma;
  // country-specific terms
  vector[n_country] gamma;
  vector[J] delta;
  // regression coefficients
  vector[K] beta[n_country];
  // scale on country priors
  real<lower = 0.> tau;
}
transformed parameters {
  vector[N] mu;
  vector[n_country] alpha;
  alpha = gamma + U * delta;
  for (i in 1:N) {
    mu[i] = alpha[country[i]] + X[i] * beta[country[i]];
  }
}
model {
  gamma ~ normal(gamma_mean, gamma_scale);
  tau ~ cauchy(0., tau_scale);
  for (k in 1:K) {
    beta[k] ~ normal(beta_mean, beta_scale);
  }
  alpha ~ normal(gamma, tau);
  y ~ normal(mu, sigma);
}
generated quantities {
}</code>
</pre>

[^corporatism-src]: Example derived from Simon Jackman, "[Corporatism: hierarchical or 'multi-level' model for economic growth in 16 OECD countries](https://web-beta.archive.org/web/20070724034043/http://jackman.stanford.edu/mcmc/corporatism.odc)", 2007-07-24.


# Resistant: Outlier-resistant regression via the Student's $t$ distribution (#resistant)


```r
library("tidyverse")
library("rstan")
```


Outlying data points can distort estimates of location, such as means or regression coefficients.  Location estimates obtained via maximizing a iid normal likelihood over heavy tailed data will be sensitive to data in the tails (outliers). A popular alternative to normal errors in regression analyses is the Student's $t$ density, with an unknown degrees of freedom parameter.  For low degrees of freedom, the t has heavier tails than the normal, but tends to the normal as the degrees of freedom parameter increases.  Treating the degrees of freedom parameter as an unknown parameter to be estimaetd thus provides a check on the appropriateness of the normal.  By embedding a model with location parameters in the Student's $t$ density, we obtain outlier-resistant estimates of location parameters.

To illustrate these ideas, I use data collected by Douglas Grob on incumbency advantage in American congressional elections, 1956-1994.   The dependent variable is the proportion of the two-party vote won by the Democratic candidate in district $i$ at election $t$ (uncontested districts are dropped from the analysis).  Indicators for Democratic and Republican incumbency are the critical explanatory variables in the analysis; coefficients on these indicators are regarded as estimates of incumbency advantage.  A series of year-specific indicators ("fixed effects") are also included in the specification.

$$
\begin{aligned}[t]
y_i &\sim \mathsf{StudentT}(\nu, \mu_i, \sigma) \\
\mu_i &= \alpha + x_i \beta
\end{aligned}
$$
The $\alpha$, $\beta$, and $\sigma$ parameters are given weakly informative priors (assuming that all $x$ are scaled to have mean 0, standard deviation 1):
$$
\begin{aligned}[t]
\alpha &\sim \mathsf{Normal}(\bar{y}, 10 s_y), \\
\beta &\sim \mathsf{Normal}(0, 2.5 s_y), \\
\sigma &\sim \mathsf{HalfCauchy}(0, 5 s_y) ,
\end{aligned}
$$
where $\bar{y}$ is the mean of $y$, and $s_y$ is the standard deviation of $y$.
The degrees of freedom parameter in the Student-t distribution is a paramter and given the weakly informative prior suggested by @JuarezSteel2010a,
$$
\nu \sim \mathsf{Gamma}(2, 0.1) .
$$
The following code operationalizes this regression model; the conditional density of the vote proportions is $t$, with unknown degrees of freedom, $\nu.$ 



```r
resistant_mod <- stan_model("stan/resistant.stan")
```
<pre>
  <code class="stan">data {
  int N;
  vector[N] y;
  int K;
  matrix[N, K] X;
  int Y;
  int<lower = 1, upper = Y> year[N];
  // priors
  real sigma_scale;
  vector[K] beta_loc;
  vector<lower = 0.>[K] beta_scale;
  real alpha_loc;
  real<lower = 0.> alpha_scale;
}
parameters {
  vector[Y] alpha;
  vector[K] beta;
  real<lower = 2.> nu;
  real<lower = 0.> sigma;
  real<lower = 0.> tau;
}
transformed parameters {
  vector[N] mu;
  for (i in 1:N) {
    mu[i] = alpha[year[i]] + X[i] * beta;
  }
}
model{
  // priors for error variance
  sigma ~ cauchy(0., sigma_scale);
  // priors for year intercepts
  alpha ~ normal(alpha_loc, alpha_scale);
	// priors for the regression coefficients
	beta ~ normal(beta_loc, beta_scale);
	// degrees of freedom
	nu ~ gamma(2, 0.1);
	// likelihood
	y ~ student_t(nu, mu, sigma);
}
generated quantities {
  real delta;
  delta = beta[3] + beta[4];
}</code>
</pre>



```r
load("data/resistant.rda")
resistant_data <- within(list(), {
  y <- resistant$y
  N <- length(y)
  X <- model.matrix(~ 0 + lagy + prvwinpty + deminc + repinc, data = resistant) %>%
    scale()
  K <- ncol(X)
  year <- resistant$year
  Y <- max(year)
  # priors
  alpha_loc <- mean(y)
  alpha_scale <- 10 * sd(y)
  beta_loc <- rep(0, K)
  beta_scale <- rep(2.5 * sd(y), K)
  sigma_scale <- 5 * sd(y)
})
```


```r
resistant_fit <- sampling(resistant_mod, data = resistant_data)
#> Warning: There were 3994 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See
#> http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> Warning: Examine the pairs() plot to diagnose sampling problems
```

```r
summary(resistant_fit, par = c("nu", "sigma", "beta", "tau"))$summary
#>               mean se_mean    sd       2.5%        25%        50%
#> nu        6.54e+00  0.0483 0.643   5.47e+00   6.05e+00   6.48e+00
#> sigma     5.82e+00  0.0074 0.102   5.63e+00   5.75e+00   5.82e+00
#> beta[1]   1.16e+01  0.0278 0.190   1.13e+01   1.15e+01   1.16e+01
#> beta[2]  -3.19e+00  0.1396 0.364  -3.92e+00  -3.44e+00  -3.15e+00
#> beta[3]   3.92e+00  0.0347 0.266   3.46e+00   3.72e+00   3.92e+00
#> beta[4]  -3.90e+00  0.0312 0.196  -4.31e+00  -4.03e+00  -3.89e+00
#> tau      8.75e+307     Inf   Inf  4.88e+306  4.23e+307  8.56e+307
#>                75%      97.5%   n_eff Rhat
#> nu        6.97e+00   7.95e+00  177.47 1.03
#> sigma     5.90e+00   6.02e+00  190.65 1.03
#> beta[1]   1.17e+01   1.20e+01   46.77 1.08
#> beta[2]  -2.94e+00  -2.57e+00    6.78 1.22
#> beta[3]   4.11e+00   4.45e+00   58.65 1.08
#> beta[4]  -3.78e+00  -3.52e+00   39.69 1.08
#> tau      1.33e+308  1.76e+308 4000.00  NaN
```

## Questions {-}

1. How does using the Student-t distribution compare to using a normal distribution for the errors?


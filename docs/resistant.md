
# Resistant: Outlier-resistant regression via the Student's $t$ distribution {#resistant}


```r
library("tidyverse")
library("rstan")
```

Outlying data points can distort estimates of location, such as means or regression coefficients.[^resistant-src]
Location estimates obtained via maximizing a iid normal likelihood over heavy tailed data will be sensitive to data in the tails (outliers).
A popular alternative to normal errors in regression analyses is the Student's $t$ density, with an unknown degrees of freedom parameter.
For low degrees of freedom, the Student's $t$ distribution has heavier tails than the normal, but tends to the normal as the degrees of freedom parameter increases.
Treating the degrees of freedom parameter as an unknown parameter to be estimated thus provides a check on the appropriateness of the normal.
By embedding a model with location parameters in the Student's $t$ density, we obtain outlier-resistant estimates of location parameters.

## Data

This example uses data collected by Douglas Grob on incumbency advantage in American congressional elections, 1956-1994 [@Jackman2000a].

```r
data("resistant", package = "bayesjackman")
glimpse(resistant)
#> List of 6
#>  $ y        : num [1:5090] 46.3 54.3 58.5 57.8 58.5 ...
#>  $ lagy     : num [1:5090] 57 46.3 54.3 58.5 70 ...
#>  $ prvwinpty: num [1:5090] 1 -1 1 1 1 1 1 1 1 1 ...
#>  $ deminc   : num [1:5090] 0 0 1 1 1 1 0 1 1 1 ...
#>  $ repinc   : num [1:5090] 0 1 0 0 0 0 0 0 0 0 ...
#>  $ year     : num [1:5090] 1 2 3 4 6 7 8 10 11 12 ...
```

The response variable is the proportion of the two-party vote won by the Democratic candidate in district $i$ at election $t$.
Indicators for Democratic and Republican incumbency are the critical explanatory variables in the analysis.
Coefficients on these indicators are regarded as estimates of incumbency advantage.
A series of year-specific indicators (*fixed effects*) are also included in the specification.

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
The degrees of freedom parameter in the Student's $t$ distribution is a parameter and given the weakly informative prior suggested by @JuarezSteel2010a,
$$
\nu \sim \mathsf{Gamma}(2, 0.1) .
$$
The following code operationalizes this regression model.
The conditional density of the vote proportions is $t$, with unknown degrees of freedom, $\nu.$


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
#> Warning: There were 3988 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See
#> http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> Warning: Examine the pairs() plot to diagnose sampling problems
```

```r
summary(resistant_fit, par = c("nu", "sigma", "beta", "tau"))$summary
#>               mean se_mean     sd       2.5%        25%        50%
#> nu        6.44e+00 0.03924 0.5612   5.44e+00   6.05e+00   6.40e+00
#> sigma     5.82e+00 0.00596 0.0989   5.62e+00   5.75e+00   5.82e+00
#> beta[1]   1.17e+01 0.02438 0.1776   1.14e+01   1.16e+01   1.17e+01
#> beta[2]  -3.35e+00 0.03354 0.3064  -3.96e+00  -3.55e+00  -3.35e+00
#> beta[3]   3.97e+00 0.02955 0.2527   3.44e+00   3.81e+00   3.97e+00
#> beta[4]  -3.96e+00 0.02970 0.2104  -4.39e+00  -4.09e+00  -3.96e+00
#> tau      8.79e+307     Inf    Inf  5.17e+306  4.28e+307  8.69e+307
#>                75%      97.5%  n_eff Rhat
#> nu        6.78e+00   7.63e+00  204.5 1.01
#> sigma     5.88e+00   6.01e+00  275.4 1.01
#> beta[1]   1.18e+01   1.20e+01   53.1 1.10
#> beta[2]  -3.17e+00  -2.75e+00   83.4 1.02
#> beta[3]   4.13e+00   4.49e+00   73.2 1.05
#> beta[4]  -3.83e+00  -3.55e+00   50.2 1.06
#> tau      1.31e+308  1.74e+308 4000.0  NaN
```

## Reparaterization: standard deviation instead of scale

In the Student's $t$ distribution, the standard deviation is a function of the degrees of freedom. 
For degrees of freedom $\nu > 2$, the variance is defined, and
$$
\sigma^* = sd(y) = \sigma \sqrt{ \frac{\nu}{\nu - 2}}
$$
This makes the sampling of $\nu$ and $\sigma$ a priori dependent.
Instead, we can place priors on the degrees of freedom $\nu$ and the standard deviation $\sigma^*$,
and treat $\sigma$ as a transformed parameter,
$$
\begin{aligned}
\sigma^* &\sim \mathsf{HalfCauchy}{(0, 5)} \\
\sigma &= \sigma^* \sqrt{\frac{\nu - 2}{\nu}} \\
\end{aligned}
$$


```r
resistant_mod2 <- stan_model("stan/resistant2.stan")
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
  real<lower = 0.> sigma_raw;
  real<lower = 0.> tau;
}
transformed parameters {
  vector[N] mu;
  real sigma;
  for (i in 1:N) {
    mu[i] = alpha[year[i]] + X[i] * beta;
  }
  // paramterization so sigma and
  sigma = sigma_raw * sqrt((nu - 2) / nu);
}
model{
  // priors for the standard deviation
  sigma_raw ~ cauchy(0., sigma_scale);
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
resistant_fit2 <- sampling(resistant_mod2, data = resistant_data)
#> Warning: There were 3987 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See
#> http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> Warning: Examine the pairs() plot to diagnose sampling problems
```

```r
summary(resistant_fit2, par = c("beta", "sigma", "nu", "tau"))$summary
#>               mean se_mean    sd       2.5%        25%        50%
#> beta[1]   1.17e+01 0.01250 0.180   1.13e+01   1.15e+01   1.17e+01
#> beta[2]  -3.32e+00 0.03435 0.313  -3.98e+00  -3.52e+00  -3.31e+00
#> beta[3]   3.94e+00 0.02430 0.236   3.45e+00   3.79e+00   3.94e+00
#> beta[4]  -3.97e+00 0.02449 0.214  -4.37e+00  -4.12e+00  -3.97e+00
#> sigma     5.82e+00 0.00512 0.101   5.62e+00   5.75e+00   5.81e+00
#> nu        6.48e+00 0.04144 0.602   5.44e+00   6.06e+00   6.43e+00
#> tau      8.93e+307     Inf   Inf  3.73e+306  4.31e+307  9.00e+307
#>                75%      97.5%  n_eff Rhat
#> beta[1]   1.18e+01   1.20e+01  208.7 1.01
#> beta[2]  -3.10e+00  -2.75e+00   83.1 1.01
#> beta[3]   4.09e+00   4.41e+00   94.6 1.03
#> beta[4]  -3.83e+00  -3.54e+00   76.7 1.05
#> sigma     5.89e+00   6.01e+00  392.1 1.01
#> nu        6.84e+00   7.77e+00  211.1 1.02
#> tau      1.36e+308  1.75e+308 4000.0  NaN
```




## Questions {-}

<<<<<<< HEAD:docs/resistant.md
1.  How does using the Student-t distribution compare to using a normal distribution for the errors?
=======
1. How does using the Student-t distribution compare to using a normal distribution for the errors? How would you evaluate it.
2. Compare the effective sample size, $\hat{R}$, and speed for $\sigma$ and $\nu$ when using the scale/degrees of freedom and standard deviation/degrees of freedom parameterizations.
>>>>>>> 2c0f6aa8a3ee03e62597a3eac350fc216e98d7a3:_book/resistant.md

[^resistant-src]: Example derived from Simon Jackman, "Resistant: Outlier-resistant regression via the t distribution," 2007-07-24, [URL](https://web-beta.archive.org/web/20070724034107/http://jackman.stanford.edu:80/mcmc/resistant.odc).

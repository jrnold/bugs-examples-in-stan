
# Truncation: How does Stan deal with truncation?

See @Stan2016a, Chapter 11 "Truncated or Censored Data" for more on how Stan handles truncation and censoring.

```r
library("tidyverse")
library("rstan")
```


Assume we have the observations, $y = 1,...,9$, from a Normal population with unknown mean and variance, subject to the constraint that $y < 10$,
$$
\begin{aligned}[t]
y &\sim \mathsf{Normal}(\mu, \sigma^2) I(-\infty, 10) .
\end{aligned}
$$

Ignoring the constraint, the MLEs for the mean and variance are 5 and 6.67; with the constraint taken into account, each observation makes likelihood contribution
$$
f (y; m, s_2)/F ((k - m)/s),
$$
where $k$ is the truncation point (in this case, 10), and the MLEs of $m, s_2$  are 5.32 and 8.28.

The posterior of this model is not well identified by the data, so the mean, $\mu$, and scale, $\sigma$, are given informative priors based on the data,
$$
\begin{aligned}[t]
\mu &\sim \mathsf{Normal}(\bar{y}, s_y) ,\\
\sigma &\sim \mathsf{HalfCauchy}(0, s_y) .
\end{aligned}
$$
where $\bar{y}$ is the mean of $y$, and $s_y$ is the standard deviation of $y$.


```r
truncation_mod <- stan_model("stan/SingleTruncation.stan")
```
<pre>
  <code class="stan">data {
  int N;
  vector[N] y;
  real U;
  real mu_mean;
  real mu_scale;
  real sigma_scale;
}
parameters {
  real mu;
  real<lower = 0.> sigma;
}
model {
  mu ~ normal(mu_mean, mu_scale);
  sigma ~ cauchy(0., sigma_scale);
  for (i in 1:N) {
    y[i] ~ normal(mu, sigma) T[, U];
  }
}</code>
</pre>


```r
truncation_data <- within(list(), {
  y <- 1:9
  N <- length(y)
  U <- 10
  mu_mean <- mean(y)
  mu_scale <- sd(y)
  sigma_scale <- sd(y)
})
```




```r
truncation_fit <- sampling(truncation_mod, data = truncation_data)
```

```r
truncation_fit
#> Inference for Stan model: SingleTruncation.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
#> 
#>         mean se_mean   sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
#> mu      5.82    0.04 1.48   3.39   4.84   5.63   6.63   9.41  1201    1
#> sigma   3.76    0.04 1.39   1.97   2.80   3.46   4.40   7.20  1250    1
#> lp__  -13.54    0.03 1.08 -16.30 -13.99 -13.22 -12.75 -12.44  1258    1
#> 
#> Samples were drawn using NUTS(diag_e) at Wed May 31 08:42:27 2017.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

We can compare these results to that of a model in which the truncation is not taken into account:
$$
\begin{aligned}[t]
y_i &\sim \mathsf{Normal}(\mu, \sigma^2), \\
\mu &\sim \mathsf{Normal}(\bar{y}, s_y) ,\\
\sigma &\sim \mathsf{HalfCauchy}(0, s_y) .
\end{aligned}
$$


```r
truncation_mod2 <- stan_model("stan/normal.stan")
```
<pre>
  <code class="stan">data {
  int N;
  vector[N] y;
  real mu_mean;
  real mu_scale;
  real sigma_scale;
}
parameters {
  real mu;
  real<lower = 0.> sigma;
}
model {
  mu ~ normal(mu_mean, mu_scale);
  sigma ~ cauchy(0., sigma_scale);
  y ~ normal(mu, sigma);
}</code>
</pre>


```r
truncation_fit2 <-
  sampling(truncation_mod2, data = truncation_data)
```

```r
truncation_fit2
#> Inference for Stan model: normal.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
#> 
#>         mean se_mean   sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
#> mu      5.00    0.02 0.93   3.17   4.43   5.01   5.59   6.91  2193    1
#> sigma   2.97    0.02 0.79   1.87   2.42   2.82   3.33   4.96  1836    1
#> lp__  -13.77    0.03 1.05 -16.69 -14.16 -13.45 -13.03 -12.75  1265    1
#> 
#> Samples were drawn using NUTS(diag_e) at Wed May 31 08:42:30 2017.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```



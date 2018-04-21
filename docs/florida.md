
# Florida: Learning About an Unknown Proportion from Survey Data {#florida}


```r
library("tidyverse")
library("rstan")
```

In this example, beliefs about an unknown proportion are updated from new survey data.
The particular example is using survey update beliefs about support for Bush in Florida in the 2000 presidential election campaign [@Jackman2004a].[^florida-src]


```r
florida_mod <- stan_model("stan/florida.stan")
```
<pre>
  <code class="stan">data {
  real<lower = 0., upper = 100.> y;
  real<lower = 0.> y_sd;
  real<lower = 0., upper = 100.> mu_mean;
  real<lower = 0.> mu_sd;
}
parameters {
  real mu;
}
model {
  mu ~ normal(mu_mean, mu_sd);
  y ~ normal(mu, y_sd);
}</code>
</pre>

The prior polls had a mean of 49.1% in support for Bush, with a standard deviation of 2.2%.
The new poll shows 55% support for Bush, with a standard deviation of 2.2%.

```r
florida_data <- list(
  mu_mean = 49.1,
  mu_sd = 2.2,
  y_sd = 2.2,
  y = 55  
)
```


```r
florida_fit <- sampling(florida_mod, data = florida_data)
```

```r
florida_fit
#> Inference for Stan model: florida.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
#> 
#>      mean se_mean   sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
#> mu   52.0    0.04 1.57 48.88 50.98 52.06 53.08  55.1  1449    1
#> lp__ -2.3    0.02 0.72 -4.41 -2.46 -2.02 -1.85  -1.8  2047    1
#> 
#> Samples were drawn using NUTS(diag_e) at Fri Apr 20 00:55:02 2018.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

After observing the new poll, the mean for the posterior is 52, with a 95% credible interval of 48.9--55.1.

[^florida-src]: This example is derived from  Simon Jackman, "Florida," *BUGS Examples,* 2007-07-24, [URL](https://web-beta.archive.org/web/20070724034219/http://jackman.stanford.edu/mcmc/florida.zip).

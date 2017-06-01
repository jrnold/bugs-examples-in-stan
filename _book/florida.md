
# Florida: Learning About an Unknown Proportion from Survey Data {#florida}


```r
library("rstan")
```

This is an learning about an unknown proportion from survey data; using survey data to update beliefs about support for Bush in Florida in the 2000 presidential election campaign [@Jackman2004a].


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
#>       mean se_mean   sd  2.5%   25% 50%   75% 97.5% n_eff Rhat
#> mu   52.04    0.04 1.52 48.99 51.05  52 53.05  55.0  1347    1
#> lp__ -2.27    0.02 0.67 -4.28 -2.41  -2 -1.85  -1.8  1604    1
#> 
#> Samples were drawn using NUTS(diag_e) at Wed May 31 05:18:36 2017.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

After observing the new poll, the mean for the posterior is 52, with a 95% credible interval of 49--55.

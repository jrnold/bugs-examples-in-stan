
# Cancer: difference in two binomial proportions {#cancer}


```r
library("tidyverse")
library("rstan")
```


Two groups chosen to be random samples from subpopulations of lung-cancer patients and cancer-free individuals.[^cancer]
The scientific question of interest is the difference in the smoking habits between two groups.
The results of the survey are:

```r
cancer <- tribble(
  ~group, ~n, ~smokers,
  "Cancer patients", 86, 82,
  "Control group", 86, 72
)
```
# A tibble: 2 x 3
            group     n smokers
            <chr> <dbl>   <dbl>
1 Cancer patients    86      82
2   Control group    86      72

[^cancer]: This example is derived from Simon Jackman, "Cancer: difference in two binomial proportions", *BUGS Examples,* 2007-07-24, http://jackman.stanford.edu:80/mcmc/cancer.odc.  [Wayback Machine](https://web-beta.archive.org/web/20070601000000*/http://jackman.stanford.edu:80/mcmc/cancer.odc). This examples comes from @JohnsonAlbert1999a, using data from @Dorn1954a.


## Two Sample Binomial Model

In implementing this model, we have just two data points (cancer patients and control group) and a binomial sampling model, in which the population proportions of smokers in each group appear as parameters.  Quantities of interest such as the difference in the population proportions and the log of the odds ratio are computed in the generated quantities section. Uniform priors on the population proportions are used in this example.

$$
\begin{aligned}[t]
r_i &\sim \mathsf{Binomial}(n_i, \pi_i)
\end{aligned}
$$
Additionally the difference,
$$
\delta = \pi_1 - \pi_2 ,
$$
and the log-odds ratio,
$$
\lambda = \log\left(\frac{\pi_1}{1 - \pi_1}\right) - \log \left( \frac{\pi_2}{1 - \pi_2} \right) ,
$$

It places uniform priors (Beta priors) are placed on $\pi$,
$$
\begin{aligned}
\pi_i &\sim \mathsf{Beta}(1, 1)
\end{aligned}
$$

The difference between and log odds ratio are defined in the `generated quantities` block.


```r
cancer_data <- list(
  r <- cancer$smokers,
  n <- cancer$n,
  # beta prior on pi
  p_a = rep(1, 2),
  p_b = rep(1, 2)
)
```

The Stan model for this is:

```r
cancer_mod1 <- stan_model("stan/cancer1.stan")
```
<pre>
  <code class="stan">data {
  int<lower = 0> r[2];
  int<lower = 1> n[2];
  // param for beta prior on p
  vector<lower = 0.>[2] p_a;
  vector<lower = 0.>[2] p_b;
}
parameters {
  vector<lower = 0., upper = 1.>[2] p;
}
model {
  p ~ beta(p_a, p_b);
  r ~ binomial(n, p);
}
generated quantities {
  real delta;
  int delta_up;
  real lambda;
  int lambda_up;

  delta = p[1] - p[2];
  delta_up = delta > 0;
  lambda = logit(p[1]) - logit(p[2]);
  lambda_up = lambda > 0;

}</code>
</pre>

Now estimate the model:

```r
cancer_fit1 <- sampling(cancer_mod1, cancer_data)
```

```r
cancer_fit1
#> Inference for Stan model: cancer1.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
#> 
#>             mean se_mean   sd   2.5%    25%    50%    75%  97.5% n_eff
#> p[1]        0.94    0.00 0.02   0.89   0.93   0.95   0.96   0.98  2857
#> p[2]        0.83    0.00 0.04   0.75   0.80   0.83   0.86   0.90  3086
#> delta       0.11    0.00 0.04   0.02   0.08   0.11   0.14   0.20  3023
#> delta_up    1.00    0.00 0.07   1.00   1.00   1.00   1.00   1.00  3652
#> lambda      1.30    0.01 0.55   0.25   0.92   1.28   1.64   2.46  2565
#> lambda_up   1.00    0.00 0.07   1.00   1.00   1.00   1.00   1.00  3652
#> lp__      -60.37    0.02 1.01 -63.09 -60.74 -60.06 -59.65 -59.40  1706
#>           Rhat
#> p[1]         1
#> p[2]         1
#> delta        1
#> delta_up     1
#> lambda       1
#> lambda_up    1
#> lp__         1
#> 
#> Samples were drawn using NUTS(diag_e) at Tue Jun  6 22:50:13 2017.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

## Binomial Logit Model of the Difference

An alternative parameterization directly models the difference in the population proportion.

$$
\begin{aligned}[t]
r_i &\sim \mathsf{Binomial}(n_i, \pi_i) \\
\pi_1 &= \frac{1}{1 + \exp(-(\alpha + \beta)} \\
\pi_2 &= \frac{1}{1 + \exp(-\alpha))}
\end{aligned}
$$
The parameters $\alpha$ and $\beta$ are given weakly informative priors on the log-odds scale,
$$
\begin{aligned}
\alpha &\sim N(0, 10)\\
\beta &\sim N(0, 2.5)
\end{aligned}
$$


```r
cancer_mod2 <- stan_model("stan/cancer2.stan")
```
<pre>
  <code class="stan">data {
  int<lower = 0> r[2];
  int<lower = 1> n[2];
  // param for beta prior on p
  real a_loc;
  real<lower = 0.> a_scale;
  real b_loc;
  real<lower = 0.> b_scale;
}
parameters {
  real a;
  real b;
}
transformed parameters {
  vector<lower = 0., upper = 1.>[2] p;
  p[1] = inv_logit(a + b);
  p[2] = inv_logit(a);
}
model {
  a ~ normal(a_loc, a_scale);
  b ~ normal(a_loc, b_scale);
  r ~ binomial(n, p);
}
generated quantities {
  real delta;
  int delta_up;
  real lambda;
  int lambda_up;

  delta = p[1] - p[2];
  delta_up = delta > 0;
  lambda = logit(p[1]) - logit(p[2]);
  lambda_up = lambda > 0;

}</code>
</pre>

Re-use `r` and `n` values from `cancer_data`, but add the appropriate values for the prior distributions.

```r
cancer_data2 <- within(cancer_data, {
  p_a <- p_b <- NULL
  a_loc <- b_loc <- 0
  a_scale <- 10
  b_scale <- 2.5
})
```

Sample from the model:

```r
cancer_fit2 <- sampling(cancer_mod2, cancer_data2)
```

```r
cancer_fit2
#> Inference for Stan model: cancer2.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
#> 
#>             mean se_mean   sd   2.5%    25%    50%    75%  97.5% n_eff
#> a           1.69    0.01 0.29   1.16   1.50   1.69   1.87   2.29  1926
#> b           1.38    0.01 0.57   0.33   0.99   1.34   1.75   2.57  1725
#> p[1]        0.95    0.00 0.02   0.90   0.94   0.95   0.97   0.98  2796
#> p[2]        0.84    0.00 0.04   0.76   0.82   0.84   0.87   0.91  1995
#> delta       0.11    0.00 0.04   0.03   0.08   0.11   0.14   0.20  1594
#> delta_up    1.00    0.00 0.07   1.00   1.00   1.00   1.00   1.00  1913
#> lambda      1.38    0.01 0.57   0.33   0.99   1.34   1.75   2.57  1725
#> lambda_up   1.00    0.00 0.07   1.00   1.00   1.00   1.00   1.00  1913
#> lp__      -55.51    0.02 1.00 -58.16 -55.83 -55.19 -54.80 -54.57  1638
#>           Rhat
#> a            1
#> b            1
#> p[1]         1
#> p[2]         1
#> delta        1
#> delta_up     1
#> lambda       1
#> lambda_up    1
#> lp__         1
#> 
#> Samples were drawn using NUTS(diag_e) at Tue Jun  6 22:50:17 2017.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```


## Questions

1. Expression the Binomial Logit model of the Difference as a regression
2. What number of success and failures is a `Beta(1,1)` prior equivalent to?

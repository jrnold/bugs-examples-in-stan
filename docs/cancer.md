
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
  group               n smokers
  <chr>           <dbl>   <dbl>
1 Cancer patients   86.     82.
2 Control group     86.     72.

## Two Sample Binomial Model

In implementing this model, we have just two data points (cancer patients and
control group) and a binomial sampling model, in which the population
proportions of smokers in each group appear as parameters.  Quantities of
interest such as the difference in the population proportions and the log of
the odds ratio are computed in the generated quantities section. Uniform priors
on the population proportions are used in this example.

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
#> p[1]        0.94    0.00 0.03   0.88   0.93   0.95   0.96   0.98  3589
#> p[2]        0.83    0.00 0.04   0.74   0.80   0.83   0.86   0.90  2325
#> delta       0.11    0.00 0.05   0.02   0.08   0.11   0.14   0.21  2565
#> delta_up    0.99    0.00 0.07   1.00   1.00   1.00   1.00   1.00  3168
#> lambda      1.31    0.01 0.58   0.24   0.92   1.27   1.66   2.51  2706
#> lambda_up   0.99    0.00 0.07   1.00   1.00   1.00   1.00   1.00  3168
#> lp__      -60.44    0.03 1.12 -63.52 -60.81 -60.08 -59.66 -59.40  1521
#>           Rhat
#> p[1]         1
#> p[2]         1
#> delta        1
#> delta_up     1
#> lambda       1
#> lambda_up    1
#> lp__         1
#> 
#> Samples were drawn using NUTS(diag_e) at Fri Apr 20 00:53:46 2018.
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
#> a           1.69    0.01 0.30   1.13   1.49   1.68   1.88   2.31  1639
#> b           1.38    0.01 0.58   0.30   0.99   1.35   1.74   2.61  1683
#> p[1]        0.95    0.00 0.02   0.90   0.94   0.95   0.97   0.98  2817
#> p[2]        0.84    0.00 0.04   0.76   0.82   0.84   0.87   0.91  1720
#> delta       0.11    0.00 0.04   0.03   0.08   0.11   0.14   0.20  1525
#> delta_up    0.99    0.00 0.07   1.00   1.00   1.00   1.00   1.00  2140
#> lambda      1.38    0.01 0.58   0.30   0.99   1.35   1.74   2.61  1683
#> lambda_up   0.99    0.00 0.07   1.00   1.00   1.00   1.00   1.00  2140
#> lp__      -55.55    0.02 1.01 -58.32 -55.94 -55.24 -54.82 -54.57  1625
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
#> Samples were drawn using NUTS(diag_e) at Fri Apr 20 00:54:23 2018.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

## Questions

1.  Expression the Binomial Logit model of the Difference as a regression
1.  What number of success and failures is a `Beta(1,1)` prior equivalent to?

[^cancer]: This example is derived from Simon Jackman,
    "[Cancer: difference in two binomial proportions](https://web-beta.archive.org/web/20070601000000*/http://jackman.stanford.edu:80/mcmc/cancer.odc)",
    *BUGS Examples,* 2007-07-24, This examples comes from @JohnsonAlbert1999a, using data from @Dorn1954a.

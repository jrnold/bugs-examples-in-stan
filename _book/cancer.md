
# Cancer: difference in two binomial proportions {#cancer}


```r
library("tidyverse")
library("rstan")
```


The following simple model is drawn from an example in @JohnsonAlbert1999a, using data collected in a @Dorn1954a.  A sample of 86 lung-cancer patients and a sample of 86 controls were questioned about their smoking habits.  The two groups were chosen to represent random samples from a subpopulation of lung-cancer patients and an otherwise similar population of cancer-free individuals.  Of the cancer patients, 83 out of 86 were smokers; among the control group 72 out of 86 were smokers.  The scientific question of interest was to assess the difference between the smoking habits in the two groups.

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
  r = c(83, 72),
  n = c(86, 86),
  # beta prior on pi
  p_a = rep(1, 2),
  p_b = rep(1, 2)
)
```



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
#> p[1]        0.95    0.00 0.02   0.90   0.94   0.96   0.97   0.99  2639
#> p[2]        0.83    0.00 0.04   0.75   0.81   0.83   0.86   0.90  2863
#> delta       0.12    0.00 0.04   0.04   0.10   0.12   0.15   0.21  2874
#> delta_up    1.00    0.00 0.05   1.00   1.00   1.00   1.00   1.00  2947
#> lambda      1.55    0.01 0.61   0.44   1.14   1.51   1.94   2.84  2290
#> lambda_up   1.00    0.00 0.05   1.00   1.00   1.00   1.00   1.00  2947
#> lp__      -57.46    0.03 1.04 -60.25 -57.83 -57.14 -56.73 -56.48  1564
#>           Rhat
#> p[1]         1
#> p[2]         1
#> delta        1
#> delta_up     1
#> lambda       1
#> lambda_up    1
#> lp__         1
#> 
#> Samples were drawn using NUTS(diag_e) at Wed May 31 05:17:34 2017.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

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


```r
cancer_data2 <- within(cancer_data, {
  p_a <- p_b <- NULL
  a_loc <- b_loc <- 0
  a_scale <- 10
  b_scale <- 2.5
})
```


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
#> a           1.68    0.01 0.29   1.15   1.49   1.67   1.87   2.29  1849
#> b           1.71    0.01 0.65   0.56   1.27   1.67   2.11   3.07  1945
#> p[1]        0.96    0.00 0.02   0.91   0.95   0.97   0.98   0.99  2751
#> p[2]        0.84    0.00 0.04   0.76   0.82   0.84   0.87   0.91  1881
#> delta       0.12    0.00 0.04   0.04   0.09   0.12   0.15   0.21  1764
#> delta_up    1.00    0.00 0.04   1.00   1.00   1.00   1.00   1.00  4000
#> lambda      1.71    0.01 0.65   0.56   1.27   1.67   2.11   3.07  1945
#> lambda_up   1.00    0.00 0.04   1.00   1.00   1.00   1.00   1.00  4000
#> lp__      -52.44    0.03 1.01 -55.22 -52.83 -52.13 -51.72 -51.47  1587
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
#> Samples were drawn using NUTS(diag_e) at Wed May 31 05:18:31 2017.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

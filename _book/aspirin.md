
# Aspirin: Borrowing Strength via Hierarchical Modeling {#aspirin}



```r
library("tidyverse")
library("rstan")
```

The following data come from a meta-analysis of heart attack data [^aspirin-src],

```r
aspirin <- 
  tibble(y = c(2.77, 2.50, 1.84, 2.56, 2.31, -1.15),
         sd = c(1.65, 1.31, 2.34, 1.67, 1.98, 0.90))
```
Each observation is the results of a study of survivorship following a heart attack (myocardial infarction).
In each study, some victims were given aspirin immediately following their heart attack, while some victims were not.
The observed values of $y$ are the differences in mean survivorship observed in each study.
Additionally each study provided a standard deviations, derived from the relative sizes of the two groups in each study.
$$
\begin{aligned}[t]
y_i &\sim \mathsf{Normal}(\theta_i, s_i) , \\
\theta_i &\sim \mathsf{Normal}(\mu, \tau) ,
\end{aligned}
$$
where $y_i$ is the mean of each study, and $s_i$ is the standard deviation for each study.
Weakly informative priors are given to the parameters $\mu$ and $\tau$,
$$
\begin{aligned}[t]
\mu &\sim \mathsf{Normal}(\bar{y}, 10 s_y) , \\
\tau &\sim \mathsf{HalfCauchy}(0, 5 s_y) ,
\end{aligned}
$$
where $\bar{y}$ is the mean of $y$, and $s_y$ is the standard deviation of $y$.

Although the data are binomial, the sample sizes are large enough in each study that the normal approximation is valid.
This simplifies the problem by reducing each study's data to the observed treatment effect and a standard deviation. 
The goal of the meta-analysis is to synthesize the six studies, in order to arrive at an overall estimate of the effects of aspirin on survivorship following a heart attack.

This is a simple example of hierarchical modeling. 
Via the exchangeability assumption, that the study-specific means have a common prior, the studies "borrow strength" from one another.
This introducing some bias, since each study's mean mean is shrunk back towards the common mean.
However, the benefit is gaining precision (smaller variance). 
We also gain a better estimate of the overall effect of aspirin on survivorship after heart attack than we would get from naively pooling the studies or using the estimate of any one study.

[^aspirin-src]: This example is derived from Simon Jackman, "Aspirin: Shrinkage (or "borrowing strength") via hierarchical modeling", 2007-07-24 [URL](https://web-beta.archive.org/web/20070724034135/http://jackman.stanford.edu/mcmc/aspirin.odc). The data and the meta-analysis is from @Draper1992a.

The Stan model for the above model is:

```r
aspirin_mod <- stan_model("stan/aspirin.stan")
```
<pre>
  <code class="stan">data {
  int N;
  vector[N] y;
  vector[N] s;
  real mu_loc;
  real<lower = 0.> mu_scale;
  real<lower = 0.> tau_scale;
  real<lower = 0.> tau_df;
}
parameters {
  vector[N] theta;
  real mu;
  real<lower = 0.> tau;
}
model {
  mu ~ normal(mu_loc, mu_scale);
  tau ~ student_t(tau_df, 0., tau_scale);
  theta ~ normal(mu, tau);
  y ~ normal(theta, s);
}
generated quantities {
  vector[N] shrinkage;
  {
    real tau2;
    tau2 = pow(tau, 2.);
    for (i in 1:N) {
      real v;
      v = pow(s[i], 2);
      shrinkage[i] = v / (v + tau2);
    }
  }
}</code>
</pre>




```r
aspirin_data <- within(list(), {
  y <- aspirin$y
  N <- nrow(aspirin)
  s <- aspirin$sd
  mu_loc <- mean(y)
  mu_scale <- 5 * sd(y)
  tau_scale <- 2.5 * sd(y)
  tau_df <- 4
})
```



```r
aspirin_fit <- sampling(aspirin_mod, data = aspirin_data)
#> Warning: There were 39 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See
#> http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> Warning: Examine the pairs() plot to diagnose sampling problems
```

```r
aspirin_fit
#> Inference for Stan model: aspirin.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
#> 
#>               mean se_mean   sd   2.5%   25%   50%   75% 97.5% n_eff Rhat
#> theta[1]      2.07    0.03 1.27  -0.25  1.21  1.98  2.88  4.72  2115 1.00
#> theta[2]      2.03    0.02 1.11  -0.14  1.28  2.00  2.81  4.25  2038 1.00
#> theta[3]      1.60    0.03 1.51  -1.29  0.62  1.56  2.56  4.73  2948 1.00
#> theta[4]      1.99    0.03 1.27  -0.42  1.12  1.92  2.84  4.57  2350 1.00
#> theta[5]      1.82    0.03 1.34  -0.82  0.94  1.79  2.69  4.51  2746 1.00
#> theta[6]     -0.43    0.03 0.92  -2.25 -1.04 -0.44  0.20  1.36   732 1.00
#> mu            1.51    0.02 1.04  -0.42  0.85  1.48  2.13  3.64  1758 1.00
#> tau           1.80    0.03 0.97   0.47  1.13  1.62  2.27  4.22   807 1.00
#> shrinkage[1]  0.52    0.01 0.22   0.13  0.35  0.51  0.68  0.92   570 1.01
#> shrinkage[2]  0.42    0.01 0.22   0.09  0.25  0.40  0.57  0.88   448 1.01
#> shrinkage[3]  0.65    0.01 0.20   0.24  0.51  0.68  0.81  0.96   714 1.00
#> shrinkage[4]  0.52    0.01 0.22   0.14  0.35  0.52  0.68  0.93   575 1.01
#> shrinkage[5]  0.59    0.01 0.21   0.18  0.43  0.60  0.75  0.95   640 1.00
#> shrinkage[6]  0.29    0.01 0.20   0.04  0.14  0.24  0.39  0.78   332 1.01
#> lp__         -7.52    0.12 2.61 -13.12 -9.09 -7.39 -5.78 -2.64   448 1.00
#> 
#> Samples were drawn using NUTS(diag_e) at Tue Jun  6 23:38:43 2017.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

Note that this model is likely to produce divergent transitions.
The reasons for this and an alternative parameterization is discussed in the next section.


## Non-centered parameterization

For few data, when there are not many groups, or when inter-group variation is high, it can be more efficient to use the non-centered parameterization. See @Stan2016a [p. 331] and @BetancourtGirolami2013a for a more detailed discussion of this.

The non-centered parameterization is
$$
\begin{aligned}[t]
\theta_i^* &\sim \mathsf{Normal}(0, 1) , \\
\theta_i &= \tau \theta^*_i + \mu .
\end{aligned}
$$


```r
aspirin_mod2 <- stan_model("stan/aspirin2.stan")
```


```r
aspirin_mod2
```

<pre>
  <code class="stan">data {
  int N;
  vector[N] y;
  vector[N] s;
  real mu_loc;
  real<lower = 0.> mu_scale;
  real<lower = 0.> tau_scale;
  real<lower = 0.> tau_df;
}
parameters {
  vector[N] theta_raw;
  real mu;
  real<lower = 0.> tau;
}
transformed parameters {
  vector[N] theta;
  theta = tau * theta_raw + mu;
}
model {
  mu ~ normal(mu_loc, mu_scale);
  tau ~ student_t(tau_df, 0., tau_scale);
  theta_raw ~ normal(0., 1.);
  y ~ normal(theta, s);
}
generated quantities {
  vector[N] shrinkage;
  {
    real tau2;
    tau2 = pow(tau, 2.);
    for (i in 1:N) {
      real v;
      v = pow(s[i], 2);
      shrinkage[i] = v / (v + tau2);
    }
  }
}</code>
</pre>


```r
aspirin_fit2 <- sampling(aspirin_mod2, data = aspirin_data,
                        control = list(adapt_delta = 0.99))
```

```r
aspirin_fit2
#> Inference for Stan model: aspirin2.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
#> 
#>               mean se_mean   sd   2.5%   25%   50%   75% 97.5% n_eff Rhat
#> theta_raw[1]  0.30    0.01 0.80  -1.27 -0.23  0.30  0.83  1.85  3096    1
#> theta_raw[2]  0.31    0.02 0.73  -1.11 -0.16  0.31  0.78  1.80  2214    1
#> theta_raw[3]  0.02    0.01 0.87  -1.63 -0.56  0.01  0.60  1.74  3497    1
#> theta_raw[4]  0.27    0.02 0.78  -1.27 -0.24  0.29  0.79  1.86  2479    1
#> theta_raw[5]  0.15    0.02 0.80  -1.43 -0.38  0.14  0.71  1.71  2600    1
#> theta_raw[6] -1.16    0.02 0.76  -2.69 -1.64 -1.16 -0.68  0.27  1670    1
#> mu            1.52    0.03 1.01  -0.41  0.87  1.46  2.10  3.71  1287    1
#> tau           1.69    0.03 0.97   0.20  1.03  1.54  2.17  4.08  1276    1
#> theta[1]      2.02    0.02 1.29  -0.33  1.11  1.93  2.83  4.85  4000    1
#> theta[2]      2.00    0.02 1.08  -0.03  1.24  1.97  2.72  4.18  4000    1
#> theta[3]      1.54    0.03 1.49  -1.36  0.62  1.49  2.47  4.55  3434    1
#> theta[4]      1.96    0.02 1.28  -0.34  1.07  1.86  2.77  4.70  4000    1
#> theta[5]      1.78    0.02 1.34  -0.64  0.84  1.70  2.64  4.68  4000    1
#> theta[6]     -0.34    0.01 0.94  -2.21 -0.97 -0.34  0.32  1.44  4000    1
#> shrinkage[1]  0.54    0.01 0.23   0.14  0.37  0.54  0.72  0.98  1274    1
#> shrinkage[2]  0.46    0.01 0.24   0.09  0.27  0.42  0.62  0.98  1189    1
#> shrinkage[3]  0.68    0.01 0.20   0.25  0.54  0.70  0.84  0.99  1322    1
#> shrinkage[4]  0.55    0.01 0.23   0.14  0.37  0.54  0.73  0.99  1276    1
#> shrinkage[5]  0.62    0.01 0.22   0.19  0.45  0.62  0.79  0.99  1302    1
#> shrinkage[6]  0.32    0.01 0.23   0.05  0.15  0.26  0.43  0.95  1076    1
#> lp__         -5.16    0.08 2.57 -11.25 -6.61 -4.77 -3.28 -1.23   969    1
#> 
#> Samples were drawn using NUTS(diag_e) at Tue Jun  6 23:39:42 2017.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

[^aspirin-example]: This example is derived from Simon Jackman, "Aspirin: Shrinkage (or `borrowing strength') via hierarchical modeling," 2007-07-24, [URL](https://web-beta.archive.org/web/20070724034135/http://jackman.stanford.edu:80/mcmc/aspirin.odc).

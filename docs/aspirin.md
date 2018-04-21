
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
#> Warning: There were 25 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See
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
#> theta[1]      2.10    0.03 1.32  -0.35  1.20  2.03  2.97  4.77  2715    1
#> theta[2]      2.04    0.02 1.12  -0.07  1.27  2.01  2.79  4.26  2065    1
#> theta[3]      1.62    0.03 1.51  -1.27  0.65  1.58  2.55  4.79  2788    1
#> theta[4]      1.98    0.03 1.30  -0.42  1.09  1.90  2.84  4.69  2601    1
#> theta[5]      1.81    0.03 1.41  -0.75  0.84  1.73  2.72  4.82  2919    1
#> theta[6]     -0.44    0.02 0.91  -2.27 -1.04 -0.43  0.16  1.28  2302    1
#> mu            1.52    0.02 1.09  -0.54  0.82  1.50  2.17  3.75  2140    1
#> tau           1.82    0.02 0.95   0.55  1.15  1.62  2.28  4.19  1495    1
#> shrinkage[1]  0.51    0.01 0.21   0.13  0.34  0.51  0.67  0.90  1270    1
#> shrinkage[2]  0.42    0.01 0.21   0.09  0.25  0.40  0.57  0.85  1181    1
#> shrinkage[3]  0.65    0.01 0.20   0.24  0.51  0.68  0.81  0.95  1394    1
#> shrinkage[4]  0.52    0.01 0.21   0.14  0.35  0.52  0.68  0.90  1275    1
#> shrinkage[5]  0.59    0.01 0.21   0.18  0.43  0.60  0.75  0.93  1336    1
#> shrinkage[6]  0.28    0.01 0.18   0.04  0.13  0.24  0.38  0.73  1046    1
#> lp__         -7.66    0.08 2.59 -13.47 -9.23 -7.44 -5.87 -3.26  1033    1
#> 
#> Samples were drawn using NUTS(diag_e) at Fri Apr 20 01:20:48 2018.
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
\theta_i^{*} &\sim \mathsf{Normal}(0, 1) , \\
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
#> theta_raw[1]  0.33    0.01 0.78  -1.21 -0.18  0.32  0.83  1.87  2757    1
#> theta_raw[2]  0.31    0.02 0.76  -1.25 -0.18  0.31  0.79  1.83  2288    1
#> theta_raw[3]  0.06    0.02 0.84  -1.58 -0.50  0.06  0.63  1.72  3065    1
#> theta_raw[4]  0.26    0.02 0.79  -1.26 -0.26  0.25  0.77  1.83  2484    1
#> theta_raw[5]  0.15    0.02 0.81  -1.47 -0.39  0.17  0.70  1.74  2776    1
#> theta_raw[6] -1.14    0.02 0.76  -2.63 -1.64 -1.11 -0.65  0.35  2097    1
#> mu            1.50    0.03 1.04  -0.40  0.84  1.45  2.12  3.71  1381    1
#> tau           1.73    0.03 1.02   0.20  1.06  1.59  2.23  4.33  1141    1
#> theta[1]      2.04    0.02 1.29  -0.27  1.13  1.97  2.89  4.69  4000    1
#> theta[2]      2.00    0.02 1.10   0.02  1.24  1.92  2.73  4.23  4000    1
#> theta[3]      1.59    0.03 1.50  -1.41  0.67  1.54  2.50  4.68  3299    1
#> theta[4]      1.94    0.02 1.29  -0.37  1.02  1.85  2.77  4.69  4000    1
#> theta[5]      1.77    0.02 1.36  -0.76  0.85  1.69  2.61  4.73  4000    1
#> theta[6]     -0.37    0.01 0.95  -2.22 -1.02 -0.37  0.30  1.42  4000    1
#> shrinkage[1]  0.54    0.01 0.23   0.13  0.35  0.52  0.71  0.99  1245    1
#> shrinkage[2]  0.45    0.01 0.24   0.08  0.26  0.40  0.60  0.98  1205    1
#> shrinkage[3]  0.67    0.01 0.21   0.23  0.52  0.68  0.83  0.99  1282    1
#> shrinkage[4]  0.54    0.01 0.23   0.13  0.36  0.52  0.71  0.99  1246    1
#> shrinkage[5]  0.61    0.01 0.22   0.17  0.44  0.61  0.78  0.99  1268    1
#> shrinkage[6]  0.32    0.01 0.24   0.04  0.14  0.24  0.42  0.95  1117    1
#> lp__         -5.13    0.08 2.62 -11.07 -6.64 -4.71 -3.21 -1.18   968    1
#> 
#> Samples were drawn using NUTS(diag_e) at Fri Apr 20 01:21:26 2018.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

## References {-}

This example is derived from Simon Jackman, "Aspirin: Shrinkage (or `borrowing strength') via hierarchical modeling," 2007-07-24, [URL](https://web-beta.archive.org/web/20070724034135/http://jackman.stanford.edu:80/mcmc/aspirin.odc).

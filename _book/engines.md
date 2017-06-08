
# Engines: right-censored failure times


```r
library("tidyverse")
library("rstan")
```


## Data

The data are 40 engines tested at various operating temptures, with the the failure time if the engine failed, or the last time of the observational period if it had not [@Tanner1996a].[^engines-src]
Of the 40 engines, 23 did not fail in their observational periods.

```r
load("data/engines.rda")
glimpse(engines)
#> Observations: 40
#> Variables: 3
#> $ y        <dbl> 3.25, 3.44, 3.54, 3.55, 3.58, 3.69, 3.72, 2.61, 2.61,...
#> $ x        <dbl> 2.26, 2.26, 2.26, 2.26, 2.26, 2.26, 2.26, 2.16, 2.16,...
#> $ censored <lgl> FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALS...
```


## Model

Let $y^*$ be the failure time of engine $i$.
The failure times are modeled as a regression with normal errors,
$$
\begin{aligned}[t]
y^*_i &\sim \mathsf{Normal}(\mu_i, \sigma) , \\
\mu_i &= \alpha + \beta x_i .
\end{aligned}
$$
However, the failure times are not always observed.
In some cases, only the last observation time is known, meaning that all is known is $y^*_i > y_i$.
Let $L$ be the set of censored observation.
$$
\begin{aligned}[t]
y_i &\sim \mathsf{Normal}(\mu_i, \sigma) & i \notin L, \\
y^*_i &\sim \mathsf{Normal}(\mu_i, \sigma) U(y_i, \infty) & i \in L, \\
\mu_i &= \alpha + \beta x_i .
\end{aligned}
$$


$$
\begin{aligned}[t]
\log L(y_i, \dots, y_N | \alpha, \beta, \sigma) &=  \sum_{i \notin L} \log \mathsf{Normal}(y_i; \mu_i, \Sigma) \\
&\quad + \sum_{i \in L} \log \int_{y_i}^{\infty} \mathsf{Normal}(y^*; \mu_i, \Sigma) d\,y^* ,
\end{aligned}
$$
where
$$
\mu_i = \alpha + \beta x .
$$


```r
mod_engines <- stan_model("stan/engines.stan")
```

```r
mod_engines
```

<pre>
  <code class="stan">data {
  // number of observations
  int<lower = 0> N;
  // observed data
  int<lower = 0, upper = N> N_obs;
  vector[N_obs] y_obs;
  // censored data
  int<lower = 0, upper = N> N_cens;
  vector[N_cens] y_cens;
  // covariates
  int<lower = 0> K;
  matrix[N_obs, K] X_obs;
  matrix[N_cens, K] X_cens;
  // priors
  real alpha_loc;
  real<lower = 0.> alpha_scale;
  vector[K] beta_loc;
  vector<lower = 0.>[K] beta_scale;
  real<lower = 0.> sigma_scale;
}
parameters {
  real alpha;
  vector[K] beta;
  real<lower = 0.> sigma;
}
transformed parameters {
  vector[N_obs] mu_obs;
  vector[N_cens] mu_cens;
  mu_obs = alpha + X_obs * beta;
  mu_cens = alpha + X_cens * beta;
}
model {
  sigma ~ cauchy(0, sigma_scale);
  alpha ~ normal(alpha_loc, alpha_scale);
  beta ~ normal(beta_loc, beta_scale);
  y_obs ~ normal(mu_obs, sigma);
  target += normal_lccdf(y_cens | mu_cens, sigma);
}</code>
</pre>


## Estimation

For the input data to the Stan model, the observations that are observed and censored have to be provided in separate vectors.

```r
X <- scale(engines$x)

engines_data <- within(list(), {
  N <- nrow(engines)
  # observed obs
  y_obs <- engines$y[!engines$censored]
  N_obs <- length(y_obs)
  X_obs <- X[!engines$censored, , drop = FALSE]
  K <- ncol(X_obs)
  # censored obs
  y_cens <- engines$y[engines$censored]  
  N_cens <- length(y_cens)
  X_cens <- X[engines$censored, , drop = FALSE]
  # priors
  # use the mean and sd of y to roughly scale the weakly informative
  # priors -- these don't account for  need to exact
  alpha_loc <- mean(engines$y)
  alpha_scale <- 10 * sd(engines$y)
  beta_loc <- array(0)
  beta_scale <- array(2.5 * sd(engines$y))
  sigma_scale <- 5 * sd(y_obs)
})
```


```r
sampling(mod_engines, data = engines_data,
         chains = 1, init = list(list(alpha = mean(engines$y))))
#> 
#> SAMPLING FOR MODEL 'engines' NOW (CHAIN 1).
#> 
#> Gradient evaluation took 4e-05 seconds
#> 1000 transitions using 10 leapfrog steps per transition would take 0.4 seconds.
#> Adjust your expectations accordingly!
#> 
#> 
#> Iteration:    1 / 2000 [  0%]  (Warmup)
#> Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Iteration: 2000 / 2000 [100%]  (Sampling)
#> 
#>  Elapsed Time: 0.067124 seconds (Warm-up)
#>                0.061251 seconds (Sampling)
#>                0.128375 seconds (Total)
#> The following numerical problems occurred the indicated number of times on chain 1
#>                                                                                  count
#> Exception thrown at line 36: normal_lpdf: Scale parameter is 0, but must be > 0!     1
#> When a numerical problem occurs, the Hamiltonian proposal gets rejected.
#> See http://mc-stan.org/misc/warnings.html#exception-hamiltonian-proposal-rejected
#> If the number in the 'count' column is small, there is no need to ask about this message on stan-users.
#> Inference for Stan model: engines.
#> 1 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=1000.
#> 
#>              mean se_mean   sd  2.5%   25%   50%  75% 97.5% n_eff Rhat
#> alpha        3.51    0.00 0.07  3.38  3.45  3.50 3.55  3.66   413 1.01
#> beta[1]      0.55    0.00 0.07  0.43  0.51  0.55 0.59  0.70   462 1.00
#> sigma        0.31    0.00 0.07  0.21  0.26  0.30 0.34  0.46   316 1.00
#> mu_obs[1]    3.75    0.00 0.09  3.61  3.69  3.74 3.80  3.95   371 1.01
#> mu_obs[2]    3.75    0.00 0.09  3.61  3.69  3.74 3.80  3.95   371 1.01
#> mu_obs[3]    3.75    0.00 0.09  3.61  3.69  3.74 3.80  3.95   371 1.01
#> mu_obs[4]    3.75    0.00 0.09  3.61  3.69  3.74 3.80  3.95   371 1.01
#> mu_obs[5]    3.75    0.00 0.09  3.61  3.69  3.74 3.80  3.95   371 1.01
#> mu_obs[6]    3.75    0.00 0.09  3.61  3.69  3.74 3.80  3.95   371 1.01
#> mu_obs[7]    3.75    0.00 0.09  3.61  3.69  3.74 3.80  3.95   371 1.01
#> mu_obs[8]    3.32    0.00 0.07  3.20  3.27  3.31 3.36  3.47   494 1.00
#> mu_obs[9]    3.32    0.00 0.07  3.20  3.27  3.31 3.36  3.47   494 1.00
#> mu_obs[10]   3.32    0.00 0.07  3.20  3.27  3.31 3.36  3.47   494 1.00
#> mu_obs[11]   3.32    0.00 0.07  3.20  3.27  3.31 3.36  3.47   494 1.00
#> mu_obs[12]   3.32    0.00 0.07  3.20  3.27  3.31 3.36  3.47   494 1.00
#> mu_obs[13]   2.74    0.00 0.10  2.54  2.67  2.73 2.79  2.96   660 1.00
#> mu_obs[14]   2.74    0.00 0.10  2.54  2.67  2.73 2.79  2.96   660 1.00
#> mu_obs[15]   2.74    0.00 0.10  2.54  2.67  2.73 2.79  2.96   660 1.00
#> mu_obs[16]   2.74    0.00 0.10  2.54  2.67  2.73 2.79  2.96   660 1.00
#> mu_obs[17]   2.74    0.00 0.10  2.54  2.67  2.73 2.79  2.96   660 1.00
#> mu_cens[1]   4.22    0.01 0.13  4.00  4.13  4.21 4.29  4.52   369 1.01
#> mu_cens[2]   4.22    0.01 0.13  4.00  4.13  4.21 4.29  4.52   369 1.01
#> mu_cens[3]   4.22    0.01 0.13  4.00  4.13  4.21 4.29  4.52   369 1.01
#> mu_cens[4]   4.22    0.01 0.13  4.00  4.13  4.21 4.29  4.52   369 1.01
#> mu_cens[5]   4.22    0.01 0.13  4.00  4.13  4.21 4.29  4.52   369 1.01
#> mu_cens[6]   4.22    0.01 0.13  4.00  4.13  4.21 4.29  4.52   369 1.01
#> mu_cens[7]   4.22    0.01 0.13  4.00  4.13  4.21 4.29  4.52   369 1.01
#> mu_cens[8]   4.22    0.01 0.13  4.00  4.13  4.21 4.29  4.52   369 1.01
#> mu_cens[9]   4.22    0.01 0.13  4.00  4.13  4.21 4.29  4.52   369 1.01
#> mu_cens[10]  4.22    0.01 0.13  4.00  4.13  4.21 4.29  4.52   369 1.01
#> mu_cens[11]  3.75    0.00 0.09  3.61  3.69  3.74 3.80  3.95   371 1.01
#> mu_cens[12]  3.75    0.00 0.09  3.61  3.69  3.74 3.80  3.95   371 1.01
#> mu_cens[13]  3.75    0.00 0.09  3.61  3.69  3.74 3.80  3.95   371 1.01
#> mu_cens[14]  3.32    0.00 0.07  3.20  3.27  3.31 3.36  3.47   494 1.00
#> mu_cens[15]  3.32    0.00 0.07  3.20  3.27  3.31 3.36  3.47   494 1.00
#> mu_cens[16]  3.32    0.00 0.07  3.20  3.27  3.31 3.36  3.47   494 1.00
#> mu_cens[17]  3.32    0.00 0.07  3.20  3.27  3.31 3.36  3.47   494 1.00
#> mu_cens[18]  3.32    0.00 0.07  3.20  3.27  3.31 3.36  3.47   494 1.00
#> mu_cens[19]  2.74    0.00 0.10  2.54  2.67  2.73 2.79  2.96   660 1.00
#> mu_cens[20]  2.74    0.00 0.10  2.54  2.67  2.73 2.79  2.96   660 1.00
#> mu_cens[21]  2.74    0.00 0.10  2.54  2.67  2.73 2.79  2.96   660 1.00
#> mu_cens[22]  2.74    0.00 0.10  2.54  2.67  2.73 2.79  2.96   660 1.00
#> mu_cens[23]  2.74    0.00 0.10  2.54  2.67  2.73 2.79  2.96   660 1.00
#> lp__        -0.46    0.08 1.40 -3.89 -1.09 -0.09 0.56  1.08   325 1.00
#> 
#> Samples were drawn using NUTS(diag_e) at Wed Jun  7 20:47:27 2017.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```


[^engines-src]: This example is derived from Simon Jackman, "Engines: right-censored failure times - the I(,) construct contrasted with other approaches", 2007-07-24,
[URL](https://web-beta.archive.org/web/20070724034205/http://jackman.stanford.edu:80/mcmc/engines.odc)

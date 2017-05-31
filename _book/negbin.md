
# Negative Binomial


```r
library("tidyverse")
library("rstan")
library("rstanarm")
```


The data  are  from the 1990 United States Census for the city of St. Louis, Missouri for Census Tracts, and from records of the St. Louis City Metropolitan Police Department for the years 1980 through 1994. For each Census Tract (with a population), N=111, an observation includes
- the median household income in 1990
- the percentage unemployed (base of labor force)
- a count of the number of homicide incidents.

The number of homicides in this 15 year period totals 2,815.  The average size of a Census Tract is 3,571 with a range of 249--8,791. Income has been rescaled by dividing by 1,000 which produces a range similar to that of percentage unemeployed and standard deviations that are very close.  Tract homicide counts range from 0 through 99 with a median of 16 (mean is 25.+).  An enhanced set of linear, predictors does better than this two predictor example.

$$
\begin{aligned}[t]
y_i &\sim \mathsf{NegBinomial2}(\mu_i,\phi) \\
\mu_i &= \frac{1}{1 + e^{-\eta_i}} \\
\eta_i &= x_i \beta
\end{aligned}
$$
The negative binomial distribution is parameterized so that $\mu \in \mathbb{R}^+$ is the location parameter, and $\phi \in \mathbb{R}^+$ is the reciprocal overdisperson parameter, such that the mean and variance of a random variable $Y$ distributed negative binomial is
$$
\begin{aligned}[t]
E[Y] &= \mu , \\
V[Y] &= \mu + \frac{\mu^2}{\phi} .
\end{aligned}
$$
As $\phi \to \infty$, the negative binomial approaches the Poisson distribution.

The parameters are given weakly informative priors,
$$
\begin{aligned}[t]
\alpha &\sim \mathsf{Normal}(0, 10), \\
\beta_k &\sim \mathsf{Normal}(0, 2.5), \\
\phi^{-1} &\sim \mathsf{HalfCauchy}(0, 5).
\end{aligned}
$$


```r
negbin_mod <- stan_model("stan/negbin.stan")
```
<pre>
  <code class="stan">data {
  int N;
  int y[N];
  int K;
  matrix[N, K] X;
  // priors
  real alpha_mean;
  real<lower = 0.> alpha_scale;
  vector[K] beta_mean;
  vector<lower = 0.>[K] beta_scale;
  real<lower = 0.> reciprocal_phi_scale;
}
parameters {
  real alpha;
  vector[K] beta;
  real<lower = 0.> reciprocal_phi;
}
transformed parameters {
  vector[N] eta;
  real<lower = 0.> phi;
  eta = alpha + X * beta;
  phi = 1. / reciprocal_phi;
}
model {
  reciprocal_phi ~ cauchy(0., reciprocal_phi_scale);
  alpha ~ normal(alpha_mean, alpha_scale);
  beta ~ normal(beta_mean, beta_scale);
  y ~ neg_binomial_2_log(eta, phi);
}
generated quantities {
  vector[N] mu;
  vector[N] log_lik;
  vector[N] y_rep;
  mu = exp(eta);
  for (i in 1:N) {
    log_lik[i] = neg_binomial_2_log_lpmf(y[i] | eta[i], phi);
    y_rep[i] = neg_binomial_2_rng(mu[i], phi);
  }
}</code>
</pre>


```r
load("data/st_louis_census.rda")
negbin_data <- within(list(), {
  y <- st_louis_census$i8094
  N <- length(y)
  X <- model.matrix(~ 0 + pcunemp9 + incrs, data = st_louis_census) %>% scale()
  K <- ncol(X)
  beta_mean <- rep(0, K)
  beta_scale <- rep(2.5, K)  
  alpha_mean <- 0
  alpha_scale <- 10
  reciprocal_phi_scale <- 5
})
```


```r
negbin_fit <- sampling(negbin_mod, data = negbin_data)
```

```r
summary(negbin_fit, par = c("alpha", "beta", "phi"))$summary
#>           mean se_mean     sd   2.5%    25%    50%    75%  97.5% n_eff
#> alpha    2.927 0.00126 0.0722  2.786  2.879  2.925  2.975  3.070  3300
#> beta[1]  0.692 0.00193 0.1100  0.478  0.618  0.692  0.763  0.913  3253
#> beta[2] -0.347 0.00171 0.1015 -0.550 -0.414 -0.346 -0.279 -0.151  3539
#> phi      1.967 0.00584 0.3143  1.428  1.752  1.943  2.154  2.666  2894
#>          Rhat
#> alpha   1.000
#> beta[1] 1.000
#> beta[2] 0.999
#> phi     1.001
```

We could also fit the model using the **[rstanarm](https://cran.r-project.org/package=rstanarm)** function [rstanarm](https://www.rdocumentation.org/packages/rstanarm/topics/stan_glm.nb),

```r
negbin_fit2 <- stan_glm.nb(i8094 ~ pcunemp9 + incrs, data = st_louis_census)
```

```r
negbin_fit2
#> stan_glm.nb
#>  family:  neg_binomial_2 [log]
#>  formula: i8094 ~ pcunemp9 + incrs
#> ------
#> 
#> Estimates:
#>                       Median MAD_SD
#> (Intercept)            2.8    0.4  
#> pcunemp9               0.1    0.0  
#> incrs                 -0.1    0.0  
#> reciprocal_dispersion  2.0    0.3  
#> 
#> Sample avg. posterior predictive 
#> distribution of y (X = xbar):
#>          Median MAD_SD
#> mean_PPD 32.0    5.5  
#> 
#> ------
#> For info on the priors used see help('prior_summary.stanreg').
```

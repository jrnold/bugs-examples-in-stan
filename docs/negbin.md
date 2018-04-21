
# Negative Binomial: Estimating Homicides in Census Tracks {#negbin}


```r
library("tidyverse")
library("rstan")
library("rstanarm")
```

The data  are  from the 1990 United States Census for the city of St. Louis,
Missouri for Census Tracts, and from records of the St. Louis City Metropolitan
Police Department for the years 1980 through 1994. For each Census Tract (with
a population), N=111, an observation includes

-   the median household income in 1990
-   the percentage unemployed (base of labor force)
-   a count of the number of homicide incidents.

The number of homicides in this 15 year period totals 2,815.  The average size
of a Census Tract is 3,571 with a range of 249--8,791. Income has been rescaled
by dividing by 1,000 which produces a range similar to that of percentage
unemployed and standard deviations that are very close.  Tract homicide counts
range from 0 through 99 with a median of 16 (mean is 25.+).  An enhanced set of
linear, predictors does better than this two predictor example.

$$
\begin{aligned}[t]
y_i &\sim \mathsf{NegBinomial2}(\mu_i,\phi) \\
\mu_i &= \frac{1}{1 + e^{-\eta_i}} \\
\eta_i &= x_i \beta
\end{aligned}
$$
The negative binomial distribution is parameterized so that $\mu \in \mathbb{R}^+$ is the location parameter, and $\phi \in \mathbb{R}^+$ is the reciprocal overdispersion parameter, such that the mean and variance of a random variable $Y$ distributed negative binomial is
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
data("st_louis_census", package = "bayesjackman")
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
#> alpha    2.926 0.00114 0.0719  2.787  2.878  2.926  2.973  3.064  4000
#> beta[1]  0.691 0.00197 0.1122  0.471  0.615  0.689  0.766  0.912  3245
#> beta[2] -0.349 0.00171 0.1007 -0.551 -0.415 -0.348 -0.280 -0.154  3481
#> phi      1.968 0.00488 0.3088  1.424  1.751  1.949  2.166  2.639  4000
#>         Rhat
#> alpha      1
#> beta[1]    1
#> beta[2]    1
#> phi        1
```

We could also fit the model using the **rstanarm** function `stan_glm.nb` (or `stan_glm`):

```r
negbin_fit2 <- stan_glm.nb(i8094 ~ pcunemp9 + incrs, data = st_louis_census)
```

```r
negbin_fit2
#> stan_glm.nb
#>  family:       neg_binomial_2 [log]
#>  formula:      i8094 ~ pcunemp9 + incrs
#>  observations: 111
#>  predictors:   3
#> ------
#>                       Median MAD_SD
#> (Intercept)            2.8    0.4  
#> pcunemp9               0.1    0.0  
#> incrs                 -0.1    0.0  
#> reciprocal_dispersion  1.9    0.3  
#> 
#> Sample avg. posterior predictive distribution of y:
#>          Median MAD_SD
#> mean_PPD 32.4    5.6  
#> 
#> ------
#> For info on the priors used see help('prior_summary.stanreg').
```

Example derived from Simon Jackman, "negative binomial using the ones trick with log link", 2005-10-27, [URL](https://web-beta.archive.org/web/20051027082311/http://jackman.stanford.edu:80/mcmc/negbineg.odc).

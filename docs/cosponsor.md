
# Cosponsorship: computing auxiliary quantities from MCMC output {#cosponsorship}


```r
library("tidyverse")
library("rstan")
library("rstanarm")
```

Typically, MCMC output consists of samples from the posterior density of model
parameters.  But note that other quantities can be generated as well, say,
imputations for missing data points, predictions, residuals, or goodness-of-fit
summary statistics.  In fact, any function of the parameters can be calculated
and output.

I demonstrate these ideas in the context of a generalized linear model for
binary data.  The specific application is @Krehbiel1995a, a study of
legislative behavior. [^cosponsor-data] [^cosponsor-src] The response variable
is a binary indicator, which is equal to 1 if the member of the U.S. House of
Representatives chose to cosponsor bill [H.R.
3266](https://www.congress.gov/bill/103rd-congress/house-bill/3266), and 0
otherwise. Of the 434 legislators for which data is available, 228 legislators
cosponsored this bill. H.R.3266 was a wide-ranging spending bill designed to
circumvent the usual budget-making process, that was considered by the [103rd
House of
Representatives](https://en.wikipedia.org/wiki/103rd_United_States_Congress) in
1993--1994. Seven covariates are used in the analysis:

-   liberalism as measured by the interest group Americans for Democratic Action (ADA)
-   fiscal conservatism published by the National Taxpayers' Union (NTU)
-   Democratic Party membership
-   Congressional seniority, measured by years since first election
-   the electoral margin of the member
-   membership of the House Appropriations Committee
-   membership of the House Budget Committee.

@Krehbiel1995a finds that, conditional on legislators' policy preferences, as measured with the ADA and NTU scores, Democrats were more likely to support H.R. 3266 than Republicans.
Seniority is also a key predictor; junior were members more likely to cosponsor this legislation, conditional on the other covariates.

## Model

$$
\begin{aligned}[t]
y_i &= \mathsf{Bernoulli}(\mu_i) \\
\mu_i &= \mathsf{Logit}^{-1}(x_i \beta)
\end{aligned}
$$

Several auxiliary quantities are estimated in the following Stan program:
the *percent correctly predicted (PCP)*,
$$
\mathrm{PCP} = \frac{1}{N} \sum_{i = 1}^N \left( y_i (\pi_i >= 0.5) + (1 - y_i) (\pi_i < 0.5) \right) ,
$$
and the expected percent correctly predicted (ePCP)* [@Herron1999a],
$$
\mathrm{PCP} = \frac{1}{N} \sum_{i = 1}^N \left( y_i \pi_i + (1 - y_i) (1 - \pi_i) \right) .
$$

Uncertainty in the parameters, and propagates to uncertainty in these quantities, and thus it is easy to calculate posterior distributions for these values.


```r
data("a2z", package = "bayesjackman")
```


```r
cosponsor_formula <- cosp ~ ada + ntu + democrat + firstelected + margin + appromember + budgetmember
```


```r
cosponsor_data <- within(list(), {
  y <- a2z$cosp
  N <- length(y)
  X <- model.matrix(update(cosponsor_formula, . ~ 0 + .), data = a2z) %>% scale()
  K <- ncol(X)
  alpha_loc <- 0
  alpha_scale <- 10
  beta_loc <- rep(0, K)
  beta_scale <- rep(2.5, K)
})
```


```r
mod_logit2 <- stan_model("stan/logit2.stan")
```
<pre>
  <code class="stan">functions {
  real pct_correct_pred(int[] y, vector mu) {
    real out;
    int N;
    N = num_elements(mu);
    out = 0.;
    for (i in 1:N) {
      if (y[i]) {
        out = out + int_step(mu[i] >= 0.5);
      } else {
        out = out + int_step(mu[i] < 0.5);
      }
    }
    out = out / N;
    return out;
  }
  real expected_pct_correct_pred(int[] y, vector mu) {
    real out;
    int N;
    N = num_elements(mu);
    out = 0.;
    for (i in 1:N) {
      if (y[i]) {
        out = out + mu[i];
      } else {
        out = out + (1. - mu[i]);
      }
    }
    out = out / N;
    return out;
  }
}
data {
  // response
  int N;
  int y[N];
  // covariates
  int K;
  matrix[N, K] X;
  // priors
  real alpha_loc;
  real<lower = 0.> alpha_scale;
  vector[K] beta_loc;
  vector<lower = 0.>[K] beta_scale;
}
parameters {
  real alpha;
  vector[K] beta;
}
transformed parameters {
  // linear predictor
  vector[N] eta;
  eta = alpha + X * beta;
}
model {
  alpha ~ normal(alpha_loc, alpha_scale);
  beta ~ normal(beta_loc, beta_scale);
  // y ~ bernoulli(inv_logit(eta));
  // this is faster and more numerically stable
  y ~ bernoulli_logit(eta);
}
generated quantities {
  // log-likelihood of each obs
  vector[N] log_lik;
  // probability
  vector[N] mu;
  // simulated outcomes
  int y_rep[N];
  // percent correctly predicted
  real PCP;
  real ePCP;
  for (i in 1:N) {
    mu[i] = inv_logit(eta[i]);
    log_lik[i] = bernoulli_logit_lpmf(y[i] | eta[i]);
    y_rep[i] = bernoulli_rng(mu[i]);
  }
  PCP = pct_correct_pred(y, mu);
  ePCP = expected_pct_correct_pred(y, mu);
}</code>
</pre>

Sample from this model:

```r
cosponsor_fit <- sampling(mod_logit2, data = cosponsor_data)
```

```r
summary(cosponsor_fit, par = c("beta", "alpha", "PCP", "ePCP"))$summary
#>           mean  se_mean      sd   2.5%    25%    50%     75%  97.5% n_eff
#> beta[1] -0.467 5.40e-03 0.34131 -1.134 -0.700 -0.464 -0.2323 0.1901  4000
#> beta[2]  4.205 1.15e-02 0.66477  2.938  3.751  4.192  4.6482 5.5554  3341
#> beta[3]  0.492 7.26e-03 0.45902 -0.431  0.195  0.496  0.8074 1.3767  4000
#> beta[4]  0.874 4.00e-03 0.25293  0.400  0.704  0.870  1.0411 1.3858  4000
#> beta[5] -0.115 3.30e-03 0.20876 -0.532 -0.253 -0.117  0.0226 0.2971  4000
#> beta[6] -0.321 4.26e-03 0.26969 -0.890 -0.498 -0.309 -0.1286 0.1749  4000
#> beta[7] -0.342 3.63e-03 0.22944 -0.807 -0.493 -0.341 -0.1802 0.0977  4000
#> alpha    1.059 4.31e-03 0.27239  0.552  0.874  1.048  1.2322 1.6194  4000
#> PCP      0.909 9.52e-05 0.00544  0.899  0.906  0.910  0.9124 0.9194  3267
#> ePCP     0.873 1.15e-04 0.00725  0.858  0.869  0.874  0.8786 0.8859  4000
#>          Rhat
#> beta[1] 1.000
#> beta[2] 1.000
#> beta[3] 0.999
#> beta[4] 1.000
#> beta[5] 0.999
#> beta[6] 0.999
#> beta[7] 0.999
#> alpha   1.000
#> PCP     1.000
#> ePCP    1.000
```

This model can be fit using the **rstanarm** function `stan_glm`.
This is a binomial model, so it uses `family = binomial()`:


```r
cosponsor_fit2 <- stan_glm(cosponsor_formula, data = a2z, family = binomial())
```


```r
cosponsor_fit2
#> stan_glm
#>  family:       binomial [logit]
#>  formula:      cosp ~ ada + ntu + democrat + firstelected + margin + appromember + 
#> 	   budgetmember
#>  observations: 434
#>  predictors:   8
#> ------
#>              Median MAD_SD
#> (Intercept)  -14.4    3.1 
#> ada            0.0    0.0 
#> ntu            0.2    0.0 
#> democrat       0.9    0.9 
#> firstelected   0.1    0.0 
#> margin         0.0    0.0 
#> appromember   -0.9    0.7 
#> budgetmember  -1.0    0.7 
#> 
#> Sample avg. posterior predictive distribution of y:
#>          Median MAD_SD
#> mean_PPD 0.5    0.0   
#> 
#> ------
#> For info on the priors used see help('prior_summary.stanreg').
```

-   probability of co-sponsorship for a representative member (median x values)
    ```
  	## clarify calculations
    ## probability of co-sponsorship for "average" member (median x vals)
  	probit(pbar) <- beta[1] + beta[2]*0.55 + beta[3]*0.32
  	              + beta[4] + beta[5]*86 + beta[6]*0.26;
  	```

-   party affiliation - holding other values at their representative values
    ```
  	## party affiliation, effect size
  	probit(p.dem) <- beta[1] + beta[2]*0.55 + beta[3]*0.32
  	               + beta[4] + beta[5]*86 + beta[6]*0.26;
    probit(p.rep) <- beta[1] + beta[2]*0.55 + beta[3]*0.32
  	               + beta[5]*86 + beta[6]*0.26;
  	d.party <- p.dem - p.rep;
  	```

-   party *attributable effect* from all sources
    ```
    ## "attributable effect", all sources, due to party affiliation
  	probit(p.dem.all) <- beta[1] + beta[2]*0.8 + beta[3]*0.2
  	                   + beta[4] + beta[5]*86 + beta[6]*0.28;
  	probit(p.rep.all) <- beta[1] + beta[2]*0.1 + beta[3]*0.74
  	                             + beta[5]*86 + beta[6]*0.24;
  	d.party.all <- p.dem.all - p.rep.all
    ```

TODO: Calculate latent residuals [@GelmanGoegebeurTuerlinckxEtAl2000a; @AlbertChib1995a]. Do they make sense and have any meaning outside of Gibbs sampling? They don't seem useful relative to LOO measures of outliers.

[^cosponsor-data]: Replication data from @Herron2010a.

[^cosponsor-src]: Example derived from Simon Jackman, "[Cosponsorship:  computing auxiliary quantities from MCMC output](https://web-beta.archive.org/web/20040501194620/http://jackman.stanford.edu:80/mcmc/kk.odc)", 2004-05-01.

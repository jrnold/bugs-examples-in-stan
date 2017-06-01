
# Bimodal: Extreme missingness in bivariate normal data {#bimodal}


```r
library("rstan")
#> Loading required package: ggplot2
#> Loading required package: StanHeaders
#> rstan (Version 2.15.1, packaged: 2017-04-19 05:03:57 UTC, GitRev: 2e1f913d3ca3)
#> For execution on a local, multicore CPU with excess RAM we recommend calling
#> rstan_options(auto_write = TRUE)
#> options(mc.cores = parallel::detectCores())
library("tidyverse")
#> Loading tidyverse: tibble
#> Loading tidyverse: tidyr
#> Loading tidyverse: readr
#> Loading tidyverse: purrr
#> Loading tidyverse: dplyr
#> Conflicts with tidy packages ----------------------------------------------
#> extract(): tidyr, rstan
#> filter():  dplyr, stats
#> lag():     dplyr, stats
library("stringr")
```


Simple methods for dealing with missing data can run into trouble given pernicious patterns of missingness.  A famous artificial data set designed to highlight this point was created by Gordon Murray, to show how an EM algorithm can run into problems [@Murray1977a,@DempsterLairdRubin1977a].

```
x1:	1	1	-1	-1	2	2	-2 -2	*	*	*	*
x2:	1	-1	1	-1	*	*	*	*	2	2	-2	-2
```

Assume bivariate normality, and that the means of the two variables are both zero, but the variances and covariance are unknown.  Inference about the correlation coefficient  $r$  between these two variables is not trivial in this instance.  The marginal complete-data likelihood for $r$  is not unimodal, and has a saddle-point at zero, and two local maxima close to -1 and 1.  A Bayesian analysis (with uninformative priors) similarly recovers a bimodal posterior density for the correlation coefficient; e.g.,
[@Tanner1996a, @Congdon2007a].


```r
bimodal_mod <- stan_model("stan/bimodal.stan")
```
<pre>
  <code class="stan">data {
  int N;
  int<lower = 0, upper=N> N_obs;
  int<lower = 0, upper=N> N_miss;
  vector[N_obs] x_obs;
  int<lower = 1, upper = N> x_obs_idx[N_obs, 2];
  int<lower = 1, upper = N> x_miss_idx[N_miss, 2];
  vector[2] mu;
}
parameters {
  cov_matrix[2] Sigma;
  vector[N_miss] x_miss;
}
transformed parameters {
  // using an array of vectors is more convenient when sampling
  // multi_normal than using an matrix
  vector[2] X[N];
  for (i in 1:N_obs) {
    X[x_obs_idx[i, 1], x_obs_idx[i, 2]] = x_obs[i];
  }
  for (i in 1:N_miss) {
    X[x_miss_idx[i, 1], x_miss_idx[i, 2]] = x_miss[i];
  }
}
model{
  for (i in 1:N) {
    X[i] ~ multi_normal(mu, Sigma);
  }
}</code>
</pre>

You can ignore the **rstan** warning,
```
DIAGNOSTIC(S) FROM PARSER:
Warning (non-fatal):
Left-hand side of sampling statement (~) may contain a non-linear transform of a parameter or local variable.
If it does, you need to include a target += statement with the log absolute determinant of the Jacobian of the transform.
Left-hand-side of sampling statement:
    X[i] ~ multi_normal(...)
```
since the left hand side is a simple linear relationship and no 
Jacobian adjustment is needed.
All we did was replace index values in the transformed parameter.


```r
X_mat <- matrix(c(1, 1, -1, -1, 2, 2, -2, -2, NA, NA, NA, NA,
       1, -1, 1, -1, NA, NA, NA, NA, 2, 2, -2, -2), ncol = 2)
X_mat <- matrix(rnorm(12), ncol = 2)
X_mat[1, 1] <- NA
X_mat[3, 2] <- NA
#       1, -1, 1, -1, NA, NA, NA, NA, 2, 2, -2, -2), ncol = 2)
X <- X_mat %>%
  as_data_frame() %>%
  mutate(.row = row_number()) %>%
  gather(.col, value, -.row) %>%
  mutate(.col = as.integer(str_replace(.col, "V", "")))

X_obs <- filter(X, !is.na(value))
X_miss <- filter(X, is.na(value))
  
bimodal_data <- within(list(), {
  N <- nrow(X)
  x_obs <- X_obs$value
  x_obs_idx <- as.matrix(X_obs[ , c(".row", ".col")])
  N_obs <- nrow(X_obs)
  x_miss_idx <- as.matrix(X_miss[ , c(".row", ".col")])
  N_miss <- nrow(X_miss)
})
```


```r
bimodal_fit <- sampling(bimodal_mod, data = bimodal_data,
                        chains = 1)
#> Warning in is.na(x): is.na() applied to non-(list or vector) of type 'NULL'
#> Warning in FUN(X[[i]], ...): data with name mu is not numeric and not used
```

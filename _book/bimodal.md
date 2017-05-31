
# Bimodal: Extreme missingness in bivariate normal data


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


Simple methods for dealing with missing data can run into trouble given pernicious patterns of missingness.  A famous artificial data set designed to highlight this point was created by Gordon Murray, to show how an EM algorithm can run into problems (see the Journal of the Royal Statistical Society Series B, 39:27, 1977; this example appears in the discussion to Dempster, Laird and Rubin's much-cited EM paper):

```
x1:	1	1	-1	-1	2	2	-2 -2	*	*	*	*
x2:	1	-1	1	-1	*	*	*	*	2	2	-2	-2
```

Assume bivariate normality, and that the means of the two variables are both zero, but the variances and covariance are unknown.  Inference about the correlation coefficient  $r$  between these two variables is not trivial in this instance.  The marginal complete-data likelihood for $r$  is not unimodal, and has a saddle-point at zero, and two local maxima close to -1 and 1.  A Bayesian analysis (with uninformative priors) similarly recovers a bimodal posterior density for the correlation coefficient; e.g., see Tanner, Tools for Statistical Inference, 3rd edition, pp95-96 or Congdon, Bayesian Statistical Modelling, p46.


```r
bimodal_mod <- stan_model("stan/bimodal.stan")
```

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
X <- matrix(c(1, 1, -1, -1, 2, 2, -2, -2, NA, NA, NA, NA,
       1, -1, 1, -1, NA, NA, NA, NA, 2, 2, -2, -2), ncol = 2) %>%
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
  mu <- rep(0, 2)
  Sigma_cov <- diag(2)
  Sigma_df <- 2
})
```


```r
bimodal_fit <- sampling(bimodal_mod, data = bimodal_data,
                        chains = 1)
```


# Multivariate Missing Data {#multivarmissing}

$$
\DeclareMathOperator{diag}{diag}
$$

```r
library("tidyverse")
library("rstan")
```

This example shows how to impute missing data. See @Stan2016a, Chapter 10 "Missing Data & Partially Known Parameters" for more discussion.[^multivarmissing-src]

Consider a data set of 10 observations on 3 variables
Only one of the variables, $z$, is completely observed.
The other two variables, x$ and $y$, have a non-overlapping pattern of missing data.


```r
multivarmissing <- tribble(
  ~x, ~y, ~z,
1,	NA,	NA,
2,	NA,	4,
3,	NA,	3,
4,	NA,	5,
5,	NA,	7,
NA,	7,	9,
NA,	8,	8,
NA,	9,	11,
NA,	8,	10,
NA,	9,	8)
```

The missing elements of $x$ and $y$ are parameters, and the observed elements of $x$, $y$, and $z$ are data.
These are combined in the `transformed parameters` block, and modeled.

## Separate Regressions

We use $z$ to predict $x$,
and $z$ and $x$ (both observed and imputed) to impute $y$.

$$
\begin{aligned}[t]
x_i &\sim \mathsf{Normal}(\mu_{x,i}, \sigma_x) \\
\mu_{x,i} &= \gamma_1 + \gamma_2 z_i \\
y_i &\sim \mathsf{Normal}(\mu_{y,i}, \sigma_y) \\
\mu_{y,i} &= \beta_1 + \beta_2 y_i + \beta_3 z_i \\
z_i &\sim \mathsf{Normal}(\mu_z, \sigma_z)
\end{aligned}
$$

The parameters are given weakly informative parameters:
$$
\begin{aligned}[t]
\sigma_x,\sigma_y,\sigma_z &\sim \mathsf{HalfCauchy}(0, 5) \\
\gamma_1, \beta_1 &\sim \mathsf{Normal}(0, 10) \\
\gamma_2, \beta_2, \beta_3 &\sim \mathsf{Normal}(0, 2.5)
\end{aligned}
$$
Note that this assumes that $x$, $y$, and $z$ are standardized to have zero mean and unit variance. 


```r
data_multivarmissing <- within(list(), {
  N <- nrow(multivarmissing)
  x_obs <- multivarmissing$x[!is.na(multivarmissing$x)] %>%
    scale() %>% as.numeric()
  x_obs_idx <- array(which(!is.na(multivarmissing$x)))
  N_x_obs <- length(x_obs_idx)  
  x_miss_idx <- array(which(is.na(multivarmissing$x)))
  N_x_miss <- length(x_miss_idx)
  y_obs <- multivarmissing$y[!is.na(multivarmissing$y)] %>%
    scale() %>% as.numeric()    
  y_obs_idx <- array(which(!is.na(multivarmissing$y)))
  N_y_obs <- length(y_obs_idx)  
  y_miss_idx <- array(which(is.na(multivarmissing$y)))
  N_y_miss <- length(y_miss_idx)
  z_obs <- multivarmissing$z[!is.na(multivarmissing$z)] %>%
    scale() %>% as.numeric()
  z_obs_idx <- array(which(!is.na(multivarmissing$z)))
  N_z_obs <- length(z_obs_idx)
  z_miss_idx <- array(which(is.na(multivarmissing$z)))
  N_z_miss <- length(z_miss_idx)
  alpha_loc <- 0
  alpha_scale <- 10
  beta_loc <- rep(0, 3)
  beta_scale <- c(10, 2.5, 2.5)
  gamma_loc <- rep(0, 2)
  gamma_scale <- c(10, 2.5)
  sigma_x_scale <- 5
  sigma_y_scale <- 5
  sigma_z_scale <- 5
})
```


```r
mod_multivarmissing <- stan_model("stan/multivarmissing2.stan")
#> In file included from filef7896b8e0b75.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:12:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/gevv_vvv_vari.hpp:5:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/var.hpp:7:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/math/tools/config.hpp:13:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/config.hpp:39:
#> /Users/jrnold/Library/R/3.4/library/BH/include/boost/config/compiler/clang.hpp:196:11: warning: 'BOOST_NO_CXX11_RVALUE_REFERENCES' macro redefined [-Wmacro-redefined]
#> #  define BOOST_NO_CXX11_RVALUE_REFERENCES
#>           ^
#> <command line>:6:9: note: previous definition is here
#> #define BOOST_NO_CXX11_RVALUE_REFERENCES 1
#>         ^
#> In file included from filef7896b8e0b75.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:42:
#> /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/set_zero_all_adjoints.hpp:14:17: warning: unused function 'set_zero_all_adjoints' [-Wunused-function]
#>     static void set_zero_all_adjoints() {
#>                 ^
#> In file included from filef7896b8e0b75.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:43:
#> /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/set_zero_all_adjoints_nested.hpp:17:17: warning: 'static' function 'set_zero_all_adjoints_nested' declared in header file should be declared 'static inline' [-Wunneeded-internal-declaration]
#>     static void set_zero_all_adjoints_nested() {
#>                 ^
#> In file included from filef7896b8e0b75.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:11:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat.hpp:59:
#> /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/autocorrelation.hpp:17:14: warning: function 'fft_next_good_size' is not needed and will not be emitted [-Wunneeded-internal-declaration]
#>       size_t fft_next_good_size(size_t N) {
#>              ^
#> In file included from filef7896b8e0b75.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:11:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat.hpp:298:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/arr.hpp:39:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/arr/functor/integrate_ode_rk45.hpp:13:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/numeric/odeint.hpp:61:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/numeric/odeint/util/multi_array_adaption.hpp:29:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/multi_array.hpp:21:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/multi_array/base.hpp:28:
#> /Users/jrnold/Library/R/3.4/library/BH/include/boost/multi_array/concept_checks.hpp:42:43: warning: unused typedef 'index_range' [-Wunused-local-typedef]
#>       typedef typename Array::index_range index_range;
#>                                           ^
#> /Users/jrnold/Library/R/3.4/library/BH/include/boost/multi_array/concept_checks.hpp:43:37: warning: unused typedef 'index' [-Wunused-local-typedef]
#>       typedef typename Array::index index;
#>                                     ^
#> /Users/jrnold/Library/R/3.4/library/BH/include/boost/multi_array/concept_checks.hpp:53:43: warning: unused typedef 'index_range' [-Wunused-local-typedef]
#>       typedef typename Array::index_range index_range;
#>                                           ^
#> /Users/jrnold/Library/R/3.4/library/BH/include/boost/multi_array/concept_checks.hpp:54:37: warning: unused typedef 'index' [-Wunused-local-typedef]
#>       typedef typename Array::index index;
#>                                     ^
#> 8 warnings generated.
```

```r
mod_multivarmissing
```

<pre>
  <code class="stan">data {
  // number of obs
  int<lower = 1> N;
  // number of variables
  int<lower = 1> K;
  // X
  int<lower = 0, upper = N * K> N_obs;
  vector[N_obs] X_obs;
  int<lower = 1, upper = N> X_obs_row[N_obs];
  int<lower = 1, upper = K> X_obs_col[N_obs];
  int<lower = 0, upper = N * K> N_miss;
  int<lower = 1, upper = N> X_miss_row[N_miss];
  int<lower = 1, upper = N> X_miss_col[N_miss];
  // priors
  vector<lower = 0.>[K] Sigma_scale_scale;
  real<lower = 0.> Sigma_corr_L_eta;
  vector[K] mu_loc;
  vector<lower = 0.>[K] mu_scale;
}
parameters {
  vector[K] mu;
  vector<lower = 0.>[K] Sigma_scale;
  cholesky_factor_corr[K] Sigma_corr_L;
  vector[N_miss] X_miss;
}
transformed parameters {
  vector[K] X[N];
  for (i in 1:N_obs) {
   X[X_obs_row[i], X_obs_col[i]] = X_obs[i];
  }
  for (i in 1:N_miss) {
    X[X_miss_row[i], X_miss_col[i]] = X_miss[i];
  }
}
model {
  Sigma_corr_L ~ lkj_corr_cholesky(Sigma_corr_L_eta);
  Sigma_scale ~ cauchy(0., Sigma_scale_scale);
  for (i in 1:N) {
    X[i] ~ multi_normal_cholesky(mu, Sigma_corr_L);
  }
}</code>
</pre>


```r
fit_multivarmissing <- 
  sampling(mod_multivarmissing, data = data_multivarmissing)
#> Warning in is.na(x): is.na() applied to non-(list or vector) of type 'NULL'
#> Warning in FUN(X[[i]], ...): data with name K is not numeric and not used
#> Warning in is.na(x): is.na() applied to non-(list or vector) of type 'NULL'
#> Warning in FUN(X[[i]], ...): data with name N_obs is not numeric and not
#> used
#> Warning in is.na(x): is.na() applied to non-(list or vector) of type 'NULL'
#> Warning in FUN(X[[i]], ...): data with name X_obs is not numeric and not
#> used
#> Warning in is.na(x): is.na() applied to non-(list or vector) of type 'NULL'
#> Warning in FUN(X[[i]], ...): data with name X_obs_row is not numeric and
#> not used
#> Warning in is.na(x): is.na() applied to non-(list or vector) of type 'NULL'
#> Warning in FUN(X[[i]], ...): data with name X_obs_col is not numeric and
#> not used
#> Warning in is.na(x): is.na() applied to non-(list or vector) of type 'NULL'
#> Warning in FUN(X[[i]], ...): data with name N_miss is not numeric and not
#> used
#> Warning in is.na(x): is.na() applied to non-(list or vector) of type 'NULL'
#> Warning in FUN(X[[i]], ...): data with name X_miss_row is not numeric and
#> not used
#> Warning in is.na(x): is.na() applied to non-(list or vector) of type 'NULL'
#> Warning in FUN(X[[i]], ...): data with name X_miss_col is not numeric and
#> not used
#> Warning in is.na(x): is.na() applied to non-(list or vector) of type 'NULL'
#> Warning in FUN(X[[i]], ...): data with name Sigma_scale_scale is not
#> numeric and not used
#> Warning in is.na(x): is.na() applied to non-(list or vector) of type 'NULL'
#> Warning in FUN(X[[i]], ...): data with name Sigma_corr_L_eta is not numeric
#> and not used
#> Warning in is.na(x): is.na() applied to non-(list or vector) of type 'NULL'
#> Warning in FUN(X[[i]], ...): data with name mu_loc is not numeric and not
#> used
#> Warning in is.na(x): is.na() applied to non-(list or vector) of type 'NULL'
#> Warning in FUN(X[[i]], ...): data with name mu_scale is not numeric and not
#> used
#> failed to create the sampler; sampling not done
```

```r
fit_multivarmissing
#> Stan model 'multivarmissing2' does not contain samples.
```

## Multivariate Normal

Alternatively, $x$, $y$, and $z$ could be modeled as coming from a multivariate normal distribution.
$$
\begin{bmatrix}
x_i \\
y_i \\
z_i
\end{bmatrix} \sim 
\mathsf{Normal}(\mu, \Sigma)
$$
where $\mu$ and $\Sigma$ are given weakly informative priors,
$$
\begin{aligned}[t]
\mu_{i,k} &\sim \mathsf{Normal}(0, 10) & k \in \{1, 2, 3\}, \\
\Sigma &= \diag{\sigma} R \diag{sigma}, \\
\sigma &\sim \mathsf{HalfCauchy}(0, 5), \\
R &\sim \mathsf{LkjCorr}(2) .
\end{aligned}
$$


```r
data_multivarmissing2 <- within(list(), {
  N <- nrow(multivarmissing)
  K <- ncol(multivarmissing)
  mu_loc <- rep(0, 3)
  mu_scale <- rep(0, 10)
  Sigma_scale_scale <- 5
  Sigma_corr_L_eta <- 2
})
```


```r
mod_multivarmissing2 <- stan_model("stan/multivarmissing2.stan")
```

```r
mod_multivarmissing2
```

<pre>
  <code class="stan">data {
  // number of obs
  int<lower = 1> N;
  // number of variables
  int<lower = 1> K;
  // X
  int<lower = 0, upper = N * K> N_obs;
  vector[N_obs] X_obs;
  int<lower = 1, upper = N> X_obs_row[N_obs];
  int<lower = 1, upper = K> X_obs_col[N_obs];
  int<lower = 0, upper = N * K> N_miss;
  int<lower = 1, upper = N> X_miss_row[N_miss];
  int<lower = 1, upper = N> X_miss_col[N_miss];
  // priors
  vector<lower = 0.>[K] Sigma_scale_scale;
  real<lower = 0.> Sigma_corr_L_eta;
  vector[K] mu_loc;
  vector<lower = 0.>[K] mu_scale;
}
parameters {
  vector[K] mu;
  vector<lower = 0.>[K] Sigma_scale;
  cholesky_factor_corr[K] Sigma_corr_L;
  vector[N_miss] X_miss;
}
transformed parameters {
  vector[K] X[N];
  for (i in 1:N_obs) {
   X[X_obs_row[i], X_obs_col[i]] = X_obs[i];
  }
  for (i in 1:N_miss) {
    X[X_miss_row[i], X_miss_col[i]] = X_miss[i];
  }
}
model {
  Sigma_corr_L ~ lkj_corr_cholesky(Sigma_corr_L_eta);
  Sigma_scale ~ cauchy(0., Sigma_scale_scale);
  for (i in 1:N) {
    X[i] ~ multi_normal_cholesky(mu, Sigma_corr_L);
  }
}</code>
</pre>


```r
fit_multivarmissing <- 
  sampling(mod_multivarmissing2, data = data_multivarmissing2)
#> Warning in is.na(x): is.na() applied to non-(list or vector) of type 'NULL'
#> Warning in FUN(X[[i]], ...): data with name N_obs is not numeric and not
#> used
#> Warning in is.na(x): is.na() applied to non-(list or vector) of type 'NULL'
#> Warning in FUN(X[[i]], ...): data with name X_obs is not numeric and not
#> used
#> Warning in is.na(x): is.na() applied to non-(list or vector) of type 'NULL'
#> Warning in FUN(X[[i]], ...): data with name X_obs_row is not numeric and
#> not used
#> Warning in is.na(x): is.na() applied to non-(list or vector) of type 'NULL'
#> Warning in FUN(X[[i]], ...): data with name X_obs_col is not numeric and
#> not used
#> Warning in is.na(x): is.na() applied to non-(list or vector) of type 'NULL'
#> Warning in FUN(X[[i]], ...): data with name N_miss is not numeric and not
#> used
#> Warning in is.na(x): is.na() applied to non-(list or vector) of type 'NULL'
#> Warning in FUN(X[[i]], ...): data with name X_miss_row is not numeric and
#> not used
#> Warning in is.na(x): is.na() applied to non-(list or vector) of type 'NULL'
#> Warning in FUN(X[[i]], ...): data with name X_miss_col is not numeric and
#> not used
#> failed to create the sampler; sampling not done
```

```r
fit_multivarmissing
#> Stan model 'multivarmissing2' does not contain samples.
```

[^multivarmissing-src]: This example is derived from Simon Jackan, "Multivariate Missing Data," 2002-06-18, [URL](https://web-beta.archive.org/web/20020618183148/http://jackman.stanford.edu:80/mcmc/multivarmissing.odc).

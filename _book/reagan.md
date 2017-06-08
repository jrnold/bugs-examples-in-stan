
# Reagan: linear regression with AR(1) disturbances {#reagan}


```r
library("tidyverse")
library("rstan")
```

Ninety-six monthly observations on presidential job approval ratings for Ronald Reagan are modeled via linear regression, with a correction for first-order serial correlation among the disturbances.[^reagan]
Note the marginal model for the first observation, and the conditioning on the lagged observation for months 2 through 96.
A uniform prior over the stationary (-1,1) interval is employed for the residual AR(1) parameter.

$$
\begin{aligned}[t]
y_i &= \mu_i + \epsilon_i + \theta \epsilon_{i - 1}  ,\\
\mu_i &= \alpha + x_i' \beta , \\
\epsilon_i &\sim \mathsf{Normal}(0, \sigma^2) ,
\end{aligned}
$$
for $i \in 1, \dots, N$.
Weakly informative priors for each parameter are used,
$$
\begin{aligned}[t]
\alpha &\sim \mathsf{Normal}(0, 10), \\
\beta_k &\sim \mathsf{Normal}(0, 2.5), & k \in 1, \dots, K, \\
\sigma &\sim \mathsf{HalfCauchy}(0, 5), \\
\theta &= 2 \theta^*  - 1 , \\
\theta^* &\sim \mathsf{Beta}(1, 1)  .
\end{aligned}
$$


```r
(load("data/ReaganApproval.rda"))
#> [1] "ReaganApproval"
ReaganApproval
#> # A tibble: 96 x 3
#>     app  infl unemp
#>   <dbl> <dbl> <dbl>
#> 1    51 11.79   7.5
#> 2    55 11.39   7.4
#> 3    60 10.61   7.4
#> 4    67 10.14   7.2
#> 5    68  9.79   7.5
#> 6    59  9.70   7.5
#> # ... with 90 more rows
```



```r
reagan_data <- within(list(), {
  y <- ReaganApproval$app
  N <- length(y)
  X <- model.matrix(~ 0 + infl + unemp, data = ReaganApproval) %>% scale()
  K <- ncol(X)
  alpha_loc <- 0
  alpha_scale <- 10
  beta_loc <- rep(0, K)
  beta_scale <- rep(2.5 * sd(y), K)
  sigma_scale <- 5 * sd(y)
  theta_a <- 1
  theta_b <- 1
})
```


```r
mod_regar1 <- stan_model("stan/regar1.stan")
#> In file included from filedc2be966059.cpp:8:
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
#> In file included from filedc2be966059.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:42:
#> /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/set_zero_all_adjoints.hpp:14:17: warning: unused function 'set_zero_all_adjoints' [-Wunused-function]
#>     static void set_zero_all_adjoints() {
#>                 ^
#> In file included from filedc2be966059.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:43:
#> /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/set_zero_all_adjoints_nested.hpp:17:17: warning: 'static' function 'set_zero_all_adjoints_nested' declared in header file should be declared 'static inline' [-Wunneeded-internal-declaration]
#>     static void set_zero_all_adjoints_nested() {
#>                 ^
#> In file included from filedc2be966059.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:11:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat.hpp:59:
#> /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/autocorrelation.hpp:17:14: warning: function 'fft_next_good_size' is not needed and will not be emitted [-Wunneeded-internal-declaration]
#>       size_t fft_next_good_size(size_t N) {
#>              ^
#> In file included from filedc2be966059.cpp:8:
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
mod_regar1
```

<pre>
  <code class="stan">data {
  // number of observations
  // need at least two to estimates
  int<lower = 2> N;
  // response
  vector[N] y;
  // regression design matrix
  int<lower = 1> K;
  matrix[N, K] X;
  // priors
  real alpha_loc;
  real<lower = 0.> alpha_scale;
  vector[K] beta_loc;
  vector<lower = 0.>[K] beta_scale;
  real<lower = 0.> sigma_scale;
  real<lower = 0.> theta_a;
  real<lower = 0.> theta_b;
}
parameters {
  // regression coefficients
  real alpha;
  vector[K] beta;
  // error scale
  real<lower=0> sigma;
  // lag coefficients
  real<lower = 0, upper = 1> theta_raw;
}
transformed parameters {
  // observation means
  vector[N] mu;
  // error terms
  vector[N] epsilon;
  // lag coefficient;
  real<lower = -1, upper = 1> theta;
  // convert range of theta from (0, 1) to (-1, 1)
  theta = (2. * theta_raw - 1.);
  // regression
  mu = alpha + X * beta;
  // construct errors
  epsilon[1] = y[1] - mu[1];
  for (i in 2:N) {
    epsilon[i] = y[i] - mu[i] - theta * epsilon[i - 1];
  }
}
model {
  alpha ~ cauchy(alpha_loc, alpha_scale);
  beta ~ cauchy(beta_loc, beta_scale);
  theta_raw ~ beta(theta_a, theta_b);
  sigma ~ cauchy(0, sigma_scale);
  for (i in 2:N) {
    y[i] ~ normal(mu[i] + theta * epsilon[i - 1], sigma);
  }
}</code>
</pre>


```r
reagan_fit <- sampling(mod_regar1, data = reagan_data)
```


```r
summary(reagan_fit, par = c("alpha", "beta", "theta", "sigma"))$summary
#>           mean  se_mean     sd   2.5%     25%    50%    75%  97.5% n_eff
#> alpha   53.196 0.013150 0.8317 51.604 52.6292 53.187 53.765 54.801  4000
#> beta[1]  0.642 0.014547 0.9200 -1.145  0.0178  0.641  1.271  2.405  4000
#> beta[2] -3.721 0.013359 0.8449 -5.362 -4.2776 -3.729 -3.175 -2.066  4000
#> theta    0.676 0.000903 0.0571  0.556  0.6398  0.679  0.716  0.777  4000
#> sigma    4.910 0.005838 0.3692  4.252  4.6588  4.886  5.140  5.728  4000
#>          Rhat
#> alpha   0.999
#> beta[1] 1.000
#> beta[2] 1.000
#> theta   1.001
#> sigma   1.000
```


## Cochrane-Orcutt/Prais-Winsten

An AR(1) error model can also be estimated via partial-differencing as 
in Cochrane-Orcutt/Prais-Winsten estimation.

[^reagan]: Example derived from Simon Jackman, "Reagan: linear regression with AR(1) disturbances," *BUGS Examples,* 2007-07-24, [URL](https://web-beta.archive.org/web/20070724034151/http://jackman.stanford.edu:80/mcmc/reagan.odc).

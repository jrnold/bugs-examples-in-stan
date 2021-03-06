
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
\theta &= 2 \theta^{*} - 1 , \\
\theta^{*} &\sim \mathsf{Beta}(1, 1)  .
\end{aligned}
$$


```r
data("ReaganApproval", package = "bayesjackman")
ReaganApproval
#> # A tibble: 96 x 3
#>     app  infl unemp
#>   <dbl> <dbl> <dbl>
#> 1   51. 11.8   7.50
#> 2   55. 11.4   7.40
#> 3   60. 10.6   7.40
#> 4   67. 10.1   7.20
#> 5   68.  9.79  7.50
#> 6   59.  9.70  7.50
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
#> In file included from filed1046087152.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:12:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/gevv_vvv_vari.hpp:5:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/var.hpp:7:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/math/tools/config.hpp:13:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/config.hpp:39:
#> /Users/jrnold/Library/R/3.4/library/BH/include/boost/config/compiler/clang.hpp:200:11: warning: 'BOOST_NO_CXX11_RVALUE_REFERENCES' macro redefined [-Wmacro-redefined]
#> #  define BOOST_NO_CXX11_RVALUE_REFERENCES
#>           ^
#> <command line>:6:9: note: previous definition is here
#> #define BOOST_NO_CXX11_RVALUE_REFERENCES 1
#>         ^
#> In file included from filed1046087152.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:14:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/matrix_vari.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat/fun/Eigen_NumTraits.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Dense:1:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Core:531:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1046087152.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:14:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/matrix_vari.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat/fun/Eigen_NumTraits.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Dense:2:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/LU:47:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1046087152.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:14:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/matrix_vari.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat/fun/Eigen_NumTraits.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Dense:3:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Cholesky:12:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Jacobi:29:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1046087152.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:14:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/matrix_vari.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat/fun/Eigen_NumTraits.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Dense:3:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Cholesky:43:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1046087152.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:14:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/matrix_vari.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat/fun/Eigen_NumTraits.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Dense:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/QR:17:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Householder:27:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1046087152.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:14:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/matrix_vari.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat/fun/Eigen_NumTraits.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Dense:5:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/SVD:48:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1046087152.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:14:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/matrix_vari.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat/fun/Eigen_NumTraits.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Dense:6:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Geometry:58:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1046087152.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:14:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/matrix_vari.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat/fun/Eigen_NumTraits.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Dense:7:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Eigenvalues:58:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1046087152.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:36:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/operator_unary_plus.hpp:7:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/scal/fun/constants.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/math/constants/constants.hpp:13:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/math/tools/convert_from_string.hpp:15:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/lexical_cast.hpp:32:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/lexical_cast/try_lexical_convert.hpp:42:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/lexical_cast/detail/converter_lexical.hpp:52:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/container/container_fwd.hpp:61:
#> /Users/jrnold/Library/R/3.4/library/BH/include/boost/container/detail/std_fwd.hpp:27:1: warning: inline namespaces are a C++11 feature [-Wc++11-inline-namespace]
#> BOOST_MOVE_STD_NS_BEG
#> ^
#> /Users/jrnold/Library/R/3.4/library/BH/include/boost/move/detail/std_ns_begin.hpp:18:34: note: expanded from macro 'BOOST_MOVE_STD_NS_BEG'
#>    #define BOOST_MOVE_STD_NS_BEG _LIBCPP_BEGIN_NAMESPACE_STD
#>                                  ^
#> /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__config:439:52: note: expanded from macro '_LIBCPP_BEGIN_NAMESPACE_STD'
#> #define _LIBCPP_BEGIN_NAMESPACE_STD namespace std {inline namespace _LIBCPP_NAMESPACE {
#>                                                    ^
#> In file included from filed1046087152.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:12:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat.hpp:83:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/csr_extract_u.hpp:6:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Sparse:26:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/SparseCore:66:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1046087152.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:12:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat.hpp:83:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/csr_extract_u.hpp:6:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Sparse:27:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/OrderingMethods:71:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1046087152.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:12:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat.hpp:83:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/csr_extract_u.hpp:6:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Sparse:29:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/SparseCholesky:43:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1046087152.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:12:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat.hpp:83:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/csr_extract_u.hpp:6:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Sparse:32:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/SparseQR:35:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1046087152.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:12:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat.hpp:83:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/csr_extract_u.hpp:6:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Sparse:33:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/IterativeLinearSolvers:46:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1046087152.cpp:686:
#> In file included from /Users/jrnold/Library/R/3.4/library/rstan/include/rstan/rstaninc.hpp:3:
#> In file included from /Users/jrnold/Library/R/3.4/library/rstan/include/rstan/stan_fit.hpp:36:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/services/optimize/bfgs.hpp:11:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/optimization/bfgs.hpp:9:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/optimization/lbfgs_update.hpp:6:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/circular_buffer.hpp:54:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/circular_buffer/details.hpp:20:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/move/move.hpp:30:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/move/iterator.hpp:27:
#> /Users/jrnold/Library/R/3.4/library/BH/include/boost/move/detail/iterator_traits.hpp:29:1: warning: inline namespaces are a C++11 feature [-Wc++11-inline-namespace]
#> BOOST_MOVE_STD_NS_BEG
#> ^
#> /Users/jrnold/Library/R/3.4/library/BH/include/boost/move/detail/std_ns_begin.hpp:18:34: note: expanded from macro 'BOOST_MOVE_STD_NS_BEG'
#>    #define BOOST_MOVE_STD_NS_BEG _LIBCPP_BEGIN_NAMESPACE_STD
#>                                  ^
#> /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__config:439:52: note: expanded from macro '_LIBCPP_BEGIN_NAMESPACE_STD'
#> #define _LIBCPP_BEGIN_NAMESPACE_STD namespace std {inline namespace _LIBCPP_NAMESPACE {
#>                                                    ^
#> In file included from filed1046087152.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:44:
#> /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/set_zero_all_adjoints.hpp:14:17: warning: unused function 'set_zero_all_adjoints' [-Wunused-function]
#>     static void set_zero_all_adjoints() {
#>                 ^
#> In file included from filed1046087152.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:45:
#> /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/set_zero_all_adjoints_nested.hpp:17:17: warning: 'static' function 'set_zero_all_adjoints_nested' declared in header file should be declared 'static inline' [-Wunneeded-internal-declaration]
#>     static void set_zero_all_adjoints_nested() {
#>                 ^
#> In file included from filed1046087152.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:12:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat.hpp:58:
#> /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/autocorrelation.hpp:17:14: warning: function 'fft_next_good_size' is not needed and will not be emitted [-Wunneeded-internal-declaration]
#>       size_t fft_next_good_size(size_t N) {
#>              ^
#> 19 warnings generated.
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
#>           mean  se_mean     sd   2.5%      25%    50%    75%  97.5% n_eff
#> alpha   53.195 0.013180 0.8336 51.576 52.62510 53.199 53.758 54.784  4000
#> beta[1]  0.613 0.014395 0.9104 -1.162  0.00317  0.612  1.240  2.396  4000
#> beta[2] -3.690 0.013022 0.8236 -5.328 -4.23172 -3.695 -3.130 -2.069  4000
#> theta    0.675 0.000881 0.0557  0.555  0.64104  0.678  0.713  0.775  4000
#> sigma    4.906 0.006004 0.3797  4.230  4.63457  4.879  5.153  5.700  4000
#>          Rhat
#> alpha   0.999
#> beta[1] 1.000
#> beta[2] 1.001
#> theta   1.000
#> sigma   1.000
```

## Cochrane-Orcutt/Prais-Winsten

An AR(1) error model can also be estimated [Prais-Winsten](https://en.wikipedia.org/wiki/Prais%E2%80%93Winsten_estimation) estimation:
$$
\begin{aligned}[t]
y_1 &\sim \mathsf{Normal}\left(\alpha + x_1' \beta, \frac{\sigma ^ 2}{1 - \theta ^ 2} \right), \\
y_i &\sim \mathsf{Normal}\left(\theta y_{i - 1} + \alpha (1 - \theta) + \beta (X_i - \theta X_{i - 1}), \sigma ^ 2 \right) & i = 2, \dots, N
\end{aligned}
$$


```r
mod_pw <- stan_model("stan/pw.stan")
#> In file included from filed1062414e11.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:12:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/gevv_vvv_vari.hpp:5:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/var.hpp:7:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/math/tools/config.hpp:13:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/config.hpp:39:
#> /Users/jrnold/Library/R/3.4/library/BH/include/boost/config/compiler/clang.hpp:200:11: warning: 'BOOST_NO_CXX11_RVALUE_REFERENCES' macro redefined [-Wmacro-redefined]
#> #  define BOOST_NO_CXX11_RVALUE_REFERENCES
#>           ^
#> <command line>:6:9: note: previous definition is here
#> #define BOOST_NO_CXX11_RVALUE_REFERENCES 1
#>         ^
#> In file included from filed1062414e11.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:14:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/matrix_vari.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat/fun/Eigen_NumTraits.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Dense:1:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Core:531:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1062414e11.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:14:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/matrix_vari.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat/fun/Eigen_NumTraits.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Dense:2:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/LU:47:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1062414e11.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:14:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/matrix_vari.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat/fun/Eigen_NumTraits.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Dense:3:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Cholesky:12:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Jacobi:29:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1062414e11.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:14:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/matrix_vari.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat/fun/Eigen_NumTraits.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Dense:3:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Cholesky:43:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1062414e11.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:14:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/matrix_vari.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat/fun/Eigen_NumTraits.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Dense:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/QR:17:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Householder:27:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1062414e11.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:14:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/matrix_vari.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat/fun/Eigen_NumTraits.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Dense:5:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/SVD:48:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1062414e11.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:14:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/matrix_vari.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat/fun/Eigen_NumTraits.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Dense:6:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Geometry:58:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1062414e11.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:14:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/matrix_vari.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat/fun/Eigen_NumTraits.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Dense:7:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Eigenvalues:58:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1062414e11.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:36:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/operator_unary_plus.hpp:7:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/scal/fun/constants.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/math/constants/constants.hpp:13:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/math/tools/convert_from_string.hpp:15:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/lexical_cast.hpp:32:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/lexical_cast/try_lexical_convert.hpp:42:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/lexical_cast/detail/converter_lexical.hpp:52:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/container/container_fwd.hpp:61:
#> /Users/jrnold/Library/R/3.4/library/BH/include/boost/container/detail/std_fwd.hpp:27:1: warning: inline namespaces are a C++11 feature [-Wc++11-inline-namespace]
#> BOOST_MOVE_STD_NS_BEG
#> ^
#> /Users/jrnold/Library/R/3.4/library/BH/include/boost/move/detail/std_ns_begin.hpp:18:34: note: expanded from macro 'BOOST_MOVE_STD_NS_BEG'
#>    #define BOOST_MOVE_STD_NS_BEG _LIBCPP_BEGIN_NAMESPACE_STD
#>                                  ^
#> /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__config:439:52: note: expanded from macro '_LIBCPP_BEGIN_NAMESPACE_STD'
#> #define _LIBCPP_BEGIN_NAMESPACE_STD namespace std {inline namespace _LIBCPP_NAMESPACE {
#>                                                    ^
#> In file included from filed1062414e11.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:12:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat.hpp:83:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/csr_extract_u.hpp:6:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Sparse:26:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/SparseCore:66:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1062414e11.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:12:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat.hpp:83:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/csr_extract_u.hpp:6:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Sparse:27:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/OrderingMethods:71:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1062414e11.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:12:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat.hpp:83:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/csr_extract_u.hpp:6:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Sparse:29:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/SparseCholesky:43:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1062414e11.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:12:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat.hpp:83:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/csr_extract_u.hpp:6:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Sparse:32:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/SparseQR:35:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1062414e11.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:12:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat.hpp:83:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/csr_extract_u.hpp:6:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/Sparse:33:
#> In file included from /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/IterativeLinearSolvers:46:
#> /Users/jrnold/Library/R/3.4/library/RcppEigen/include/Eigen/src/Core/util/ReenableStupidWarnings.h:10:30: warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
#>     #pragma clang diagnostic pop
#>                              ^
#> In file included from filed1062414e11.cpp:638:
#> In file included from /Users/jrnold/Library/R/3.4/library/rstan/include/rstan/rstaninc.hpp:3:
#> In file included from /Users/jrnold/Library/R/3.4/library/rstan/include/rstan/stan_fit.hpp:36:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/services/optimize/bfgs.hpp:11:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/optimization/bfgs.hpp:9:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/optimization/lbfgs_update.hpp:6:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/circular_buffer.hpp:54:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/circular_buffer/details.hpp:20:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/move/move.hpp:30:
#> In file included from /Users/jrnold/Library/R/3.4/library/BH/include/boost/move/iterator.hpp:27:
#> /Users/jrnold/Library/R/3.4/library/BH/include/boost/move/detail/iterator_traits.hpp:29:1: warning: inline namespaces are a C++11 feature [-Wc++11-inline-namespace]
#> BOOST_MOVE_STD_NS_BEG
#> ^
#> /Users/jrnold/Library/R/3.4/library/BH/include/boost/move/detail/std_ns_begin.hpp:18:34: note: expanded from macro 'BOOST_MOVE_STD_NS_BEG'
#>    #define BOOST_MOVE_STD_NS_BEG _LIBCPP_BEGIN_NAMESPACE_STD
#>                                  ^
#> /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__config:439:52: note: expanded from macro '_LIBCPP_BEGIN_NAMESPACE_STD'
#> #define _LIBCPP_BEGIN_NAMESPACE_STD namespace std {inline namespace _LIBCPP_NAMESPACE {
#>                                                    ^
#> In file included from filed1062414e11.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:44:
#> /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/set_zero_all_adjoints.hpp:14:17: warning: unused function 'set_zero_all_adjoints' [-Wunused-function]
#>     static void set_zero_all_adjoints() {
#>                 ^
#> In file included from filed1062414e11.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core.hpp:45:
#> /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/core/set_zero_all_adjoints_nested.hpp:17:17: warning: 'static' function 'set_zero_all_adjoints_nested' declared in header file should be declared 'static inline' [-Wunneeded-internal-declaration]
#>     static void set_zero_all_adjoints_nested() {
#>                 ^
#> In file included from filed1062414e11.cpp:8:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/src/stan/model/model_header.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math.hpp:4:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/rev/mat.hpp:12:
#> In file included from /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat.hpp:58:
#> /Users/jrnold/Library/R/3.4/library/StanHeaders/include/stan/math/prim/mat/fun/autocorrelation.hpp:17:14: warning: function 'fft_next_good_size' is not needed and will not be emitted [-Wunneeded-internal-declaration]
#>       size_t fft_next_good_size(size_t N) {
#>              ^
#> 19 warnings generated.
```

```r
mod_pw
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
  // lag coefficient;
  real<lower = -1, upper = 1> theta;
  // convert range of theta from (0, 1) to (-1, 1)
  theta = (2. * theta_raw - 1.);
  // regression
  mu[1] = alpha + dot_product(beta, X[1, ]);
  mu[2:N] = alpha * (1 - theta) + (X[2:N, ] - theta * X[1:(N - 1), ]) * beta;
}
model {
  alpha ~ cauchy(alpha_loc, alpha_scale);
  beta ~ cauchy(beta_loc, beta_scale);
  theta_raw ~ beta(theta_a, theta_b);
  sigma ~ cauchy(0, sigma_scale);
  y[1] ~ normal(mu[1], sigma / sqrt(1 + theta ^ 2));
  y[2:N] ~ normal(mu[2:N], sigma);
}</code>
</pre>


```r
reagan_fit2 <- sampling(mod_pw, data = reagan_data)
```


```r
summary(reagan_fit2, par = c("alpha", "beta", "theta", "sigma"))$summary
#>           mean se_mean    sd   2.5%    25%    50%     75%  97.5% n_eff
#> alpha   48.366 0.16737 6.476 36.705 43.689 48.130 52.6559 61.579  1497
#> beta[1]  0.561 0.01305 0.682 -0.806  0.139  0.566  1.0225  1.848  2733
#> beta[2] -3.473 0.01755 0.780 -5.133 -3.970 -3.433 -2.9177 -2.066  1978
#> theta   -0.119 0.00407 0.152 -0.450 -0.218 -0.103 -0.0102  0.136  1394
#> sigma    6.938 0.01048 0.514  6.016  6.582  6.915  7.2597  8.043  2402
#>         Rhat
#> alpha      1
#> beta[1]    1
#> beta[2]    1
#> theta      1
#> sigma      1
```

[^reagan]: Example derived from Simon Jackman, "Reagan: linear regression with AR(1) disturbances," *BUGS Examples,* 2007-07-24, [URL](https://web-beta.archive.org/web/20070724034151/http://jackman.stanford.edu:80/mcmc/reagan.odc).

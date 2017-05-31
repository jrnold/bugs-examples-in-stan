
# Generalized Beetles: generalizing link functions for binomial GLMs


```r
library("rstan")
library("tidyverse")
```

GLMs rely on link functions, linking the linear predictors and the response probability, $\pi$. Logit and probit are perhaps the most familiar link functions,
mapping from the unit probability interval to the real line using the inverse CDFs of the logistic and standard Normal distributions, respectively. The logit and probit link functions have the interesting property that they are symmetric
about $\pi = 0.5$, and guarantee the effects of $x_i$ on $\pi$ to be greatest when $\pi  = 0.5$. To see this, recall that in GLMs for binomial data the effects of $x_i$ on $\pi$ are not constant, but vary over $\pi$.  For logit and probit, with link functions symmetric around zero, the effect of $x_i$  on $\pi$  is at its greatest when $f(x_i \beta)$ is its maximum, which for logit and probit occurs at $x_i \beta = 0$. In dose-response studies, this means that responsiveness to dose is at its greatest when subjects are on the cusp of a response, at, that is, when $E(\pi) = 0.5$.  In a study of voter turnout, ordinary logit or probit is estimated subject to the constraint that the effects of the covariates are at their greatest for citizens who are indifferent between turning out and abstaining [Nagler1994a]. Furthermore, for logit/probit, these marginal effects diminish in magnitude symmetrically as we move away from $E(\pi) = 0.5$; this symmetry follows from the symmetry of the logistic and normal PDFs/CDFs.

One can easily envisage situations where the researcher would not want to impose these features of the logit or probit link functions on their data. In many settings, knowledge of exactly where the marginal impact of the covariates is maximized is of tremendous practical importance, with implications for targeting policy interventions, resource allocation, and so on. For example, how to distribute resources for educational or health improvements? Given that the effects of interventions are not constant across a set of baseline probabilities,
knowing where proposed interventions are likely to have bigger or smaller effects is valuable information for policy-makers. As we have seen, logit/probit constrain these effects to be at their greatest at $E(\pi) = .5$, via their symmetric S-shaped link functions. Ceteris paribus we would prefer to estimate the shape of the link function from the data.

A relatively straightforward way to let the data be informative as to the shape of the link function is via a simple one-parameter transformation of the logit link [Prentice1976a]:
$$
\pi = \frac{1}{(1 + \exp(-x_i \beta))^m}
$$
where $m > 0$ is a parameter that skews the logit link. The standard logit model is a special case, where $m = 1$.

Estimating $m$ and $b$ by maximum likelihood is relatively straightforward, although there is little reason to believe the frequentist sampling distribution for $m$ is likely to be well approximated by the normal in a finite sample. Notice that $m$ enters the model in a highly non-linear fashion, and that different ranges of $m$ imply quite different relationships between the linear predictors and $\pi$.  In Bayesian terms, we can reasonably expect the posterior density of $m$ to be non-normal, and probably log-normal. Likewise, the posterior densities of related quantities of interest such as the LD50 could well be non-normal. Inferences for these quantities could well be misleading if we were to rely only on point estimates and asymptotic normal approximations; instead, a Bayesian approach via MCMC offers a way for us to obtain arbitrarily precise approximations to the posterior densities of these quantities.

To demonstrate the use of MCMC methods in this context, I use the famous beetles data of @Bliss1935a. 
These data have been extensively used by statisticians in studies generalized link functions [@Prentice1976a; @Stukel1988a], and are used by @SpiegelhalterBestGilks1996a to demonstrate how BUGS handles GLMs for binomial data. @CarlinLouis2000a use these data in an MCMC implementation of the one-parameter generalization used here; they
made use of re-parameterizations and problem-specific computer programming to handle the non-conjugacies introduced via the skew parameter $m$. Here I show how WinBUGS can handle this problem relatively simply.

$$
\begin{aligned}[t]
r_i &\sim \mathsf{Binomial}(n_i, \pi_i) \\
\pi_i &= 1 - \frac{1}{(1 + e^{(\alpha + x_i' \beta_i)})^\nu}
\end{aligned}
$$

I give $\mu$ a Gamma prior with a of 1.0 (corresponding to the standard logit model), and a standard deviation of 2,
$$
\nu \sim \mathsf{Gamma}(0.25, 0.25)
$$
The regression coefficients are given weakly informative priors,
$$
\begin{aligned}[t]
\alpha &\sim N(0, 10) \\
\beta &\sim N(0, 2.5)
\end{aligned}
$$


```r
genbeetles_mod <- stan_model("stan/genbeetles.stan")
```
<pre>
  <code class="stan">data {
  int N;
  int r[N];
  int n[N];
  vector[N] x;
}
parameters {
  real alpha;
  real beta;
  real<lower = 0.> nu;
}
transformed parameters {
  vector[N] mu;
  for (i in 1:N) {
    mu[i] = pow(inv_logit(alpha + beta * x[i]), nu) ;
  }
}
model {
  alpha ~ normal(0., 10.);
  beta ~ normal(0., 2.5);
  nu ~ gamma(0.25, 0.25);
  r ~ binomial(n, mu);
}
generated quantities {
  // probability where the maximum marginal effect
  real pdot;
  pdot = pow(inv_logit(nu), nu);
}</code>
</pre>


```r
genbeetles_data <- dget("data/genbeetles.R") %>%
  within({
    x <- as.numeric(scale(x))
  })
```



```r
genbeetles_fit <- sampling(genbeetles_mod, data = genbeetles_data)
```

```r
genbeetles_fit
#> Inference for Stan model: genbeetles.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
#> 
#>          mean se_mean   sd    2.5%     25%     50%     75%   97.5% n_eff
#> alpha   -1.34    0.03 0.90   -3.31   -1.91   -1.29   -0.73    0.24   809
#> beta     4.07    0.03 0.88    2.59    3.45    3.97    4.60    6.06   780
#> nu       0.34    0.00 0.14    0.17    0.25    0.32    0.40    0.69   775
#> mu[1]    0.10    0.00 0.03    0.06    0.08    0.10    0.12    0.16  1922
#> mu[2]    0.19    0.00 0.03    0.13    0.17    0.19    0.21    0.25  2856
#> mu[3]    0.33    0.00 0.03    0.27    0.31    0.33    0.35    0.39  4000
#> mu[4]    0.54    0.00 0.04    0.47    0.52    0.54    0.57    0.62  1713
#> mu[5]    0.77    0.00 0.03    0.70    0.75    0.77    0.80    0.83  2300
#> mu[6]    0.92    0.00 0.02    0.88    0.91    0.92    0.94    0.96  3217
#> mu[7]    0.98    0.00 0.01    0.95    0.97    0.98    0.99    0.99  1329
#> mu[8]    0.99    0.00 0.00    0.98    0.99    0.99    1.00    1.00  1058
#> pdot     0.84    0.00 0.04    0.76    0.81    0.84    0.87    0.90   815
#> lp__  -185.43    0.04 1.24 -188.64 -186.03 -185.11 -184.50 -183.97  1143
#>       Rhat
#> alpha    1
#> beta     1
#> nu       1
#> mu[1]    1
#> mu[2]    1
#> mu[3]    1
#> mu[4]    1
#> mu[5]    1
#> mu[6]    1
#> mu[7]    1
#> mu[8]    1
#> pdot     1
#> lp__     1
#> 
#> Samples were drawn using NUTS(diag_e) at Tue May 30 22:19:31 2017.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

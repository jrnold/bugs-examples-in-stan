
# Undervoting for President, by Race: Difference in Two Binomial Proportions {#undervote}


```r
library("tidyverse")
library("rstan")
```

Does undervoting for the US president differ by race?
Intentional undervoting is when a voter chooses not to cast vote for
an item on a ballot.

@TomzHouweling2003a analyze this phenomenon using two surveys:

-   Voter News Service (VNS) exit poll for the 1992 election
-   American National Election Studies (ANES) for the 1964--2000 elections

Each of these surveys asked voters whether they voted for president, as well as the race of the respondents.
The results of these surveys is contained in the `undervote` data frame.
The column `undervote` is the number of respondents who reported voting but not voting for president.


```r
undervote <- tribble(
  ~survey, ~race, ~n, ~undervote,
  "VNS", "black", 6537, 26,
  "VNS", "white", 44531, 91,
  "ANES", "black", 1101, 10,
  "ANES", "white", 9827, 57
  )
```


survey   race         n   undervote  Survey   Race     No. Voted   Didn't vote for president
-------  ------  ------  ----------  -------  ------  ----------  --------------------------
VNS      black     6537          26  VNS      black         6537                          26
VNS      white    44531          91  VNS      white        44531                          91
ANES     black     1101          10  ANES     black         1101                          10
ANES     white     9827          57  ANES     white         9827                          57

We are interested in analyzing the difference in proportions for each of these surveys independently.
We will model the proportions of each race and survey,
$$
\begin{aligned}[t]
y_i &\sim \mathsf{Binomial}(n_i, \pi_i) ,
\end{aligned}
$$
where
$$
i \in \{ (\text{VNS},\text{black}), (\text{VNS},\text{white}),  (\text{ANES},\text{black}),  (\text{ANES},\text{white}) \} .
$$

We will model the proportions independently by assigning them identical independent uninformative priors,

$$
\begin{aligned}[t]
\pi_i &\sim \mathsf{Beta}(1, 1) .
\end{aligned}
$$
The racial differences in undervoting in each survey are auxiliary quantities,
$$
\begin{aligned}[t]
\delta_{\text{VNS}} &= \pi_{\text{VNS},\text{black}} - \pi_{\text{VNS},\text{white}} ,\\
\delta_{\text{ANES}} &= \pi_{\text{ANES},\text{black}} - \pi_{\text{ANES},\text{white}} . \\
\end{aligned}
$$
We are also interested in the probability that black undervoting is greater than white undervoting, $\Pr(\delta_j) > 0$, in each survey.


```r
undervote_mod <- stan_model("stan/undervote.stan")
```
<pre>
  <code class="stan">data {
  int n[4];
  int y[4];
  vector[4] pi_a;
  vector[4] pi_b;
}
parameters {
  vector<lower = 0., upper = 1.>[4] pi;
}
model {
  y ~ binomial(n, pi);
  pi ~ beta(pi_a, pi_b);
}
generated quantities {
  vector[2] delta;
  int good[2];
  delta[1] = pi[2] - pi[1];
  delta[2] = pi[4] - pi[3];
  good[1] = int_step(delta[1]);
  good[2] = int_step(delta[2]);
}</code>
</pre>


```r
# this analysis depends on the order of the data frame
undervote_data <-
  list(y = undervote$undervote,
       n = undervote$n,
       N = nrow(undervote),
       pi_a = rep(1, 4),
       pi_b = rep(1, 4))
```


```r
undervote_fit <- sampling(undervote_mod, data = undervote_data)
```

```r
undervote_fit
#> Inference for Stan model: undervote.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
#> 
#>              mean se_mean   sd     2.5%      25%      50%      75%
#> pi[1]        0.00    0.00 0.00     0.00     0.00     0.00     0.00
#> pi[2]        0.00    0.00 0.00     0.00     0.00     0.00     0.00
#> pi[3]        0.01    0.00 0.00     0.01     0.01     0.01     0.01
#> pi[4]        0.01    0.00 0.00     0.00     0.01     0.01     0.01
#> delta[1]     0.00    0.00 0.00     0.00     0.00     0.00     0.00
#> delta[2]     0.00    0.00 0.00    -0.01    -0.01     0.00     0.00
#> good[1]      0.00    0.00 0.03     0.00     0.00     0.00     0.00
#> good[2]      0.08    0.00 0.26     0.00     0.00     0.00     0.00
#> lp__     -1254.94    0.03 1.38 -1258.39 -1255.60 -1254.63 -1253.94
#>             97.5% n_eff Rhat
#> pi[1]        0.01  4000    1
#> pi[2]        0.00  4000    1
#> pi[3]        0.02  4000    1
#> pi[4]        0.01  4000    1
#> delta[1]     0.00  4000    1
#> delta[2]     0.00  4000    1
#> good[1]      0.00  4000    1
#> good[2]      1.00  2850    1
#> lp__     -1253.23  2020    1
#> 
#> Samples were drawn using NUTS(diag_e) at Fri Apr 20 00:53:03 2018.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

## References {-}

Simon Jackman, "[Undervoting for President, by Race: difference in two binomial proportions](https://web-beta.archive.org/web/20070724034102/http://jackman.stanford.edu:80/mcmc/undervote.odc)", *BUGS Examples* 2007-07-24.

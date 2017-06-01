
# Undervoting for President, by Race: difference in two binomial proportions {#undervote}


```r
library("tidyverse")
library("rstan")
```


Does undervoting for the US president differ by race? 
Intentional undervoting is when a voter chooses not to cast vote for
an item on a ballot.

@TomzHouweling2003a analyze this phenomenon using two surveys:

- Voter News Service (VNS) exit poll for the 1992 election
- Americn National Election Studies (ANES) (1964-2000)

Each of these surveys asked voters whether they voted for president,
as well as the race of the respondents.

Of 6,537 black voters, 26 said they did not vote for president; of 44,531 white voters, 91 said they did not vote for president.
In the American National Election Studies (1964-2000), of 1,101 black voters, 10 report not voting for president, while 57 of 9,827 white voters report not voting for president.  


```r
undervote <- 
  tribble(
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
We'll model the proportions of each race and survey,
$$
\begin{aligned}[t]
y_i &\sim \mathsf{Binomial}(n_i, \pi_i) ,
\end{aligned}
$$
where
$$
i \in \{ (\text{VNS},\text{black}), (\text{VNS},\text{white}),  (\text{ANES},\text{black}),  (\text{ANES},\text{white}) \} .
$$
We'll model the proportions independently by assigning them identical independent uninformative priors,
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
#> good[1]      0.00    0.00 0.04     0.00     0.00     0.00     0.00
#> good[2]      0.07    0.00 0.26     0.00     0.00     0.00     0.00
#> lp__     -1255.02    0.03 1.47 -1258.73 -1255.72 -1254.66 -1253.94
#>             97.5% n_eff Rhat
#> pi[1]        0.01  4000    1
#> pi[2]        0.00  4000    1
#> pi[3]        0.02  4000    1
#> pi[4]        0.01  4000    1
#> delta[1]     0.00  4000    1
#> delta[2]     0.00  4000    1
#> good[1]      0.00  4000    1
#> good[2]      1.00  3143    1
#> lp__     -1253.24  1868    1
#> 
#> Samples were drawn using NUTS(diag_e) at Wed May 31 07:45:56 2017.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

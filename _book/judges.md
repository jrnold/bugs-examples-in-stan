
# Judges: estimating the ideological locations of Supreme Court justices {#judges}


```r
library("pscl")
library("tidyverse")
library("rstan")
```

This program implements an ideal-point model (similar to the legislators example), estimating both the locations of the justices on a latent ideological dimension, and two parameters specific to each case (corresponding to the item difficulty and item discrimination parameters of a two-parameter IRT model).[^judges-src]
The data consist of the decisions of Justices Rehnquist, Stevens, O'Connor, Scalia, Kennedy, Souter, Thomas, Ginsberg and Bryer, in that order, $i = 1, \dots , 9$.
The decisions are coded 1 for votes with the majority, and 0 for votes against the majority, and `NA` for abstentions.

In these models, the only observed data are votes, and the analyst wants to model those votes as a function of legislator- ($\theta_i$), and vote-specific ($\alpha_i$, $\lambda_i$) parameters.
The vote of legislator $i$ on roll-call $j$ ($y_{i,j}$) is a function of a the legislator's ideal point ($\theta_i$), the vote's difficulty parameter and the vote's discrimination ($\beta_j$):
$$
\begin{aligned}[t]
y_{i,j} &\sim \mathsf{Bernoulli}(\pi_i) \\
\pi_i &= \frac{1}{1 + \exp(-\mu_{i,j})} \\
\mu_{i,j} &= \beta_j \theta_i - \alpha_j
\end{aligned}
$$

$$
\begin{aligned}[t]
\beta_j &\sim \mathsf{Normal}(0, 2.5) \\
\alpha_j &\sim \mathsf{Normal}(0, 5) \\
\theta_i &\sim \mathsf{Normal}(0, 1) \\
\end{aligned}
$$



```r
data("sc9497", package = "pscl")
```
To simplify the analysis, the outcomes will be aggregated to "Yes", "No", and missing values (which 

```r
sc9497_vote_data <- tibble(vote = colnames(sc9497$votes)) %>%
  mutate(.vote_id = row_number())
                                 
sc9497_legis_data <- as.data.frame(sc9497$legis.names) %>%
  rownames_to_column("judge") %>%
  mutate(.judge_id = row_number())

sc9497_votes <- sc9497$votes %>%
  as.data.frame() %>%
  rownames_to_column("judge") %>%
  gather(vote, yea, -judge) %>%
  filter(!is.na(yea)) %>%
  inner_join(dplyr::select(sc9497_vote_data, vote, .vote_id), by = "vote") %>%
  inner_join(dplyr::select(sc9497_legis_data, judge, .judge_id), by = "judge")
```


```r
# mod_ideal_point <- stan_model("ideal_point.stan")
```

```r
# mod_ideal_point
```

[^judges-src]: This example is derived from Simon Jackman, "Judges: estimating the ideological locations of Supreme Court justices", *BUGS Examples*, 2007-07-24, [URL](https://web-beta.archive.org/web/20070724034049/http://jackman.stanford.edu:80/mcmc/judges.odc).



# Political Sophistication: item-response modeling with mixed data types {#sophistication}


```r
library("tidyverse")
library("rstan")
```

## Data


```r
data("PoliticalSophistication", package = "bayesjackman")
```

As part of a survey of French public opinion, 2,148 respondents were
administered a series of 19 items assessing their knowledge of political
leaders, political institutions, constitutional provisions, and the policies of
the political parties [@GrunbergMayerSniderman2002a].[^sophistication-src] Each
response is coded "correct" (1) or "incorrect" (0), and is modeled via a
two-parameter item-response model, with a logistic link function; each
respondent's level of political sophistication is the latent trait.

In addition, at the end of the thirty minute phone interview, the interviewer
assigned a score for each respondent's level of political information (based on
their impressions of the respondents formed over the course of the entire
interview) on a zero to twenty scale. These responses are modeled via a linear
regression, with each respondent's latent trait appearing as an unobserved
predictor, and an intercept specific to each interviewer (modeled
hierarchically in the code below). To uniquely orient the latent trait (higher
values corresponding to more political sophistication),  the interviewer
ratings are constrained to positively discriminate with respect to the latent
trait (see the constraint on the prior for gamma).

## Model

The survey data consists of 20 items, $y_1, \dots, y_20$. The first 19 items, $y_1, \dots, y_19$ binary responses to political information questions.
The final item, $y_20$, is a political sophistication score (0--20) assigned by the interviewer.

Let $y_{i,j}$ be the response of respondent $i \in 1, \dots, N$ to question $j \in 1, \dots, 20$.
$$
\begin{aligned}[t]
y_{i,j} &\sim \mathsf{Bernoulli}(\mathsf{Logit}^{-1}(\beta_j \xi_i - \alpha_j))
\end{aligned}
$$
for $j \in 1, \dots, j$.
Item 20 is modeled as
$$
\begin{aligned}[t]
y_{20,j} &\sim \mathsf{Normal}(\theta_i, \sigma^2) , \\
\theta_i &= \gamma \xi_i + \nu_{m[i]} .
\end{aligned}
$$
Since the question is assigned by the interviewer, $\theta_i$ is a linear function of a the latent score of the respondent ($\xi_i$) and an interviewer specific random effect, $\nu_{m[i]}$, where $m[i]$ means that $i$ was interviewed by interviewer $m \in 1, \dots, M$.
The interviewer effects are given a prior distribution,
$$
\nu_m \sim \mathsf{Normal}(0, \tau)  .
$$

To fix scale and location invariances, the respondents' abilities are given a standard normal distribution,
$$
\xi_m \sim \mathsf{Normal}(0, 1) .
$$
Since higher interviewer assessments should be associated with a higher latent political knowledge score, the rotation invariance is resolved by restricting the coefficient for the respondents to be positive,
$$
\gamma \sim \mathsf{HalfNormal}(0, 2.5 s_z) .
$$
The remaining parameters are assigned weakly informative priors.
$$
\begin{aligned}[t]
\tau &\sim \mathsf{HalfCauchy}(0, 5 s_{y_{20}}) , \\
\sigma &\sim \mathsf{HalfCauchy}(0, 5 s_{y_{20}})  , \\
\delta &\sim \mathsf{Normal}(10, 10 s_{y_{20}}) , \\
\beta_k &\sim \mathsf{Normal}(0, 2.5) , \\
\alpha_k &\sim \mathsf{Normal}(0, 10) .
\end{aligned}
$$
where $s_{y_{20}}$ is the scale for $y_20$.
We could use the empirical standard deviation of $y_20$, or an a-priori measure.
A value of $s_{y_{20}} = 21 / 4$, would place 95% of the mass of a normal distribution between 0 and 20.


```r
mod_sophistication <- stan_model("stan/sophistication.stan")
```

<pre>
  <code class="stan">data {
  // number of respondents
  int N;
  // number of items
  int K;
  // binary responses
  int<lower = 0, upper = 1> y_bern[K, N];
  // interviewer overall rating
  vector[N] y_norm;
  // interviewers
  int J;
  int<lower = 1, upper = J> interviewer[N];
  // priors
  real alpha_loc;
  real<lower = 0.> alpha_scale;
  real beta_loc;
  real<lower = 0.> beta_scale;
  real<lower = 0.> gamma_scale;
  real<lower = 0.> sigma_scale;
  real<lower = 0.> tau_scale;
  real delta_loc;
  real delta_scale;
}
parameters {
  // respondent latent score
  vector[N] xi_raw;
  // item discrimination
  vector[K] beta;
  // item difficulty
  vector[K] alpha;
  // coefficient in interviewer rating
  real<lower = 0.> gamma;
  // error in interviewer rating
  real<lower = 0.> sigma;
  // interviewer random effects
  vector[J] nu;
  // location of interviewer random effects
  real delta;
  // scale of interviewer random effects
  real<lower = 0.> tau;
}
transformed parameters {
  // interviewer rating
  vector[N] theta;
  // abilities
  vector[N] xi;
  xi = (xi_raw - mean(xi_raw));
  // respondent latent score
  for (i in 1:N) {
    theta[i] = gamma * xi[i] + nu[interviewer[i]];
  }
}
model {
  // priors
  xi_raw ~ normal(0., 1.);
  beta ~ normal(beta_loc, beta_scale);
  alpha ~ normal(alpha_loc, alpha_scale);
  gamma ~ normal(0., gamma_scale);
  sigma ~ cauchy(0., sigma_scale);
  tau ~ cauchy(0., tau_scale);
  delta ~ normal(delta_loc, delta_scale);
  // binary responses
  for (k in 1:K) {
    y_bern[k] ~ bernoulli_logit(beta[k] * xi - alpha[k]);
  }
  // interviewer random effects
  nu ~ normal(delta, tau);
  // interviewer score
  y_norm ~ normal(theta, sigma);
}</code>
</pre>

## Estimation


```r
y_scale <- 21 / 4

data_sophistication <- within(list(), {
  y_bern <- t(as.matrix(select(PoliticalSophistication, Q1:Q19)))
  N <- ncol(y_bern)
  K <- nrow(y_bern)
  y_norm <- PoliticalSophistication$Q20
  interviewer <- PoliticalSophistication$interviewer
  J <- max(interviewer)
  # priors
  sigma_scale <- 5 * y_scale
  xi_loc <- 0
  xi_scale <- 1
  alpha_loc <- 0
  alpha_scale <- 5
  beta_loc <- 0
  beta_scale <- 2.5
  gamma_scale <- 2.5 * y_scale
  # priors for interviewer effects
  tau_scale <- y_scale
  delta_loc <- 10
  delta_scale <- y_scale
})
```


```r
fit_sophistication <- sampling(mod_sophistication, data = data_sophistication, init = 0, chains = 1)
#> 
#> SAMPLING FOR MODEL 'sophistication' NOW (CHAIN 1).
#> 
#> Gradient evaluation took 0.003267 seconds
#> 1000 transitions using 10 leapfrog steps per transition would take 32.67 seconds.
#> Adjust your expectations accordingly!
#> 
#> 
#> Iteration:    1 / 2000 [  0%]  (Warmup)
#> Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Iteration: 2000 / 2000 [100%]  (Sampling)
#> 
#>  Elapsed Time: 238.869 seconds (Warm-up)
#>                97.5362 seconds (Sampling)
#>                336.405 seconds (Total)
```

## Questions / Extensions

1.  An alternative parameterization would place the political sophistication on the same 0-20 scale as `Q20`.
1.  Model `Q20` as an ordinal variable instead of a continuous variable.

[^sophistication-src]: Simon Jackman, "[Political Sophistication: item-response modeling with mixed data types](https://web-beta.archive.org/web/*/http://jackman.stanford.edu:80/mcmc/sophistication2002.odc)", *BUGS Examples*.

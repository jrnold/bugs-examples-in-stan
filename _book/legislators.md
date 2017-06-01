
# Legislators: estimating legislators' ideal points from voting histories (roll call data) {#legislators}


```r
library("tidyverse")
library("rstan")
```

Recorded votes in legislative settings (roll calls) are often used to recover the underlying preferences of legislators.  Political scientists analyze roll call data using the Euclidean spatial voting model: each legislator ($i = 1, \dots, n$) has a preferred policy position ($\xi$, a point in low-dimensional Euclidean space), and each vote ($j = 1, \dots, m$) amounts to a choice between an "Aye" ($\psi_y$) and a "Nay" location ($q_j$) respectively. Legislators are assumed to choose on the basis of utility maximization, with utilities (in one-dimension),
$$
\begin{aligned}[t]
U_i (q_j) &= -(\xi_i - q_j)^2 + \nu_{i,j} , \\ 
U_i (\psi_j) &= -(\xi_i - \psi_j)^2 + \omega_{i,j} ,
\end{aligned}
$$
where $\nu_{i,j}$ and $\omega_{i,j}$ are iid disturbances.
The utility differential $y_{i,j}^* = U_i(q_j) - U_i(\psi_j)$ is observed only in terms of sign.
Assuming of utility maximization, the legislator votes "Nay" when $y_{i,j}^* < 0$, and "Yea" when $y_{i,j}^* > 0$.
Simple algebra yields a hierarchical binary response model in which the legislators' ideal points appear as unobserved covariates: i.e.,
$$
y_{i,j}^* = \xi_i \beta_{j} - \alpha_j + \epsilon_{i,j}
$$
where the $\alpha_j$ and $\beta_j$ are functions of the vote-specific parameters $q_j$  and $\psi_j$.
Subject to identifying restrictions, the legislators' preferred positions can be estimated from the roll call data, via a logit or probit model.
In addition, with one-dimensional models, the "Yea" and "Nay" locations can be estimated
In higher dimensions, only the hyperplane separating the "Aye" and "Nay" locations can be recovered.

There is an extremely close correspondence between the statistical analysis of roll call data and item-response (IRT) models used in educational testing.  A two-parameter item-response model is equivalent to the statistical operationalization of the model described above, with the unobserved ideal point taking the part of the latent trait, and the item-discrimination parameters tapping ideological discrimination. @JohnsonAlbert1999a show that the two parameter IRT model is well suited to estimation and inference via Bayesian simulation.  The probability of an "Aye" vote is the probability that $y_{i,j}^* > 0$, if $\nu_{i,j}$ and $\omega_{i,j}$ have Type-1 extreme value distributions, then
$$
\begin{aligned}[t]
y_{i,j} &\sim \mathsf{Bernoulli}(\pi_{i,j}) \\
\mathsf{Logit}(\pi_{i,j}) &= \xi_i \beta_j - \alpha_j
\end{aligned}
$$

In the implementation below, an iid $\mathsf{Normal}(0,1)$ prior is used for the unobserved ideal points; vague normal priors are used for the $\beta$ and $\alpha$ parameters. For identification, two legislators are constrained to have fixed ideal points: here I set the ideal point of a legislator known to be liberal to -1.0, and the ideal point of a conservative to 1.0. Alternatively, we might employ restrictions or tight priors on particular $\beta$ parameters.

The following example uses data from the 106th U.S. Senate, which sat from January 1997 through October 1998, with n = 102 legislators.  Unanimous roll calls were dropped from the analysis. ...


- @Albert1992a
- @AlbertChib1993a
- @ClintonJackmanRivers2004a
- @EnelowHinich1984a
- @Jackman2000a
- @Jackman2001a
- @JohnsonAlbert1999a
- @Londregan2007a
- @PooleRosenthal2000a


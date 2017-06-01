
# Cosponsorship:  computing auxiliary quantities from MCMC output {#cosponsorship}

Typically, MCMC output consists of samples from the posterior density of model parameters.  But note that other quantities can be generated as well, say, imputations for missing data points, predictions, residuals, or goodness-of-fit summary statistics.  In fact, any function of the parameters can be calculated and output.

I demonstrate these ideas in the context of a generalized linear model for binary data.  The specific application is @Krehbiel1995a, a study of legislative behavior.  The dependent variable is a binary indicator, coded 1 if members of the U.S. House of Representatives chose to cosponsor bill HR 3266, and 0 otherwise (228 legislators cosponsored HR 3266, out of the 434 legislators for which data is available).  HR3266 was a wide-ranging spending bill designed to circumvent the usual budget-making process, that was considered by the 103rd House of Representatives in 1993/94.  Seven covariates are used in the analysis: a measure of each member's liberalism as measured by the interest group Americans for Democratic Action (ADA), a measure of fiscal conservatism published by the National Taxpayers' Union (NTU), an indicator variable for Democratic Party membership, a measure of Congressional seniority (years since first election, and thus inversely proportional to seniority), the electoral margin of the member, an indicator for membership of the House Appropriations Committee, and an indicator for membership of the House Budget Committee.  @Krehbiel1995a multivariate analysis finds that after controlling for legislators' policy preferences (as measured with the ADA and NTU scores), Democrats were actually more likely to support H.R. 3266 than Republicans.  Seniority is also a key predictor, with junior members more likely to cosponsor this legislation than members with greater seniority.

Several auxiliary quantities are estimated in the following BUGS program.  Percent-correctly-predicted (PCP) and a related quantity (expected PCP, or ePCP; see @Herron1999a) are easily coded using the BUGS "equals" and "step" constructs.   @Herron1999a also investigated the effect sizes of certain key predictors, by holding other covariates at fixed values and then comparing the predicted probabilities under various hypothetical scenarios.  Again, these are easily programmed and computed in WinBUGS; uncertainty in these quantities follows from the uncertainty in the parameters, and it is extremely easy to use WinBUGS to produce samples from the posterior density of these quantities.

I also compute latent residuals [@GelmanGoegebeurTuerlinckxEtAl2000a; @AlbertChib1995a].   These quantities are simply the difference between the latent y* and their fitted values.  In a typical logit or probit analysis, y* is completely ignored by the analyst.  But for Bayesian simulation, y* is central to the recovery of the parameters themselves: recall that Bayesian simulation in this context
largely amounts to treating the y* as missing data.  Hence not only are the y* recovered, but with very little additional effort we can also recover the latent residuals, at each iteration t.  Moreover, these latent residuals have uncertainty associated with them, stemming from uncertainty in b, and so we recover not just point estimates but their posterior densities.

$$
\begin{aligned}[t]
y_i &= \mathsf{Bernoulli}(\mu_i) \\
\mu_i &= \mathsf{Logit}^{-1}(x_i \beta)
\end{aligned}
$$

model{
	for (i in 1:N){	                   # loop over observations
		mu[i] <- beta[1] +
			beta[2]*ada[i] +
			beta[3]*ntu[i] +
			beta[4]*dem[i] +
			beta[5]*frosh[i] +
			beta[6]*margin[i] +
			beta[7]*approp[i] +
			beta[8]*budget[i];

		ystar[i] ~ dnorm(mu[i],1)I(lo[y[i]+1], up[y[i]+1]);  # trunc Normal

		## auxiliary quantities, individual specific
		e[i] <-  ystar[i] - mu[i];        # latent residuals
		probit(p[i]) <- mu[i];            # predicted probs

		## goodness of fit summaries
		llh[i] <- (equals(y[i],1)*log(p[i])) + (equals(y[i],0)*log(1-p[i]));
		epcp[i] <- (equals(y[i],1)*p[i]) + (equals(y[i],0)*(1-p[i]));
		pcp[i] <- equals(y[i],1)*step(p[i]-.5)
		        + equals(y[i],0)*step(.5-p[i]);
	}

	## aggregate auxiliary quantities
	sumllh <- sum(llh[]);     # sum likelihood contributions
	PCP <- sum(pcp[])/N;      # sum pcp contributions

	## clarify calculations
  ## probability of co-sponsorship for "average" member (median x vals)
	probit(pbar) <- beta[1] + beta[2]*0.55 + beta[3]*0.32
	              + beta[4] + beta[5]*86 + beta[6]*0.26;

	## party affiliation, effect size
	probit(p.dem) <- beta[1] + beta[2]*0.55 + beta[3]*0.32
	               + beta[4] + beta[5]*86 + beta[6]*0.26;
  probit(p.rep) <- beta[1] + beta[2]*0.55 + beta[3]*0.32
	               + beta[5]*86 + beta[6]*0.26;
	d.party <- p.dem - p.rep;

  ## "attributable effect", all sources, due to party affiliation
	probit(p.dem.all) <- beta[1] + beta[2]*0.8 + beta[3]*0.2
	                   + beta[4] + beta[5]*86 + beta[6]*0.28;
	probit(p.rep.all) <- beta[1] + beta[2]*0.1 + beta[3]*0.74
	                             + beta[5]*86 + beta[6]*0.24;
	d.party.all <- p.dem.all - p.rep.all

  ## bounds for truncated normal sampling
	lo[1] <- -50; lo[2] <-  0;
	up[1] <-   0; up[2] <- 50;

	## priors, multivariate normal, mean and precision in data file
	beta[1:k] ~ dmnorm(b0[ ] , B0[ , ]);
}

Initial Values (MLEs):
list(beta=c(-7.94, -.74, 8.82, 1.00, .052, -.051, -.24, -.858))

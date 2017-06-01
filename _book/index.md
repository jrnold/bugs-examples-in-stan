
--- 
title: "Bayesian Model Examples"
author: "Jeffrey B. Arnold and Simon Jackman"
date: "2017-05-31"
site: "bookdown::bookdown_site"
output:
  bookdown::gitbook: default
documentclass: book
bibliography: ["jackmanbayes.bib"]
biblio-style: apalike
link-citations: yes
---

# Preface {-}

This work contains some Bayesian model examples written by Simon Jackman and previously available on his website. 
These were originally written in WinBUGS or JAGS.
I have translated these examples into Stan and revised or edited them as appropriate.

The examples include:

1. [Undervote](undervote): difference of two independent proportions; racial differences in self-reported undervoting
2. [Cancer](cancer): difference of two independent proportions; differences in rates of lung cancer by smoking
3. [Florida](florida): learning about an unknown proportion from survey data; using survey data to update beliefs about support for Bush in Florida in the 2000 presidential election campaign
4. [Turnout](turnout2005): logit/probit models for binary response; voter turnout as a function of covariates
5. [Co-Sponsor](cosponsor): computing auxiliary quantities from MCMC output, such as residuals, goodness of fit; logit model of legislative co-sponsorship
6. [Reagan](reagan): linear regression with AR(1) disturbances; monthly presidential approval ratings for Ronald Reagan
7. [Political Sophistication](sophistication):  generalized latent variable modeling (item-response modeling with a mix of binary and ordinal responses); assessing levels of political knowledge among survey respondents in France
8. [Legislators](legislators):  generalized latent variable modeling (two-parameter item-response model); estimating legislative ideal points from roll call data
9. [Judges](judges): item response modeling; estimating ideological locations of Supreme Court justices via analysis of decisions
10. [Resistant](resistant): outlier-resistant regression via the t density; votes in U.S. Congressional elections, 1956-1994; incumbency advantage.
11. [House of Commons](uk92): analysis of compositional data; vote shares for candidates to the U.K. House of Commons
12. [Campaign](campaign): tracking a latent variable over time; support for candidates over the course of an election campaign, as revealed by polling from different survey houses.
13. [Aspirin](aspirin): meta-analysis via hierarchical modeling of treatment effects; combining numerous experimental studies of effect of aspirin on surviving myocardial infarction (heart attack)
14. [Corporatism](corporatism) hierarchical linear regression model, normal errors; joint impact of left-wing governments and strength of trade unions in structuring the determinants of economic growth
16. [Bimodal](bimodal): severe pattern of missingness in bivariate normal data; bimodal density over correlation coefficient
17. [Unidentified](unidentified): the consequences of over-parameterization; contrived example from Carlin and Louis
18. [Engines](engines): modeling truncated data; time to failure, engines being bench-tested at different operating temperatures
19. [Truncated](truncated): Example of sampling from a truncated normal distribution.
20. [Generalized Beetles](genbeetles): Generalizing link functions for binomial GLMs.
21. [Negative Binomial](negbin): Example of a negative binomial regression of homicides

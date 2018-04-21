
---
title: "Bayesian Model Examples"
author: "Jeffrey B. Arnold and Simon Jackman"
date: "2018-04-21"
site: "bookdown::bookdown_site"
output:
  bookdown::gitbook: default
documentclass: book
bibliography:
- "bayes.bib"
biblio-style: apalike
link-citations: yes
---

# Preface {-}

This work contains the Bayesian model examples written by Simon Jackman and previously available on his website.
These were originally written in WinBUGS or JAGS.
I have translated these examples into Stan and revised or edited them as appropriate.

This work is licensed under the [Creative Commons Attribution 4.0 International License](http://creativecommons.org/licenses/by/4.0/)

1.  [Undervote](undervote): difference of two independent proportions; racial differences in self-reported undervoting
1.  [Cancer](cancer): difference of two independent proportions; differences in rates of lung cancer by smoking
1.  [Florida](florida): learning about an unknown proportion from survey data; using survey data to update beliefs about support for Bush in Florida in the 2000 presidential election campaign
1.  [Turnout](turnout2005): logit/probit models for binary response; voter turnout as a function of covariates
1.  [Co-Sponsor](cosponsor): computing auxiliary quantities from MCMC output, such as residuals, goodness of fit; logit model of legislative co-sponsorship
1.  [Reagan](reagan): linear regression with AR(1) disturbances; monthly presidential approval ratings for Ronald Reagan
1.  [Political Sophistication](sophistication):  generalized latent variable modeling (item-response modeling with a mix of binary and ordinal responses); assessing levels of political knowledge among survey respondents in France
1.  [Legislators](legislators):  generalized latent variable modeling (two-parameter item-response model); estimating legislative ideal points from roll call data
1.  [Judges](judges): item response modeling; estimating ideological locations of Supreme Court justices via analysis of decisions
1.  [Resistant](resistant): outlier-resistant regression via the t density; votes in U.S. Congressional elections, 1956-1994; incumbency advantage.
1.  [House of Commons](uk92): analysis of compositional data; vote shares for candidates to the U.K. House of Commons
1.  [Campaign](campaign): tracking a latent variable over time; support for candidates over the course of an election campaign, as revealed by polling from different survey houses.
1.  [Aspirin](aspirin): meta-analysis via hierarchical modeling of treatment effects; combining numerous experimental studies of effect of aspirin on surviving myocardial infarction (heart attack)
1.  [Corporatism](corporatism) hierarchical linear regression model, normal errors; joint impact of left-wing governments and strength of trade unions in structuring the determinants of economic growth
1.  [Bimodal](bimodal): severe pattern of missingness in bivariate normal data; bimodal density over correlation coefficient
1.  [Unidentified](unidentified): the consequences of over-parameterization; contrived example from Carlin and Louis
1.  [Engines](engines): modeling truncated data; time to failure, engines being bench-tested at different operating temperatures
1.  [Truncated](truncated): Example of sampling from a truncated normal distribution.
1.  [Generalized Beetles](genbeetles): Generalizing link functions for binomial GLMs.
1.  [Negative Binomial](negbin): Example of a negative binomial regression of homicides

## Dependencies {-}

The R packages, Stan models, and datasets needed to run the code examples can be installed with

```r
# install.packages("devtools")
devtools::install_github("jrnold/jackman-bayes", subdir = "bayesjackman")
```

## Colonophon {-}


```r
sessionInfo()
#> R version 3.4.4 (2018-03-15)
#> Platform: x86_64-apple-darwin15.6.0 (64-bit)
#> Running under: macOS High Sierra 10.13.3
#> 
#> Matrix products: default
#> BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] methods   stats     graphics  grDevices utils     datasets  base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] Rcpp_0.12.16       knitr_1.20         magrittr_1.5      
#>  [4] munsell_0.4.3      colorspace_1.3-2   rlang_0.2.0       
#>  [7] stringr_1.3.0      plyr_1.8.4         tools_3.4.4       
#> [10] parallel_3.4.4     grid_3.4.4         gtable_0.2.0      
#> [13] xfun_0.1           htmltools_0.3.6    StanHeaders_2.17.2
#> [16] lazyeval_0.2.1     rprojroot_1.3-2    digest_0.6.15     
#> [19] tibble_1.4.2       rstan_2.17.3       bookdown_0.7.7    
#> [22] gridExtra_2.3      ggplot2_2.2.1      inline_0.3.14     
#> [25] evaluate_0.10.1    rmarkdown_1.9      stringi_1.1.7     
#> [28] pillar_1.2.1       compiler_3.4.4     scales_0.5.0      
#> [31] backports_1.1.2    stats4_3.4.4
```

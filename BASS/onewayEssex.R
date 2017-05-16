## one way anova
source("dataPrep.R")

## key var is math in data1
## data2 contains level2 covariates

## some soaking and poking
summary(data1$math)
table(data1$school)
length(unique(data1$school))

## classical one-way ANOVA
classical <- aov(math ~ as.factor(school),
                 data=data1)

## intra-class correlation
library(psychometric)
ICC1.lme(math,school,data1)
library(multilevel)
ICC1(classical)

library(survey)
d <- svydesign(ids=~school,data=data1)
clustered <- svymean(~math,
                     design=d,
                     deff="replace")

## lme4
## linear mixed model 
library(lme4)
l1 <- lmer(math ~ (1 | school),
           data=data1)

l1.omega <- as.numeric(summary(l1)@REmat[1,"Std.Dev."])
l1.sigma <- as.numeric(summary(l1)@REmat[2,"Std.Dev."])
l1.icc <- (l1.omega^2)/(l1.omega^2 + l1.sigma^2)

##########################################################
## MCMC via JAGS
library(rjags)
forJags <- list()
forJags$math <- data1$math
forJags$j <- match(data1$school,     ## indexes which school we're in
                   unique(data1$school))
forJags$N <- length(data1$math)      ## number of observations
forJags$J <- max(forJags$j)          ## number of schools

## initialize the Gibbs sampler
inits <- list(mu0=mean(forJags$math),
              alpha=tapply(forJags$math,
                forJags$j,mean),     ## school-specific means
              mu.y=forJags$math)

## compile JAGS
foo <- jags.model(file="oneway2.bug",
                  inits=inits,
                  data=forJags)

## get 4k iterations
out <- coda.samples(foo,4000,
                    variable.names=c("alpha","mu0",
                      "sigma.eps","sigma.mu"))

## raftery diagnostic
raftery.diag(out[,-c(1:160)])

## see summaries for all of the output
summary(out)

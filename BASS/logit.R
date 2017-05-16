## example 8.3 from BASS
library(pscl)
data(iraqVote)
summary(iraqVote)
attach(iraqVote)

## run a logistic regression
m1 <- glm(y ~ rep + gorevote,
          data=iraqVote,
          family=binomial)
summary(m1)

#######################
## subset to Dems only
## get data bundle for JAGS
forJags <- list(y=y[!rep],
                gorevote=gorevote[!rep]/100,
                n=sum(!rep),
                b0=rep(0,2),
                B0=diag(.001,2))

## initial values
inits <- list(list(beta=c(8,-15)))

## invoke JAGS
library(rjags)
foo <- jags.model(file="logitDemsOnly.bug",
                  inits=inits,
                  data=forJags)
out <- coda.samples(foo,
                    n.iter=25E3,
                    thin=1,
                    variable.names="beta")
                  

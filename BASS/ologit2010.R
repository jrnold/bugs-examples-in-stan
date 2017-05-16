#####################################
## read political information data
##
## get ready for JAGS
## run hierarchical model
##
#####################################
library(pscl)
data(politicalInformation)

forJags <- list(y=match(politicalInformation$y,
                  levels(politicalInformation$y)),
                x=cbind(politicalInformation$collegeDegree=="Yes",
                  politicalInformation$female=="Yes",
                  log(politicalInformation$age),
                  politicalInformation$homeOwn=="Yes",
                  politicalInformation$govt=="Yes",
                  log(politicalInformation$length)),
                id=match(politicalInformation$id,
                  unique(politicalInformation$id)))

## screen for missing (listwise deletion)
ok <- apply(forJags$x,1,
            function(x){
              all(!is.na(x))
            })
forJags$y <- forJags$y[ok]
forJags$x <- forJags$x[ok,]
forJags$N <- sum(ok)

## priors
forJags$b0 <- rep(0,6)
forJags$B0 <- diag(.01,6)

## initial values
inits <- list(list(beta=rep(0,6),
                   tau0=seq(-2,2,length=4))
              )
              
## JAGS
library(rjags)
foo <- jags.model(file="oLogit.bug",
                  inits=inits,
                  data=forJags)

out <- coda.samples(foo,
                    variable.names=c("beta","tau"),
                    n.iter=5e3)

library(pscl)
data(absentee)

denom <- absentee$absdem + absentee$absrep
y <- (absentee$absdem - absentee$absrep)/denom * 100
denom <- absentee$machdem + absentee$machrep
x <- (absentee$machdem - absentee$machrep)/denom *100

ols <- lm(y ~ x,
          subset=c(rep(TRUE,21),FALSE)  ## drop data point 22
          )

## set up for JAGS
forJags <- list(n=21,
                x=x[-22],
                y=y[-22],
                b0=rep(0,2),
                B0=solve(1000*diag(2)))

## poke holes in the data set
forJags$y[c(3,7,11)] <- NA
forJags$x[c(2,6,10)] <- NA

##forJags$NMISSX <- sum(is.na(forJags$x))
##forJags$indexMissingX <- which(is.na(forJags$x))

library(rjags)
foo <- jags.model(file="swissCheese.bug",
                  data=forJags)
nstore <- 5e3
thin <- 10
fun <- coda.samples(foo,
                    n.iter=nstore*thin,
                    thin=thin,
                    variable.names=c("y","x","beta"))


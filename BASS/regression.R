####################################################################
## a linear regression example
## section 6.3, pp256ff
library(pscl)
data(absentee)

attach(absentee)

## create variables for regression analysis
y <- (absdem - absrep)/(absdem + absrep)*100
x <- (machdem - machrep)/(machdem + machrep)*100

## exploratory graph, regression
plot(y~x,
		las=1,
		xlab="Democratic Margin, Machine Ballots",
		ylab="Democratic Margin, Absentee Ballots")
points(x[22],y[22],
			col="red",
			cex=3,pch=16)
abline(v=0,lty=2)
# ols linear regression, subset to 1st 21 observations
ols <- lm(y~x,
				subset=c(rep(TRUE,21),FALSE))
abline(ols,col="blue",lwd=3) ## overlay regression fit

## environment for passing to JAGS
forJags <- list(y=y[1:21],
                x=x[1:21],
                n=21,
                xstar=x[22])

## initial values, list containing one list, one chain
inits <- list(list(beta=c(0,0),
                   sigma=5,
                   ystar=0,
                   .RNG.seed=1234,   ## seed the RNG explicitly
                   .RNG.name="base::Mersenne-Twister"))

## p259 book
library(rjags)

## compilaton, no Gibbs sampler yet
foo <- jags.model(file="regression.bug",
                  n.chains=1,
                  data=forJags,
                  inits=inits)
update(foo,5E3)

# now we get Gibbs sampling
out <- coda.samples(model=foo,
                    variable.names=c("beta","sigma","ystar"),
                    n.iter=50E3)

summary(out)

## posterior predictive distribution for 
## disputed election
ystar <- out[[1]][,"ystar"]
hist(ystar,51,
     probability=TRUE,
     main="Posterior Predictive Density,\nDisputed Absentee Ballot Result")
abline(v=y[22],col="red",lwd=3)


## Bayesian p-value
table(ystar>y[22])
mean(ystar>y[22])

##################################
## difference of two binomial
## proportions
##
## simon jackman
## august 12 2010
#################################

## Example 3.2 from BASS
nsims <- 1e6
theta1 <- rbeta(nsims,2,8)
theta0 <- rbeta(nsims,3,75)
q <- theta1 - theta0
summary(q)
mean(q>0)

## now done in JAGS
data <- data.frame(r=c(1,2),
                   n=c(8,76))

## hand-off to JAGS
library(rjags)

## make JAGS compile the model with data
m <- jags.model(data=data,
                file="jas.bug")

out <- coda.samples(m,             ## the Bayesian model
                    n.iter=500e3,  ## desired number of samples
                    variable.names=c("theta","q")) ## what to store


summary(out)
plot(out)





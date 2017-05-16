## two-way ANOVA via JAGS
##
## simon jackman
## 11/17/2010
#######################################
library(pscl)

data <- presidentialElections
data <- data[data$state != "DC" & data$year > 1971,]
data$j <- match(data$state,unique(data$state))
data$k <- match(data$year,sort(unique(data$year)))

## a data bundle that will go to JAGS
forJags <- list(y=data$demVote,
                j=data$j,
                N=dim(data)[1],
                J=max(data$j))

## initial values for 2 chains
inits <- list(list(mu=20,omega=1,sigma=2,
                   .RNG.seed=1234,
                   .RNG.name="base::Mersenne-Twister"),
              list(mu=80,omega=9,sigma=9,
                   .RNG.seed=99999,
                   .RNG.name="base::Mersenne-Twister"))

## hand off to JAGS
library(rjags)
## model compilation phase
foo <- jags.model(file="madridWednesday.bug",
                  n.chains=2,
                  inits=inits,
                  data=forJags)

## Gibbs sampler is actually here
out <- coda.samples(foo,
                    n.iter=10E3,
                    thin=1,
                    variable.names=c("mu","sigma","omega"))

## convergence and run-length diagnostics
sink("geweke.out")
geweke.diag(out)
sink()
sink("raftery.out")
raftery.diag(out)
sink()

## extract alphas in separate object
out <- coda.samples(foo,
                    n.iter=5E3,
                    thin=1,
                    variable.names=c("alpha"))

## stack the output from the two chains on top of one another
alpha <- rbind(out[[1]],
               out[[2]])
dimnames(alpha) <- list(NULL,unique(data$state))

alphaBar <- colMeans(alpha)   ## column means of matrix alpha
names(alphaBar) <- unique(data$state) 

indx <- order(alphaBar)   ## permutation of indices that order alphaBar

alphaCI <- apply(alpha,2,function(x)HPDinterval(as.mcmc(x)))
alphaCI <- t(alphaCI)

## compute posterior probability that Utah < Wyoming
mean(alpha[,"Utah"] < alpha[,"Wyoming"])

## "caterpillar" plot of alphas (a la league table)
par(mar=c(3,5,3,1))

plot(x=alphaBar[indx],
     xlim=range(alphaCI),
     y=1:50,
     pch=16,
     xlab=expression(alpha),
     ylab="",
     axes=FALSE)

for(i in 1:50){
  lines(y=c(i,i),
        x=alphaCI[indx[i],])
}

text(x=alphaCI[indx,1],
     y=1:50,
     pos=2,
     unique(data$state)[indx])


###############################################
## which state is the median (occupies rank 26)?
r <- apply(alpha,1,rank)
r <- t(r)
medianPostProb <- apply(r,2,function(x)mean(x==26))

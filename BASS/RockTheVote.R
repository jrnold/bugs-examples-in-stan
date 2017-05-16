## RockTheVote example
library(pscl)
data(RockTheVote)
attach(RockTheVote)

## look at the data
length(r)                  ## 85 cable systems
length(unique(strata))     ## 40 of these, obviously unbalanced

sum(treated)               ## 42 treated units
tapply(treated,strata,sum) ## number of treated units per strata

## MLE of average treatment effect
## "complete pooling"
m1 <- glm(cbind(r,n-r) ~ treated,
          data=RockTheVote,
          family=binomial)

## add fixed effects for strata
m2 <- update(m1, ~ . + as.factor(strata))
summary(m1)$coef[2,]
summary(m2)$coef[2,]

## no pooling, analysis in each strata
## estimate MLEs of treatment effects

## this function runs a logit in each strata 
deltaFunction <- function(data){
  ## logit model
  model <- glm(cbind(r,n-r)~treated,
               data=data,
               family=binomial)
  ## estimate of the treatment effect
  ## and 95% CI
  c(coef(model)[2],
    confint(model)[2,])
}

tmp <- by(RockTheVote,
          as.factor(RockTheVote$strata),
          deltaFunction,
          simplify=TRUE)

tmp <- matrix(unlist(tmp),ncol=3,byrow=TRUE)

indx <- order(tmp[,1])

plot(y=1:40,
     x=tmp[indx,1],
     pch=16,cex=1.25,
     xlim=range(tmp),
     ylab="",
     axes=FALSE,
     xlab="Estimated Treatment Effect (MLEs, Logit Scale)")
text(y=1:40,
     x=par()$usr[1],
     as.character((1:40)[indx]),
     cex=.5)
segments(x0=tmp[indx,2],
         x1=tmp[indx,3],
         y0=1:40,
         y1=1:40)
axis(1)
axis(3)
abline(v=0)
abline(v=coef(m1)[2],lty=2)  ## overlay avg effect

#####################################
## hand to JAGS
library(rjags)
inits <- list(list(alpha=rep(0,40),
                   delta=rep(0,42)))
                   
foo <- jags.model(file="binomial.bug",
                  data=RockTheVote,
                  inits=inits)
update(foo,10e3)
out <- coda.samples(foo,
                    variable.names=c("alpha","delta",
                      "delta.new","p.new",
                      "mu","sigma","pvalue"),
                    n.iter=50000,
                    thin=5)

## now extract the deltas (42 of them!)
deltaCols <- grep("delta\\[",dimnames(out[[1]])[[2]])
delta <- out[[1]][,deltaCols]

## pointwise 95 CIs
deltaSummary <- apply(delta,2,
                      function(x){
                        c(mean(x),
                          quantile(x,c(.025,.975)))
                      })
deltaSummary <- t(deltaSummary)

## now plot means and CIs
indx <- order(deltaSummary[,1])       ## permutation that sorts the means
plot(y=1:42,
     x=deltaSummary[indx,1],          ## plot means
     xlim=range(deltaSummary[,2:3]),  ## make graph wide enough
     pch=16,                          ## solid plotting symobl
     axes=FALSE,                      ## no axes for now
     xlab=expression(delta),          ## Greek letter delta in axis label
     ylab="")                         ## no label on vertical axis

## 95% CIs
segments(x0=deltaSummary[indx,2],     ## x locations
         x1=deltaSummary[indx,3],
         y0=1:42,                     ## y locations
         y1=1:42)
axis(1)
axis(3)
abline(v=0)  ## vertical reference line at zero

## labels
text(x=rep(par()$usr[1] - .05,42),
     y=1:42,
     xpd=NA,       ## allow printing out of bounds
     (1:42)[indx], ## the labels (just sorted ids)
     pos=4,        ## print to the right
     cex=.5)       ## smaller than usual size

## transform to ranks
deltaRank <- apply(delta,1,rank)
deltaRank <- t(deltaRank)

## look at posterior probability of occupying rank 1
probTheBest <- apply(deltaRank,2,function(x)mean(x==1))
indx <- order(probTheBest)
plot(probTheBest[indx],
     1:42,
     xlab="Posterior Probability of Occupying Rank 1",
     ylab="",
     axes=FALSE)
text(x=rep(par()$usr[1] - .05,42),
     y=1:42,
     xpd=NA,
     (1:42)[indx],
     pos=4,
     cex=.5)
axis(1)



     

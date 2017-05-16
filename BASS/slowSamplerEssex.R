library(pscl)
data(absentee)

denom <- absentee$absdem + absentee$absrep
y <- (absentee$absdem - absentee$absrep)/denom * 100
denom <- absentee$machdem + absentee$machrep
x <- (absentee$machdem - absentee$machrep)/denom *100

ols <- lm(y ~ x,
          subset=c(rep(TRUE,21),FALSE)  ## drop data point 22
          )

## set up for JAGS, "slow"/inefficient sampler
forJags <- list(n=21,
                x=x[-22],
                y=y[-22],
                b0=rep(0,2),
                B0=solve(1000*diag(2)))

library(rjags)
foo <- jags.model(file="separateSampler.bug",
                  data=forJags)
separate <- coda.samples(foo,n.iter=1000,
                         variable.names="beta")

## now run the "blocked" parameterization
forJags$b0 <- rep(0,2)
forJags$B0 <- solve(1000*diag(2))

foo <- jags.model(file="jointSampler2.bug",
                  data=forJags)
joint <- coda.samples(foo,n.iter=1000,
                         variable.names="beta")


pdf(file="compareSamplers.pdf",
    height=8.5,
    width=6.5)
par(mfrow=c(5,2),
    las=1,
    cex.main=.85,
    cex.axis=.85,
    mgp=c(1.25,.5,0),
    mar=c(2.5,2.25,1,.25))

########################################################
## traceplots
traceplot(as.mcmc(separate[,1]),bty="n")
text(x=mean(par()$usr[1:2]),
     y=par()$usr[4],
     xpd=TRUE,
     expression(paste("Unblocked, Trace: ",beta[1])))

traceplot(as.mcmc(separate[,2]),bty="n")
text(x=mean(par()$usr[1:2]),
     y=par()$usr[4],
     xpd=TRUE,
     expression(paste("Unblocked, Trace: ",beta[2])))

traceplot(as.mcmc(joint[,1]),bty="n")
text(x=mean(par()$usr[1:2]),
     y=par()$usr[4],
     xpd=TRUE,
     expression(paste("Blocked, Trace: ",beta[1])))

traceplot(as.mcmc(joint[,2]),bty="n")
text(x=mean(par()$usr[1:2]),
     y=par()$usr[4],
     xpd=TRUE,
     expression(paste("Blocked, Trace: ",beta[2])))

########################################################
## acfs
plot(acf(separate[[1]][,1],plot=FALSE),
     main="",
     bty="n",
     xlab="Lag",
     ylab="")
text(x=mean(par()$usr[1:2]),
     y=par()$usr[4],
     pos=1,
     expression(paste("Unblocked ACF: ",beta[1])))

plot(acf(separate[[1]][,2],plot=FALSE),
     main="",
     bty="n",
     xlab="Lag",
     ylab="")
text(x=mean(par()$usr[1:2]),
     y=par()$usr[4],
     pos=1,
     expression(paste("Unblocked ACF: ",beta[2])))

plot(acf(joint[[1]][,1],plot=FALSE),
     main="",
     bty="n",
     xlab="Lag",
     ylab="")
text(x=mean(par()$usr[1:2]),
     y=par()$usr[4],
     pos=1,
     expression(paste("Blocked ACF: ",beta[1])))

plot(acf(joint[[1]][,2],plot=FALSE),
     main="",
     bty="n",
     xlab="Lag",
     ylab="")
text(x=mean(par()$usr[1:2]),
     y=par()$usr[4],
     pos=1,
     expression(paste("Blocked ACF: ",beta[2])))

## densities
d1 <- density(unclass(joint[[1]])[,1])
d2 <- density(unclass(separate[[1]])[,1])
xlim <- range(d1$x,d2$x)
ylim <- range(d1$y,d2$y)

plot(d1,
     xlim=xlim,
     ylim=ylim,
     lwd=2,
     xlab=expression(beta[1]),
     ylab="",
     axes=FALSE,
     main="")
axis(1)
lines(d2$x,d2$y,col=gray(.45),lwd=2)
legend(x="topleft",
       bty="n",
       legend=c("Blocked","Unblocked"),
       lwd=c(2,2),
       col=c("black",gray(.45)))

d1 <- density(unclass(joint[[1]])[,2])
d2 <- density(unclass(separate[[1]])[,2])
xlim <- range(d1$x,d2$x)
ylim <- range(d1$y,d2$y)

plot(d1,
     xlim=xlim,
     ylim=ylim,
     lwd=2,
     xlab=expression(beta[2]),
     ylab="",
     axes=FALSE,
     main="")
axis(1)
lines(d2$x,d2$y,col=gray(.45),lwd=2)
legend(x="topleft",
       bty="n",
       legend=c("Blocked","Unblocked"),
       lwd=c(2,2),
       col=c("black",gray(.45)))


dev.off()
###################################################################
## long run for raftery-lewis
foo <- jags.model(file="separateSampler.bug",
                  data=forJags)
separate <- coda.samples(foo,n.iter=15000,
                         variable.names="beta")

foo <- jags.model(file="jointSampler2.bug",
                  data=forJags)
joint <- coda.samples(foo,n.iter=15000,
                         variable.names="beta")

###################################################################

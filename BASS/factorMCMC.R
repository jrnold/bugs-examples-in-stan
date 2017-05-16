## MCMC for a simple one factor model
## this uses data that is shipped as part of R, state.x77 data frame
y <- state.x77[,c("Income","Illiteracy","Life Exp","Murder","HS Grad")]
y[,1] <- y[,1]/1000

ystar <- sweep(y,2,apply(y,2,mean))

## conventional factor analysis
f1 <- factanal(y,1,scores="regression")
lambda <- unclass(f1$loadings)
psi <- diag(f1$uniquenesses)
sigma <- var(y)

gamma <- t(lambda)%*%solve(psi)%*%lambda
beta <- t(solve(1 + gamma)%*%t(lambda)%*%solve(psi))
x1 <- scale(y,TRUE,TRUE)%*%beta
x0 <- f1$scores

## look at the correlation matrix
print(cor(y))

## make a data frame for passing to JAGS
forJags <- list(y=y,n=dim(y)[1],
                g0=rep(0,2),
                G0=diag(.01,2))


######################################################
## initial values
g <- matrix(NA,5,2)
omega <- rep(NA,5)
## run regressions to get start values
for(j in 1:5){
  tmpReg <- lm(y[,j] ~ x0)
  g[j,] <- coef(tmpReg)
  omega[j] <- summary(tmpReg)$sigma
}

inits <- list(list(gamma=g,
                   tau=1/omega^2,
                   xistar=as.vector(x0),
                   .RNG.name="base::Mersenne-Twister",
                   .RNG.seed=314159))
#######################################################
## call rjags
#######################################################
require(rjags)
foo <- jags.model(file="factorMCMC.bug",
                  data=forJags,
                  inits=inits)
out <- coda.samples(foo,
                    variable.names=c("xi","gamma","omega"),
                    n.iter=5e3,
                    thin=1)

extractFromCoda <- function(coda,vname){
  dnames <- dimnames(coda[[1]])[[2]]
  theString <- paste("^",vname,sep="")
  cat(paste("searching for",theString,"in MCMC output\n"))
  theOnes <- grep(theString,
                  dnames)
  cat(paste("found",length(theOnes),"matching columns\n"))
  as.matrix(coda[[1]][,theOnes])
}

x <- extractFromCoda(out,"x")
g <- extractFromCoda(out,"gamma")
omega <- extractFromCoda(out,"omega")

nsims <- dim(x)[1]

## summaries for x
myfunc <- function(x)c(mean(x),
                       quantile(x,c(.025,.975)))
xbar <- apply(x,2,myfunc)
xsd <- apply(x,2,sd)
xbar <- t(xbar)
dimnames(xbar) <- list(dimnames(state.x77)[[1]],
                       c("Mean","2.5%","97.5%"))

gbar <- t(apply(g,2,myfunc))  ## summaries for gammas
omegabar <- t(apply(omega,2,myfunc))

## r-squared
r2 <- matrix(NA,nsims,5)
yss <- apply(y,2,function(x)sum((x-mean(x))^2))
for(iter in 1:nsims){
  rho <- cor(cbind(x[iter,],y))[2:6,1]
  r2[iter,] <- rho^2
}  
r2bar <- t(apply(r2,2,myfunc))

#########################################################################
## graphs 
indx <- order(xbar[,1])          ## sort ideal point estimates
cols <- rep("black",50)          ## colors
southernStates <- c("Louisiana","Mississippi","South Carolina",
"Alabama","Georgia","North Carolina","Kentucky",
                    "Tennessee","Texas","Arkansas","Virginia","Florida")
cols[dimnames(xbar)[[1]] %in% southernStates] <- gray(.45)

pdf(file="xLatentVariable.pdf",
    width=11.5,height=8)
par(mfrow=c(1,2),                 ## plotting options
    mar=c(2,0.1,2,0.1),
    omi=rep(0,4))
plotlims <- range(xbar)
axislims <- c(-2.5,2)
axisticks <- seq(from=axislims[1],to=axislims[2],by=.5)
axislabs <- as.character(axisticks)
textloc <- 1.01*xbar[,2]
textloc[textloc > -3] <- 1.01 * -3

plot(x=plotlims,y=c(1,25),       ## empty plotting region
     type="n",
     cex=.5,  
     ylim=c(.5,25.5),yaxs="i",
     axes=FALSE,
     xlab="",ylab="")

for(i in 1:25){    
  lines(y=c(i,i),  ## line for confidence interval
        x=c(xbar[indx[i],2],xbar[indx[i],3]),
        lwd=2.5)
  points(y=i,x=xbar[indx[i],1],  ## overlay posterior mean
         col=cols[indx[i]],
         pch=19,cex=1.25)
  text(xbar[indx[i],3],i,        ## text, state name
       dimnames(xbar)[[1]][indx[i]],
       adj=-0.1)
}
axis(1,at=axisticks,labels=axislabs,cex=.5,lwd=1.5)
axis(3,at=axisticks,labels=axislabs,cex=.5,lwd=1.5)

plot(x=plotlims,y=c(1,25),cex=.5,
     type="n",axes=F,xlab="",ylab="",
     ylim=c(0.5,25.5),yaxs="i")
for(j in 1:25){
  i <- 25 + j
  lines(y=c(j,j),
        x=c(xbar[indx[i],2],xbar[indx[i],3]),lwd=2.5)
  points(y=j,x=xbar[indx[i],1],
         col=cols[indx[i]],
         pch=19,cex=1.25)
  text(xbar[indx[i],2],j,
       dimnames(xbar)[[1]][indx[i]],
       adj=1.1)
}
axis(1,at=axisticks,labels=axislabs,cex=.5,lwd=1.5)
axis(3,at=axisticks,labels=axislabs,cex=.5,lwd=1.5)
dev.off()

##########################################################################
## compare with factor analysis
pdf("compareFactorScore.pdf")
f1 <- factanal(y,1,scores="regression")
plot(f1$scores,xbar[,1],
     xlab="Factor Scores (Regression Scoring)",
     ylab="Posterior Mean of Latent Variable",
     pch=19,
     col=cols)
abline(0,1)
legend(x="topleft",
       bty="n",
       legend=c("South","Non-South"),
       col=c("red","black"),
       pch=19)
dev.off()
#########################################################################

## trace plots
pdf(file="factorAnalysisTracePlots.pdf",
    width=11.5,height=8)
par(mfrow=c(4,1),                 ## plotting options
    omi=rep(0,4),
    mar=c(3,4,.5,1),
    mgp=c(1.5,.5,0),
    cex.axis=.65)

dimnames(x)[[2]] <- dimnames(xbar)[[1]]
plot.ts(x[,"California"],
        xlab="Iterations",
        ylab="California",
        xaxs="i")
abline(h=xbar["California",],col=gray(.45),lwd=2)
lf <- loess(x[,"California"] ~ I(1:10000),span=.10)
lines(1:10000,predict(lf),col="blue",lwd=2)

plot.ts(x[,"Louisiana"],
        xlab="Iterations",
        ylab="Louisiana",
        xaxs="i")
abline(h=xbar["Louisiana",],col=gray(.45),lwd=2)
lf <- loess(x[,"Louisiana"] ~ I(1:10000),span=.10)
lines(1:10000,predict(lf),col="blue",lwd=2)

dimnames(g)[[1]] <- NULL
plot(1:10000,
     g[,"gamma[1,2]"],
     xlab="Iterations",
     ylab="Slope,Income Indicator",
     xaxs="i",type="l")
abline(h=gbar["gamma[1,2]",],col=gray(.45),lwd=2)
lf <- loess(g[,"gamma[1,2]"] ~ I(1:10000),span=.10)
lines(1:10000,predict(lf),col="blue",lwd=2)

plot(1:10000,
     omega[,1],
     xlab="Iterations",
     ylab="Measurement Error StdDev,\nIncome Indicator",
     xaxs="i",type="l")
abline(h=omegabar[1,],col=gray(.45),lwd=2)
lf <- loess(omega[,1] ~ I(1:10000),span=.10)
lines(1:10000,predict(lf),col="blue",lwd=2)

dev.off()

################################################################
plotFunc <- function(x,lab){
  hist(x,51,probability=TRUE,
       main="",
       xlab=as.expression(parse(text=lab)),
       ylab="",
       axes=FALSE)
  axis(1)
  if(lab=="rho")
    abline(h=.5,lty=3)
  ##if(lab=="omega[12]")
  ##  lines(density(s12),lty=3)
  ##if(lab=="omega[22]")
  ##  lines(density(s2),lty=3)

  plot.ts(x,
          bty="n",
          axes=FALSE,
          ylab=as.expression(parse(text=lab)),
          xlab="Iteration")
  axis(1)
  axis(2)
          
  
  plot(acf(x,plot=FALSE),
       ci=0,
       axes=FALSE,
       xaxs="i",
       las=1,
       ylab="",
       main="")
  
  axis(1)
  axis(2,las=1)
  labExpr <- as.expression(parse(text=paste("Autocorrelations: ",lab)))
  mtext(side=3,
        line=-1,
        cex=.75,
        labExpr)
  text(x=par()$usr[2],
       y=.20,
       adj=c(1,-1),
       xpd=TRUE,
       cex=.85,
       paste("Effective sample size:",
             prettyNum(round(effectiveSize(x)),",")
             )
       )
}

pdf(file="threeColumn.pdf",
    w=6.5,
    h=7)
par(mfrow=c(4,3),
    cex.axis=.65,
    cex.lab=.85,
    cex.main=.85,
    mgp=c(1.1,.3,.0),
    las=1,
    tcl=-.2,
    mar=c(2.5,2,.5,1))

plotFunc(out[[1]][,"gamma[1,2]"],
         lab="gamma[12]")
plotFunc(out[[1]][,"omega[4]"],
         lab="omega[4]")
plotFunc(out[[1]][,"xi[5]"],
         lab="xi[5]")
plotFunc(out[[1]][,"xi[18]"],
         lab="xi[18]")
dev.off()
#################################################################

#################################################################
## acf plot
rho <- apply(unclass(out[[1]]),2,function(x)acf(x,plot=FALSE)$acf[2,1,1])
rho <- rho[order(rho)]
pdf(file="rho.pdf",
    h=6.5,
    w=4)
par(mar=c(2.5,3.5,2.5,.5),
    mgp=c(1.1,.3,0),
    tcl=-.2,
    cex.axis=.65)
plot(x=c(0,max(rho)),
     y=c(.5,length(rho)+.5),
     type="n",
     xaxs="i",
     yaxs="i",
     xlab=expression(rho),
     ylab="",
     axes=FALSE)
abline(h=1:length(rho),
       lty=3,
       col=gray(.75))

pch <- rep(16,length(rho))
pch[grep("omega",names(rho))] <- 1
pch[grep("gamma",names(rho))] <- 22

for(i in 1:length(rho)){

  points(x=rho[i],
         y=i,
         pch=pch[i],
         xpd=TRUE)

  theText <- names(rho)[i]
  theText <- sub(theText,pattern="x",replacement="xi")
  theText <- sub(theText,pattern=",",replacement="")
  text(x=rho[i]-.03,
       y=i,
       xpd=TRUE,
       as.expression(parse(text=theText)),
       cex=.5)
}
axis(1)
axis(3)
dev.off()

#################################################################


## ranks
xrank <- t(apply(x,1,rank))
## a function to compute some summary statistics
summaryFunc <- function(x,q=c(.025,.975)){
  c(mean(x),
    quantile(x,q))
}
## now use the summary function over the columns of xrank
xrank.Summary <- t(apply(xrank,2,summaryFunc))

xpost <- apply(xrank,2,function(x)mean(x==1))
round(rev(sort(xpost))[1:10],2)

pdf(file="xrankLowest.pdf",
    w=8.5,h=8.5)
par(mfrow=c(3,3),                 ## 3 by 3 grid of graphs
    mar=c(3.35,4,3.5,1),          ## axis margins
    las=1,                
    mgp=c(2.5,.75,0),            ## tighten up axis labeling
    omi=rep(0,4),                 ## no outer margins
    xpd=TRUE)
indx <- rev(order(xpost))         ## sort by prob of being most liberal
them <- indx[1:9]                 ## 12 highest probs
for(i in 1:9){
  hist(xrank[,them[i]],
       breaks=0:max(xrank[,them]) + .5,
       probability=TRUE,
       col=gray(.45),
       xlab="Rank",
       ylab="Probability",
       ylim=c(0,.30),
       axes=FALSE,
       main=dimnames(xbar)[[1]][them[i]])
  text(x=3,
       y=.29,
       paste("Probability of occupying\n",
             "lowest rank: ",
             format(xpost[them[i]],
                    nsmall=2,digits=2,
                    scientific=FALSE,
                    trim=TRUE),
             sep=""),
       adj=0)
  offset <- c(.5,.0075)
  ##   if(xpost[them[i]]<.51)
  ##     arrows(x0=1,
  ##            y0=.51,
  ##            x1=1,
  ##            y1=xpost[them[i]]+offset[2],
  ##            length=.08,lwd=2,
  ##            angle=25,
  ##            col="black")
  axis(2)
  axis(1,at=c(1,seq(5,max(xrank[,them]),by=5)))
  rug(xrank.Summary[them[i],2],lwd=3)
  rug(xrank.Summary[them[i],3],lwd=3)
}
dev.off()


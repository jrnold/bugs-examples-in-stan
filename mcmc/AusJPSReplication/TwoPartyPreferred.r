#################################################################
## assume data has been read in, see read.r
##
## dumps data and runs JAGS for two-party preferred daily track
## generates picture summarizing Gibbs sampler output
##
## simon jackman, dept of political science
## stanford university, october 2005
#################################################################

if(exists("foo"))
  rm(foo)
foo <- list()

foo$y <- data$coalition2PP/100
var <- foo$y*(1-foo$y)/data$sampleSize
foo$prec <- 1/var
foo$date <- data$date - min(data$date) + 1
foo$org <- data$org
foo$NPOLLS <- length(data$y)
foo$NPERIODS <- length(min(data$date):282)
foo$alpha <- c(rep(NA,length((min(data$date)):282)-1),
               .5274)  ## actual 2PP on last day

## write content of object foo back to top level directory
for(i in 1:length(foo))
  assign(x=names(foo)[i],
         value=foo[[i]])
dump(list=names(foo))   ## dump
rm(list=names(foo))     ## now clean-up

## run jags job in batch mode
system("jags jags.cmd")

## read JAGS output back into R
library(coda)
alpha <- read.jags()
house <- alpha[,116:120]
sigma <- alpha[,115]
alpha <- alpha[,1:114]


z <- alpha[,c(113,115:120)]
z <- unclass(z)
pdf(file="traceplots.pdf",
    width=8,
    height=6)
par(bg="black",fg="white")
plot.ts(z,xlab="Iterations")
        
dev.off()

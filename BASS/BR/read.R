load("rollcall_lulaII.Rdata")

data <- rollcall_lulaII
rm(rollcall_lulaII)

summary(data)
apply(data,2,function(x)length(unique(x)))

## vote breakdown by rollcall
tapply(data$vote,data$rollcall,table)

## O votes, by rollcall
ovoteRate <- tapply(data$vote,data$rollcall,
                    function(x){
                      100*mean(x=="O",na.rm=TRUE)
                    }
                    )
sort(ovoteRate)
hist(ovoteRate)

###############################################
## make a rollcall object from these data
library(pscl)
fullID <- paste(data$id,data$party,sep=":")
n <- length(unique(fullID))
m <- length(unique(data$rollcall))
y <- matrix(NA,n,m)

## put each vote in the right row and col
## of the rollcall matrix
theRows <- match(fullID,unique(fullID))
theCols <- match(data$rollcall,unique(data$rollcall))
for(i in 1:dim(data)[1]){
  y[theRows[i],theCols[i]] <- data$vote[i]
}
dimnames(y) <- list(unique(fullID),
                    unique(data$rollcall))

## recodes
y[y=="Y"] <- 1
y[y=="N"] <- 0
y[y %in% c(NA,"O","A")] <- NA

## how many votes per legislator (fullID)
votesPerLegislator <- apply(y,1,function(x)sum(!is.na(x)))

## data set, characteristics of legislators
ldata <- data.frame(id=data$id[match(unique(fullID),fullID)],
                    fullID=unique(fullID),
                    name=data$name[match(unique(fullID),fullID)],
                    party=data$party[match(unique(fullID),fullID)],
                    nvotes=votesPerLegislator)

## rollcall object
rc <- rollcall(data=y,
               legis.names=unique(fullID),
               legis.data=ldata,
               vote.names=unique(data$rollcall))

## run ideal (!)
id1 <- ideal(rc,d=1,
             burnin=10E3,
             thin=100,
             maxiter=510E3,
             store.item=TRUE,
             normalize=TRUE,
             verbose=TRUE)

## dump rc and id1 to file
save("rc",file="rc.rda")
sace("id1",file="id1.rda")

## add posterior mean, sd and credible intervals
## to legislator data set, ldata
ldata$x <- id1$xbar
ldata$sd <- apply(id1$x,2,sd)
ldata$up <- apply(id1$x,2,quantile,.975)
ldata$lo <- apply(id1$x,2,quantile,.025)


## look at distribution of ideal point estimates
## by party, for the 8 biggest parties
bigParties <- c("PMDB","PT","DEM","PSDB",
                "PR","PP","PTB","PSB")
library(lattice)
histogram(~x | party,
          data=ldata,
          subset=party %in% bigParties)

xByParty <- tapply(ldata$x,
                   ldata$party,
                   mean)

## party switchers
dupids <- unique(ldata$id[duplicated(ldata$id)])

dups <- ldata[ldata$id %in% dupids,]
dups <- dups[order(dups$id),]

## make before/after data set, this only gets the 1st move
tmp <- by(dups,
          dups$id,
          function(obj){
            data.frame(id=unique(obj$id),
                       name=unique(obj$name),
                       fullidBefore=obj$fullID
                       fromParty=obj$party[1],
                       toParty=obj$party[2],
                       fromPoint=obj$x[1],
                       toPoint=obj$x[2])
          })

tmp <- do.call("rbind",tmp)


















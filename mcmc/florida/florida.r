## run jags on florida setup

system("jags florida.cmd")

library(coda)
mu <- read.jags()
plot(mu)
summary(mu)

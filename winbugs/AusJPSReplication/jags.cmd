model in kalman.bug
data in dumpdata.R
compile
initialize
update 1000
monitor alpha, thin(500)
monitor sigma, thin(500)
monitor house, thin(500)
update 25000
coda *
exit

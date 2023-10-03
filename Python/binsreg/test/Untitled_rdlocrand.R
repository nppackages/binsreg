###############################################################################
# rdlocrand: illustration file
# !version 1.0 21-Jun-2022
# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
###############################################################################

rm(list = ls())
options(width=200)

#install.packages("rdlocrand")
library(rdlocrand)

###############################################################################
## Load data
###############################################################################

data <- read.csv("/Users/rmasini/Library/CloudStorage/Dropbox/rdlocrand/R/rdlocrand_senate.csv")
dim(data)
names(data)

# Select predetermined covariates to be used for window selector

X  <-  cbind(data$presdemvoteshlag1,
             data$population/1000000,
             data$demvoteshlag1,
             data$demvoteshlag2,
             data$demwinprv1,
             data$demwinprv2,
             data$dopen,
             data$dmidterm,
             data$dpresdem)

# Assign names to the covariates

colnames(X) <-  c("DemPres Vote",
                  "Population",
                  "DemSen Vote t-1",
                  "DemSen Vote t-2",
                  "DemSen Win t-1",
                  "DemSen Win t-2",
                  "Open", "Midterm",
                  "DemPres")

# Running variable and outcome variable

R <- data$demmv
Y <- data$demvoteshfor2
D <- as.numeric(R>=0)


###############################################################################
## rdwinselect
###############################################################################

# Deprecated default options (Stata Journal 2016)

tmp <- rdwinselect(R,X,obsstep=2)

# Window selection with default options

tmp <- rdwinselect(R,X)

# Window selection with default options and symmetric windows

tmp <- rdwinselect(R,X,wasymmetric=TRUE)

# Window selection with specified window length and increments (replicate CFT)

tmp <- rdwinselect(R,X,wmin=.5,wstep=.125,reps=10000)

# Window selection using large sample approximation and plotting p-values

tmp <- rdwinselect(R,X,wmin=.5,wstep=.125,approx=TRUE,nwin=80,quietly=TRUE,plot=TRUE)


###############################################################################
## rdrandinf
###############################################################################

# Randomization inference using recommended window

tmp <- rdrandinf(Y,R,wl=-.75,wr=.75)


# Randomization inference using recommended window, all statistics

tmp <- rdrandinf(Y,R,wl=-.75,wr=.75,statistic='all')


# Randomization inference using recommended window using rdwinselect

tmp <- rdrandinf(Y,R,statistic='all',covariates=X,wmin=.5,wstep=.125,rdwreps=10000)


# Randomization inference using recommended window, linear adjustment

tmp <- rdrandinf(Y,R,wl=-.75,wr=.75,statistic='all',p=1)


# Randomization inference under interference

tmp <- rdrandinf(Y,R,wl=-.75,wr=.75,interfci=.05)


###############################################################################
## rdsensitivity
###############################################################################

tmp <- rdsensitivity(Y,R,wlist=seq(.75,10,by=.25),tlist=seq(0,20,by=1))

# Replicate contour plot

xaxis <- tmp$wlist
yaxis <- tmp$tlist
zvalues <- tmp$results
filled.contour(xaxis,yaxis,t(zvalues),
               xlab='window',ylab='treatment effect',
               key.title=title(main = 'p-value',cex.main=.8),
               levels=seq(0,1,by=.01),col=gray.colors(100,1,0))

# Obtain 95 percent confidence interval for window [-.75 ; .75]

tmp <- rdsensitivity(Y,R,wlist=seq(.75,2,by=.25),tlist=seq(0,20,by=1),ci=c(-0.75,0.75))
tmp$ci

# rdsensitivity to calculate CI from within rdrandinf

tmp <- rdrandinf(Y,R,wl=-.75,wr=.75,ci=c(.05,seq(3,20,by=1)))


###############################################################################
## rdrbounds
###############################################################################

tmp <- rdrbounds(Y,R,expgamma=c(1.5,2,3),wlist=c(.5,.75,1),reps=1000)
tmp$lower.bound
tmp$upper.bound

# Bernoulli and fixed margins p-values

tmp <- rdrbounds(Y,R,expgamma=c(1.5,2,3),wlist=c(.5,.75,1),reps=1000,fmpval=TRUE)
tmp$lower.bound
tmp$upper.bound

###############################################################################
## rdrandinf with eval options
###############################################################################

ii <- (R>=-.75) & (R<=.75) & !is.na(Y) & !is.na(R)
m0 <- mean(R[ii & D==0],na.rm=TRUE)
m1 <- mean(R[ii & D==1],na.rm=TRUE)
tmp <- rdrandinf(Y,R,wl=-.75,wr=.75,p=1,evall=m0,evalr=m1)

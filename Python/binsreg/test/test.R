rm(list=ls(all=TRUE))
library(ggplot2)

# Load the all functions in the binsreg package
for (fun in dir()) source(fun)
library(matrixStats)
library(splines)
library(sandwich)
library(quantreg)

# Load the binsreg package
# library(binsreg)


setwd("~/Dropbox/binsreg/Python/binsreg/test")
data <- read.csv("binsreg_sim.csv", sep=",")
y <- data$y; x <- data$x; w <- data$w; t <- data$t; id <- data$id; d <- data$d
#summary(data)

####################################
#######TEST BINSREGSELECT ##########
####################################

# out = binsregselect(y, x)
# summary(out)

# out = binsregselect(y, x, w, savegrid = TRUE)
# summary(out)

####################################
############# BINSREG ##############
####################################

# Default syntax
# est <-  binsreg(y, x, w, nbins = 21)
# summary(out)
# print(out$data.plot$`Group Full Sample`$data.dots)

# Evaluate the estimated function at median of w rather than the mean
# est <- binsreg(y, x, w, at="median", data=data)
# summary(est)
# print(est$data.plot$`Group Full Sample`$data.dots)

# Setting quantile-spaced bins to J=20, add a linear fit
# est <- binsreg(y, x, w, data=data, nbins=20, polyreg=1)
# summary(est)

# Adding lines
#est <- binsreg(y, x, w, nbins=20, line=c(3,3))
# summary(est)
# print(out$data.plot$`Group Full Sample`$data.line)

# Adding a CI
# est <- binsreg(y, x, w, nbins=20, ci=c(3,3))
# summary(est)
# print(est$data.plot$`Group Full Sample`$data.ci)

# Adding a CB
# est <- binsreg(y, x, w, nbins=20, cb=c(3,3))
# summary(est)
# print(est$data.plot$`Group Full Sample`$data.cb)

# Adding a polyreg
# est <- binsreg(y, x, w, nbins=20, polyreg = 4)
# summary(est)
# print(est$data.plot$`Group Full Sample`$data.poly)

# Adding line, ci, cb and polyreg together
# est <- binsreg(y, x, w, data=data, nbins=20, line=c(3,3), ci=c(3,3), cb=c(3,3), polyreg=4)
# summary(est)

# Comparison by groups
#est <- binsreg(y, x, w, nbins=20, by=t, line=c(3,3), cb=c(3,3),
#                bycolors=c("blue", "red"), bysymbols=c(19,17))

####################################
############# BINSGLM ##############
####################################

# out = binsglm(d, x, w, nbins = 10, family = binomial())
# summary(out)

# Plot the function in the inverse link (logistic) function rather than the conditional probability
# ut =  binsglm(d, x, w, nbins = 17, family = binomial(), nolink = T)
# summary(out)


####################################
############# BINSQreg ##############
####################################

# out = binsqreg(y, x, w, nbins=10, quantile = 0.25)
# summary(out)

# use bootstrap-based VCE
# out = binsqreg(y, x, w, data=data, quantile=0.25, ci=c(3,3), vce="boot", R=100)
# summary(out)

########################################
############# BINSTEST ##############
########################################
# basic syntax: linearity? (default method: least squares regression)
# bstest <- binstest(y, x, w, nbins = 21, testmodelpoly=1)
# summary(bstest)

# Test many things simultaneously
# bstest <- binstest(y, x, w, nbins=20, testshaper=c(-2,0), testshapel=4,
#                  testmodelpoly=1, nsims=1000, simsgrid=30)
# summary(bstest)

# Quantile regression
# Median regression: linear?
# bstest <- binstest(y, x, w, estmethod="qreg", quantile=0.5, testmodelpoly=1)
# summary(bstest)

# Logistic Regression
# Shape restriction test: increasing?
# bstest = binstest(d, x, w, estmethod="glm", family=binomial(), deriv=1, nbins=20, testshaper=0)
# summary(bstest)
 
 
 ########################################
 ########### BINSPWC ####################
 ########################################
 
 # Basic syntax
 # bsc <- binspwc(y, x, w, by=t)
 # summary(bsc)
 
 # GLM syntax with pre-specified # of bins
 # bsc <-binspwc(d, x, w, bynbins =c(5,4), by=t, family=binomial())
 # summary(bsc)
 
 # Compare quantile regression functions
 bsc <- binspwc(y, x, w, bynbins = c(5,4), by=t, estmethod="qreg", quantile=0.2)
 summary(bsc)



################################################################################
# Binsreg: illustration file
# Authors: M. D. Cattaneo, R. Crump, M. Farrell and Y. Feng
# Last update: 19-MAR-2019
rm(list=ls(all=TRUE))
library(binsreg); library(ggplot2)


#######################################
###### Read the same data #############
####### used for STATA ################
#######################################

data <- read.csv("binsreg_sim.csv", sep=",")
y <- data$y; x <- data$x; w <- data$w; t <- data$t
summary(data)

####################################
############# BINSREG ##############
####################################

# Default syntax
est <- binsreg(y, x, w)
est$bins_plot

# Setting quantile-spaced bins to J=20, add a linear fit
est <- binsreg(y, x, w, nbins=20, polyreg=1)
est$bins_plot

# Adding lines, ci, cb, polyreg
est <- binsreg(y, x, w, nbins=20, line=c(3,3))
est$bins_plot

est <- binsreg(y, x, w, nbins=20, line=c(3,3), ci=c(3,3))
est$bins_plot

est <- binsreg(y, x, w, nbins=20, line=c(3,3), ci=c(3,3), cb=c(3,3))
est$bins_plot

est <- binsreg(y, x, w, nbins=20, line=c(3,3), ci=c(3,3), cb=c(3,3), polyreg=4)
est$bins_plot

# VCE option, factor vars, ggplot object modification
est <- binsreg(y, x, cbind(w, t), dots=c(0,0), line=c(3,3), ci=c(3,3),
               cb=c(3,3), polyreg=4, vce="HC1", cluster=data$id)

est$bins_plot + ggtitle("Binned Scatter Plot") +
             theme(plot.title=element_text(hjust=0.5, vjust=0.5, face='bold'))

# Comparision by groups
est <- binsreg(y, x, w, by=t, line=c(3,3), cb=c(3,3), 
               bycolors=c("blue", "red"), bysymbols=c(19,17))

########################################
############# BINSREGTEST ##############
########################################
# basic syntax: linearity?
bstest <- binsregtest(y, x, w, testmodelpoly=1)
summary(bstest)

# Alternative: save parametric fit in another data frame or matrix
# If not available, first create by using binsregselect
bins <- binsregselect(y,x,w, simsgrid=30, savegrid = T)
grid <- bins$data.grid
colnames(grid)[4] <- "w"
ols <- lm(y~x+w)
ols.pred <- predict(ols, newdata=grid)
bstest <- binsregtest(y, x, w, testmodelparfit=cbind(grid[1], ols.pred))
summary(bstest)

# Shape restriction test: increasing?
bstest <- binsregtest(y, x, w, deriv=1, nbins=20, testshaper=0)
summary(bstest)

# Test many things simultaneously
bstest <- binsregtest(y, x, w, nbins=20, testshaper=c(-2,0), testshapel=4,
                      testmodelpoly=1, nsims=1000, simsgrid=30)
summary(bstest)

########################################
########### BINSREGSELECT ##############
########################################
# basic syntax
bins <- binsregselect(y,x,w)
summary(bins)

# J ROT specified manually and require evenly-spaced binning
bins <- binsregselect(y,x,w, nbinsrot=20, binspos="es")
summary(bins)

# Save grid for prediction purpose
bins <- binsregselect(y,x,w, simsgrid=30, savegrid = T)
grid <- bins$data.grid

# Extrapolating the optimal number of bins to the full sample
bins <- binsregselect(y,x,w, useeffn=1000, subset=(t==0))
summary(bins)

################################################################
# All mass point checks can be shut down to speed up computation
################################################################
est <- binsreg(y,x,w,masspoints = "off")
summary(est)

##########################################################
# BINSREG, integrated with BINSREGTEST and BINSREGSELECT #
##########################################################
# Automatic bin selection, binscatter plotting, and testing
est <- binsreg(y, x, w, dots=c(0,0), line=c(3,3), ci=c(3,3), cb=c(3,3),
               polyreg=4, testmodelpoly=1, testshapel=4)
summary(est)

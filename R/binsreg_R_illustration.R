################################################################################
# Binsreg: illustration file
# Authors: M. D. Cattaneo, R. Crump, M. Farrell and Y. Feng
# Last update: 07/23/2021
rm(list=ls(all=TRUE))
library(binsreg); library(ggplot2)


#######################################
###### Read the same data #############
####### used for STATA ################
#######################################

data <- read.csv("binsreg_sim.csv", sep=",")
y <- data$y; x <- data$x; w <- data$w; t <- data$t; id <- data$id; d <- data$d
summary(data)

####################################
############# BINSREG ##############
####################################

# Default syntax
est <- binsreg(y, x, w, data=data)
est$bins_plot

# Evaluate the estimated function at median of w rather than the mean
est <- binsreg(y, x, w, at="median", data=data)
est$bins_plot

# Use a formula containing a factor variable (t), evaluate the estimated function at w=0.2 and t=1 saved in a data frame
evalcovar <- data.frame(w=0.2, t=1)
est <- binsreg(y, x, w=~w+factor(t), data=data, at=evalcovar)
est$bins_plot

# Setting quantile-spaced bins to J=20, add a linear fit
est <- binsreg(y, x, w, data=data, nbins=20, polyreg=1)
est$bins_plot

# Adding lines, ci, cb, polyreg
est <- binsreg(y, x, w, data=data, nbins=20, line=c(3,3))
est$bins_plot

est <- binsreg(y, x, w, data=data, nbins=20, line=c(3,3), ci=c(3,3))
est$bins_plot

est <- binsreg(y, x, w, data=data, nbins=20, line=c(3,3), ci=c(3,3), cb=c(3,3))
est$bins_plot

est <- binsreg(y, x, w, data=data, nbins=20, line=c(3,3), ci=c(3,3), cb=c(3,3), polyreg=4)
est$bins_plot

# Specify y, x, w variables directly (without specifying a data frame), VCE option, ggplot object modification
est <- binsreg(y, x, cbind(w, t), dots=c(0,0), line=c(3,3), ci=c(3,3),
               cb=c(3,3), polyreg=4, vce="HC1", cluster=data$id)
# Modify other ggplot features
est$bins_plot + ggtitle("Binned Scatter Plot") +
                theme(plot.title=element_text(hjust=0.5, vjust=0.5, face='bold'))

# CI and CB: alternative formula for standard errors (nonparametric component only)
est <- binsreg(y, x, cbind(w, t), dots=c(0,0), line=c(3,3), ci=c(3,3),
               cb=c(3,3), polyreg=4, vce="HC1", cluster=data$id, asyvar=T)
est$bins_plot

# Comparison by groups
est <- binsreg(y, x, w, data=data, by=t, line=c(3,3), cb=c(3,3),
               bycolors=c("blue", "red"), bysymbols=c(19,17))

# Shut down all mass point checks to speed computation
est <- binsreg(y, x, w, data=data, masspoints="off")


########################################
############# BINSQREG #################
########################################

# 0.25 quantile
binsqreg(y, x, w, data=data, quantile=0.25)

# use bootstrap-based VCE
binsqreg(y, x, w, data=data, quantile=0.25, ci=c(3,3), vce="boot", R=100)

# Estimate 0.25 and 0.75 quantiles and combine them with the results from binsreg
est.25   <- binsqreg(y, x, quantile=0.25, line=c(3, 3))
dat.25   <- est.25$data.plot$`Group Full Sample`$data.line
dat.25["id"] <- "0.25 quantile"
est.75   <- binsqreg(y, x, quantile=0.75, line=c(3, 3))
dat.75   <- est.75$data.plot$`Group Full Sample`$data.line
dat.75["id"] <- "0.75 quantile"
est.mean <- binsreg(y, x, line=c(3, 3), cb=c(3, 3))
dat.mean.dots <- est.mean$data.plot$`Group Full Sample`$data.dots
dat.mean.line <- est.mean$data.plot$`Group Full Sample`$data.line
dat.mean.cb   <- est.mean$data.plot$`Group Full Sample`$data.cb
dat.mean.dots["id"] <- dat.mean.line["id"] <- "mean"

fig <- ggplot() + geom_point(data=dat.mean.dots, aes(x=x, y=fit, colour=id), size=2) +
                  geom_line(data=dat.mean.line, aes(x=x, y=fit, colour=id)) +
                  geom_ribbon(data=dat.mean.cb, aes(x=x, ymin=cb.l, ymax=cb.r), alpha=0.2, fill="navy") +
                  geom_line(data=dat.25, aes(x=x, y=fit, colour=id)) +
                  geom_line(data=dat.75, aes(x=x, y=fit, colour=id)) +
                  theme_bw() + labs(x="X", y="Y") +
                  scale_color_manual(name="", values = c("maroon", "darkgreen","navy"),
                                     guide=guide_legend(override.aes = list(
                                     linetype=rep("solid", 3), shape=c(NA, NA, 19))))
fig

########################################
############# BINSGLM ##################
########################################

# Basic syntax: binscatter logistic regression
est <- binsglm(d, x, w, data=data, family = binomial())
est$bins_plot

# Plot the function in the inverse link (logistic) function rather than the conditional probability
est <- binsglm(d, x, w, data=data, family = binomial(), nolink = T)
est$bins_plot


########################################
############# BINSTEST ##############
########################################
# basic syntax: linearity? (default method: least squares regression)
bstest <- binstest(y, x, w, data=data, testmodelpoly=1)
summary(bstest)

# Alternative: save parametric fit in another data frame or matrix; use L2 metric rather than sup
# If not available, first create by using binsregselect
bins <- binsregselect(y,x,w, data=data, simsgrid=30, savegrid = T)
grid <- bins$data.grid
colnames(grid)[4] <- "w"
ols <- lm(y~x+w, data=data)
ols.pred <- predict(ols, newdata=grid)
bstest <- binstest(y, x, w, data=data, testmodelparfit=cbind(grid[1], ols.pred), lp=2)
summary(bstest)

# Shape restriction test: increasing?
bstest <- binstest(y, x, w, data=data, deriv=1, nbins=20, testshaper=0)
summary(bstest)

# Test many things simultaneously
bstest <- binstest(y, x, w, data=data, nbins=20, testshaper=c(-2,0), testshapel=4,
                      testmodelpoly=1, nsims=1000, simsgrid=30)
summary(bstest)

# Quantile regression
# Median regression: linear?
bstest <- binstest(y, x, w, data=data, estmethod="qreg", quantile=0.5, testmodelpoly=1)
summary(bstest)

# Logistic Regression
# Shape restriction test: increasing?
bstest <- binstest(d, x, w, data=data, estmethod="glm", family=binomial(), deriv=1, nbins=20, testshaper=0)
summary(bstest)

########################################
########### BINSPWC ####################
########################################

# Basic syntax
bsc <- binspwc(y, x, w, data=data, by=t)
summary(bsc)

# Compare quantile regression functions
bsc <- binspwc(y, x, w, data=data, by=t, estmethod="qreg", quantile=0.4)
summary(bsc)

########################################
########### BINSREGSELECT ##############
########################################
# Basic syntax
bins <- binsregselect(y, x, w, data=data)
summary(bins)

# J ROT specified manually and require evenly-spaced binning
bins <- binsregselect(y,x,w, data=data, nbinsrot=20, binspos="es")
summary(bins)

# Save grid for prediction purpose
bins <- binsregselect(y,x,w, simsgrid=30, savegrid = T)
grid <- bins$data.grid

# Extrapolating the optimal number of bins to the full sample
bins <- binsregselect(y,x,w, useeffn=1000, subset=(t==0))
summary(bins)

# Use a random subsample to select the number of bins for the full sample
bins <- binsregselect(y,x,w, randcut=0.3)
summary(bins)

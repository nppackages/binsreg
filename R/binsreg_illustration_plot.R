################################################################################
# Binsreg: illustration file for plot
# Authors: M. D. Cattaneo, R. Crump, M. Farrell, Y. Feng and Ricardo Masini
# Last update: Jul 22, 2024
################################################################################

rm(list=ls(all=TRUE))
library(binsreg); library(ggplot2)

######################################
###### Read the same data #############
####### used for STATA ################
#######################################

data <- read.csv("binsreg_sim.csv", sep=",")
summary(data)

########################################################
###### EXTERNAL PLOT USING BINSREG OUTPUT ##############
########################################################

# Run binsreg
est <- binsreg(y, x, w, data=data, line = c(3,3), ci=c(3,3), cb=c(3,3), polyreg=4)

# Extract the plotting information
result <- est$data.plot$`Group Full Sample`

# Create the figure to plot
fig <- ggplot() + labs(x='X',y ='Y')

# Add the dots
fig <- fig + geom_point(data=result$data.dots, aes(x=x, y=fit), color="blue", size=2, shape='o')

# Add the line
fig <- fig + geom_line(data=result$data.line, aes(x=x, y=fit), color="blue", size=0.5)

# Add the CI
fig <- fig + geom_errorbar(data=result$data.ci, aes(x=x, ymin=ci.l, ymax=ci.r), color="blue", size=0.5, width = 0.02, linetype='solid')

# Add the CB
fig <- fig + geom_ribbon(data=result$data.cb, aes(x=x, ymin=cb.l, ymax=cb.r), fill="blue", alpha=0.2)

# Add the polyreg
fig <- fig + geom_line(data=result$data.poly, aes(x=x, y=fit), color="red", size=0.5)

# Display the plot
print(fig)

############################
## Create a plot for CATE ##
## group indicator: t ######

# Approach 1: use binspwc() to plot point estimates and confidence bands directly
# Note: point estimates ("dots") and confidence bands are constructed on
#       an evenly-spaced grid over the common support of the different groups
est <- binspwc(y, x, w, by=t, data=data, plot=T, dotsngrid = 20)

# Approach 2: combine results from binsreg() (binsqreg(), or binsglm()) with those from binspwc()

# Step 1: get "dots" and "lines" based on 20 bins via binsreg()
# Note: here the same grid of evaluation points are used to get fits for two groups
#       by setting samebinsby=T, dotsgridmean=F and dotsgrid=1 
#       (so each dot is placed at the midpoint of the bin)
est <- binsreg(y, x, w, data=data, by=t, line=c(3,3),
               dotsgrid=1, dotsgridmean=F, nbins=20, samebinsby=T)
# extract the data for dots and lines 
dots.0 <- est$data.plot$`Group 0`$data.dots
dots.1 <- est$data.plot$`Group 1`$data.dots
line.0 <- est$data.plot$`Group 0`$data.line
line.1 <- est$data.plot$`Group 1`$data.line

# Step 2: prepare the data frames for plotting CATE
# point estimates (dots and line) of CATE
cate.dots <- data.frame(x=dots.0$x, fit=dots.1$fit-dots.0$fit)
cate.line <- data.frame(x=line.0$x, fit=line.1$fit-line.0$fit)

# Step 3: get "confidence band" via binspwc()
pwc <- binspwc(y, x, w, by=t, data=data, pwc=c(3,3), plot=T, simsgrid=100, nbins=20, samebinsby=T)
cate.cb <- pwc$data.plot$`Group 1 - Group 0`$data.cb 

# Step 4: generate the plot for CATE ("dots", "line" and "band")
fig <- ggplot() + labs(x='X',y ='Y')

# Add the dots
fig <- fig + geom_point(data=cate.dots, aes(x=x, y=fit), color="blue", size=2, shape='o')

# Add the line
fig <- fig + geom_line(data=cate.line, aes(x=x, y=fit), color="blue", size=0.5)

# Add the CB
fig <- fig + geom_ribbon(data=cate.cb, aes(x=x, ymin=cb.l, ymax=cb.r), fill="blue", alpha=0.2)
print(fig)

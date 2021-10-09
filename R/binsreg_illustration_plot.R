################################################################################
# Binsreg: illustration file for plot
# Authors: M. D. Cattaneo, R. Crump, M. Farrell, Y. Feng and Ricardo Masini
# Last update: Oct 9, 2021
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

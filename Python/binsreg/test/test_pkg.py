# Test suite for funs AFTER the package is built

from binsreg import *
import pandas as pd

# Inputs for testing
data = pd.read_csv('/Users/rmasini/Dropbox/binsreg/Python/binsreg/test/binsreg_sim.csv')
y = data.y
x = data.x
w = data.w
t = data.t
id = data.id
d = data.d

# ISSUES:
# 1) R quantile type =2 has no analog in Python
# 2) vce in python accpets only const, HC0-4 OR cluster but not both as in R
# 3) DPI and DPI-UKNOT in Binsregselect # bins, bias^2 and var does not match R exactly. Probably due to 1) and 2) ?!?!
# 4) vce = "boot" is binsqreg is not implemented. It appears in the description and illustration files
# 5) savegrid option in binsregselect append a matrix of zeros under wname. Why?
# 6) binsregselect uses cluster = x as default. Is that expected?
# 7) Vcov of quantile regrssion does not change with vce changes. Probably not yet implemented in statmodels

# TO DO:
# 1) Fix the link to the papers in PyPi description
# 2) Fix the links to the illustration file and simulated data

#############################
#### Test Binsregselect #####
#############################

out = binsregselect(y, x)
print(out)

out = binsregselect(y, x, w)
print(out)

out = binsregselect(y, x, w, savegrid = True)
print(out)

#############################
#### Test Binsreg ###########
#############################

# Default syntax
est =  binsreg(y, x, w, nbins = 21)
print(est)
print(est.data_plot[0].dots)

# Default syntax with data
est =  binsreg(y='y', x = 'x', w = 'w', data = data, nbins = 21)
print(est)
print(est.data_plot[0].dots)

# Evaluate the estimated function at median of w rather than the mean
est = binsreg(y, x, w, nbins = 21, at="median")
print(est)
print(est.data_plot[0].dots)

# Setting quantile-spaced bins to J=20, add a linear fit
est = binsreg(y, x, w, nbins=20, polyreg=1)
print(est)

# Adding line
est = binsreg(y, x, w, nbins=20, line=(3,3))
print(est)
print(est.data_plot[0].line)

# Adding a CI
est = binsreg(y, x, w, nbins=20, ci=(3,3))
print(est)
print(est.data_plot[0].ci)

# Adding a CB
est = binsreg(y, x, w, nbins=20,cb=(3,3))
print(est)
print(est.data_plot[0].cb)

# Adding a polyreg
est = binsreg(y, x, w, nbins=20, polyreg = 4)
print(est)
print(est.data_plot[0].poly)

# Adding line, ci, cb and polyreg together
est = binsreg(y, x, w, nbins=20, line = (3,3), ci =(3,3), cb =(3,3), polyreg = 4)
print(est)

# Comparison by groups
est = binsreg(y, x, w, nbins=20, by=t, line=(3,3), cb=(3,3),
                 bycolors=("blue", "red"), bysymbols=('o','^'))

est = binsreg(y, 'x', w, nbins=20, by='t', data=data, line=(3,3), cb=(3,3),
               bycolors=("blue", "red"), bysymbols=('o','^'))

#############################
#### Test Binsglm ###########
#############################

out = binsglm(d, x, w, nbins=10, dist = 'Binomial')
out = binsglm(d, x, w, nbins=10, dist = 'Binomial', link = 'logit')
print(out)

# Plot the function in the inverse link (logistic) function rather than the conditional probability
out = binsglm(d, x, w, nbins = 17, dist = 'Binomial', nolink = True)
print(out)

##############################
#### Test Binsqreg ###########
##############################

out = binsqreg(y, x, w, nbins=10, quantile =0.25)
print(out)

out = binsqreg(y, x, w, data=data, quantile=0.25, ci=(3,3))
print(out)

####################################
############ BINSTEST ##############
####################################

# basic syntax: linearity? (default method: least squares regression)
out = binstest(y, x, w,  nbins=21, testmodelpoly=1)
print(out)

# Test many things simultaneously
out = binstest(y, x, w, nbins=20, testshaper=(-2,0), testshapel=4,
                    testmodelpoly=1, nsims=1000, simsgrid=30)
print(out)

# Quantile regression - Median regression: linear?
out = binstest(y, x, w, estmethod="qreg", quantile=0.5, testmodelpoly=1)
print(out)

# Logistic Regression - Shape restriction test: increasing?
out = binstest(d, x, w, estmethod="glm", dist= 'Binomial', deriv=1, nbins=20, testshaper=0)
print(out)

#######################################
########## BINSPWC ####################
#######################################

# Basic syntax
bsc = binspwc(y, x, w, by=t)
print(bsc)

# Basic syntax with pre-specified # of bins
bsc = binspwc(y, x, w, bynbins =(5,4), by=t)
print(bsc)

# Compare GLM  with pre-specified # of bins
bsc = binspwc(d, x, w, bynbins =(5,4), by=t, estmethod="glm", dist='Binomial')
print(bsc)

# Compare quantile regression functions
bsc =  binspwc(y, x, w, bynbins = (5,4), by=t, estmethod="qreg", quantile=0.2)
print(bsc)
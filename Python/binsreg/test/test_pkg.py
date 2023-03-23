#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TEST SUIT FOR BINSREG

@author: rmasini
"""

import numpy as np
import pandas as pd
import time

####### To install the package from PyPi ##########
# (in terminal) pip install binsreg
# from binsreg import *
###################################################

####### To install the package from PyPi TEST ###############################
# (in terminal) pip install -i https://test.pypi.org/simple/ binsreg==X.X.X
# from binsreg import *
#############################################################################

####### To install the package locally ###############################
# (in terminal) pip uninstall binsreg
# (in terminal) python3 -m pip install ./dist/binsreg-X.X.X.tar.gz
# from binsreg import *
#######################################################################

########## To run from the source ##################################################
# (in terminal) pip uninstall binsreg
import sys
import numpy as np
sys.path.insert(0, '/Users/rmasini/Dropbox/binsreg/Python/binsreg/src/binsreg')
from binsregselect import binsregselect
from binsreg import binsreg
from binsglm import binsglm
from binsqreg import binsqreg
from binstest import binstest
from binspwc import binspwc
from funs import *
###################################################################################


# NEW Issues:
# Condition on the line, CI and CB output of binsreg will never evaluate to False because it is a matrix on NAs


# New To Do
# Pass optional arguments to the optimizers (reg, qreg and glm)
# Check covariance options (why qrea ones change the name?)

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
# 4) vce = "boot" in binsqreg is not implemented. It appears in the description and illustration files
# 5) savegrid option in binsregselect append a matrix of zeros under wname. Why?
# 6) binsregselect uses cluster = x as default. Is that expected?
# 7) Vcov of quantile regrssion does not change with vce changes. Probably not yet implemented in statmodels

# TO DO:
# 1) Fix the link to the papers in PyPi description
# 2) Fix the links to the illustration file and simulated data

#############################
#### Test Binsregselect #####
#############################

# Basic Syntax
# out = binsregselect(y, x)
# print(out)
# out.summary()

# Add covariates
# out = binsregselect(y, x, w)
# print(out)
# out.summary()

# Save the grid
# out = binsregselect(y, x, w, savegrid = True)
# print(out)
# out.summary()

# #############################
# #### Test Binsreg ###########
# #############################

# # Default syntax
# out  =  binsreg(y, x)
# print(out)
# out.summary()
# print(out.data_plot[0].dots)

# Default syntax with data
# out =  binsreg(y='y', x = 'x', w = 'w', data = data, nbins = 21)
# print(out)
# out.summary
# print(out.data_plot[0].dots)

# # Evaluate the estimated function at median of w rather than the mean
# out = binsreg(y, x, w, nbins = 21, at="median")
# print(out)
# out.summary
# print(out.data_plot[0].dots)

# # Setting quantile-spaced bins to J=20, add a linear fit
# est = binsreg(y, x, w, nbins=20, polyreg=1)
# print(est)

# # Adding line
# est = binsreg(y, x, w, nbins=20, line=(3,3))
# print(est)
# print(est.data_plot[0].line)

# # Adding a CI
# est = binsreg(y, x, w, nbins=20, ci=(3,3))
# print(est)
# print(est.data_plot[0].ci)

# # Adding a CB
# est = binsreg(y, x, w, nbins=20,cb=(3,3))
# print(est)
# print(est.data_plot[0].cb)

# # Adding a polyreg
# est = binsreg(y, x, w, nbins=20, polyreg = 4)
# print(est)
# print(est.data_plot[0].poly)

# # Adding line, ci, cb and polyreg together
# out = binsreg(y, x, w, nbins=20, line = (3,3), ci =(3,3), cb =(3,3), polyreg = 4)
# print(out)
# out.summary()

# Using multiple nbins
# out = binsreg(y,x,w, nbins = np.arange(20), pselect = 1, sselect=1)
# print(out)
# out.summary()

# Using pselect/sselect
# out = binsreg(y,x,w, nbins = 10, pselect = np.arange(5), sselect=np.arange(5))
# print(out)
# out.summary()

# # Comparison by groups
# out = binsreg(y, x, w, by=t, line=(3,3), cb=(3,3),
#                  bycolors=("blue", "red"), bysymbols=('o','^'))
# print(out)
# out.summary()

# out = binsreg(y, 'x', w, nbins=20, by='t', data=data, line=(3,3), cb=(3,3),
#               bycolors=("blue", "red"), bysymbols=('o','^'))
# out.summary()

# #############################
# #### Test Binsglm ###########
# #############################

# Basic syntax
# out = binsglm(d, x, dist = 'Binomial')
# print(out)
# out.summary()

# Altering optimizer parameters
# out = binsglm(d, x, dist = 'Binomial', maxiter = 500, tol = 1e-6)
# print(out)
# out.summary()

# Explict nbins and link
# out = binsglm(d, x, w, nbins=10, dist = 'Binomial', link = 'logit')
# print(out)
# out.summary()

# Plot the function in the inverse link (logistic) function rather than the conditional probability
# out = binsglm(d, x, w, nbins = 17, dist = 'Binomial', nolink = True)
# print(out)
# out.summary()

# ##############################
# #### Test Binsqreg ###########
# ##############################

# Basic Syntax
# out = binsqreg(y, x)
# print(out)
# out.summary()

# out = binsqreg(y, x, w, data=data, quantile=0.25, ci=(3,3))
# print(out)
# out.summary()

# Change Optimizer parameters
# out = binsqreg(y, x, max_iter = 5000, p_tol = 1e-6)
# print(out)
# out.summary()

# ####################################
# ############ BINSTEST ##############
# ####################################

# # basic syntax: linearity? (default method: least squares regression)
# out = binstest(y, x)
# print(out)
# out.summary()


# # Test many things simultaneously
# out = binstest(y, x, w, nbins=20, testshaper=(-2,0), testshapel=4,
#                     testmodelpoly=1, nsims=1000, simsgrid=30)
# print(out)
# out.summary()

# # Quantile regression - Median regression: linear?
# out = binstest(y, x, w, estmethod="qreg", quantile=0.5, testmodelpoly=1)
# print(out)
# out.summary()

# # Logistic Regression - Shape restriction test: increasing?
# out = binstest(d, x, w, estmethod="glm", dist= 'Binomial', deriv=1, nbins=20, testshaper=0)
# print(out)
# out.summary()

# #######################################
# ########## BINSPWC ####################
# #######################################

# Basic syntax
# out = binspwc(y, x)
# print(out)
# out.summary()

# Basic syntax 2
# out = binspwc(y, x, by=t)
# print(out)
# out.summary()

# Basic syntax with pre-specified # of bins
# out = binspwc(y, x, w, bynbins =(5,4), by=t)
# print(out)
# out.summary()

# Compare GLM  with pre-specified # of bins
# out = binspwc(d, x, w, bynbins =(5,4), by=t, estmethod="glm", dist='Binomial')
# print(out)
# out.summary()

#  Compare quantile regression functions
# out =  binspwc(y, x, w, bynbins = (5,4), by=t, estmethod="qreg", quantile=0.2)
# print(out)
# out.summary()
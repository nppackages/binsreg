################################################################################
# Binsreg: illustration file for Python
# Authors: M. D. Cattaneo, R. Crump, M. Farrell, Y. Feng and Ricardo Masini
# Last update: March 17, 2023
################################################################################

import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from plotnine import *
from binsreg import *


######################################
###### Read the same data #############
####### used for STATA ################
#######################################

data = pd.read_csv("/Users/yjfeng/Dropbox/binscatter/Python/binsreg_sim.csv")
data.describe().T

####################################
############# BINSREG ##############
####################################

# Default syntax
est = binsreg('y', 'x', 'w', data=data)
est.bins_plot

# Alternative: Specify y, x, w variables directly (without specifying a data frame)
y = data.y
x = data.x
w = data.w
est = binsreg(y, x, w)
est.bins_plot

# Evaluate the estimated function at median of w rather than the mean
est = binsreg('y', 'x', 'w', at="median", data=data)
est.bins_plot

# Setting quantile-spaced bins to J=20, add a linear fit
est = binsreg('y', 'x', 'w', data=data, nbins=20, polyreg=1)
est.bins_plot

# Adding lines, ci, cb, polyreg
est = binsreg('y', 'x', 'w', data=data, nbins=20, line=(3,3))
est.bins_plot

est = binsreg('y', 'x', 'w', data=data, nbins=20, line=(3,3), ci=(3,3))
est.bins_plot

est = binsreg('y', 'x', 'w', data=data, nbins=20, line=(3,3), ci=(3,3), cb=(3,3))
est.bins_plot

est = binsreg('y', 'x', 'w', data=data, nbins=20, line=(3,3), ci=(3,3), cb=(3,3), polyreg=4)
est.bins_plot

#  VCE option, ggplot object modification
est = binsreg('y', 'x', ['w', 't'], data=data, dots=(0,0), line=(3,3), ci=(3,3),
               cb=(3,3), polyreg=4, vce='HC1', cluster=data.id)

# Modify other ggplot features
est.bins_plot + ggtitle('Binned Scatter Plot') + theme(plot_title=element_text(hjust=0.5, vjust=0.5, face='bold'))

# CI and CB: alternative formula for standard errors (nonparametric component only)
est = binsreg('y', 'x', ['w', 't'], data=data, dots=(0,0), line=(3,3), ci=(3,3),
               cb=(3,3), polyreg=4, vce="HC1", cluster=data.id, asyvar=True)
est.bins_plot

# Comparison by groups
est = binsreg('y', 'x', 'w', data=data, by='t', line=(3,3), cb=(3,3),
               bycolors=("blue", "red"), bysymbols=('o','^'))

# Shut down all mass point checks to speed computation
est = binsreg('y', 'x', 'w', data=data, masspoints="off")

# Select the degree p and smoothness given the number of bins J
# Note: The selected p and s are used for point estimation;
#       p+1 and s+1 are used for confidence intervals/bands
est = binsreg('y', 'x', 'w', data=data, pselect=(1,2,3,4), nbins=20)

########################################
############# BINSQREG #################
########################################

# Estimate 0.25 quantile
est = binsqreg('y', 'x', 'w', data=data, quantile=0.25)
est.bins_plot

# Estimate 0.25  and line = (3,3)
est_25   = binsqreg('y', 'x', data = data, quantile=0.25, line=(3,3))

# Increase the number of iterations
out = binsqreg('y', 'x', 'w', data=data, quantile=0.25, line=(3,3), max_iter = 5000)

# Estimate 0.25 and 0.75 quantiles and combine them with the results from binsreg
dat_25   = est_25.data_plot[0].line
dat_25.insert(0,"id","0.25 quantile")
est_75   = binsqreg('y', 'x', data=data, quantile=0.75, line=(3,3), max_iter = 5000)
dat_75   = est_75.data_plot[0].line
dat_75.insert(0,"id","0.75 quantile")
est_mean = binsreg('y', 'x', data=data, line=(3,3), cb=(3,3))
dat_mean_dots = est_mean.data_plot[0].dots
dat_mean_line = est_mean.data_plot[0].line
dat_mean_cb   = est_mean.data_plot[0].cb
dat_mean_dots.insert(0,"id", "mean") 
dat_mean_line.insert(0,"id", "mean") 

fig = ggplot() + theme_bw() + labs(x="X", y="Y")
fig += theme(legend_position = (0.77,0.23),
             legend_title = element_blank(),
             legend_background = element_rect(fill = 'white'),
             legend_key = element_blank())
fig += geom_point(data=dat_mean_dots, mapping = aes(x ='x', y='fit', colour='id'), size=2)
fig += geom_line(data=dat_mean_line, mapping = aes(x ='x', y='fit', colour='id'))
fig += geom_ribbon(data=dat_mean_cb, mapping = aes(x ='x', ymin='cb_l', ymax='cb_r'), alpha=0.2, fill="navy")
fig += geom_line(data=dat_25, mapping = aes(x ='x', y='fit', colour='id'))
fig += geom_line(data=dat_75, mapping = aes(x ='x', y='fit', colour='id'))
fig += scale_color_manual(name="", values = ("navy", "maroon","darkgreen"),
                                     guide=guide_legend(override_aes = {    
                                        'linetype':["solid"]*3, 'shape':('o', 'None', 'None')}))
fig

# Control for the fitting process, e.g., the maximum number of iterations and tolerance
binsqreg('y', 'x', 'w', data=data, quantile=0.25, max_iter = 5000, p_tol = 1e-6)


########################################
############# BINSGLM ##################
########################################

# Basic syntax: binscatter logistic regression
est = binsglm('d', 'x', 'w', data=data, dist = 'Binomial')
est.bins_plot

# Plot the function in the inverse link (logistic) function rather than the conditional probability
est = binsglm('d', 'x', 'w', data=data, dist = 'Binomial', nolink = True)
est.bins_plot

# Control for the fitting process, e.g., the maximum number of iterations and tolerance
est = binsglm('d', 'x', 'w', data=data, dist = 'Binomial', maxiter = 100, tol = 1e-6)
est.bins_plot


########################################
############# BINSTEST ##############
########################################

# basic syntax: linearity? (default method: least squares regression)
bstest = binstest('y', 'x', 'w', data=data, testmodelpoly=1)
bstest.summary()

# Recommended strategy: test if 1st deriv=const
bstest = binstest('y', 'x', 'w', data=data, testmodelpoly=1, deriv=1)
bstest.summary()

# Alternative: save parametric fit in another data frame or matrix; use L2 metric rather than sup
# If not available, first create by using binsregselect
bins = binsregselect('y','x','w', data=data, simsgrid=30, savegrid = True)
grid = bins.data_grid
grid.insert(0,'w', np.zeros(grid.shape[0]))
ols = smf.ols('y ~ x + w', data).fit()
ols_pred = ols.predict(grid)
model = np.column_stack((grid.x, ols_pred))

bstest = binstest('y', 'x', 'w', data=data, testmodelparfit=model, lp=2, deriv=1)
bstest.summary()

# Shape restriction test: increasing?
bstest = binstest('y', 'x', 'w', data=data, deriv=1, nbins=20, testshaper=0)
print(bstest)

# Test many things simultaneously
bstest = binstest('y', 'x', 'w', data=data, nbins=20, testshaper=(-2,0), testshapel=4,
                      testmodelpoly=1, nsims=1000, simsgrid=30)
print(bstest)

# Quantile regression
# Median regression: linear?
bstest = binstest('y', 'x', 'w', data=data, estmethod="qreg", quantile=0.5, testmodelpoly=1)
print(bstest)

# Logistic Regression
# Shape restriction test: increasing?
bstest = binstest('d', 'x', 'w', data=data, estmethod="glm", dist='Binomial',
                     deriv=1, nbins=20, testshaper=0)
print(bstest)

########################################
########### BINSPWC ####################
########################################

# Basic syntax
bsc = binspwc('y', 'x', 'w', data=data, by='t')
print(bsc)

# Compare quantile regression functions
bsc = binspwc('y', 'x', 'w', data=data, by='t', estmethod="qreg", quantile=0.4)
print(bsc)

#######################################
########## BINSREGSELECT ##############
#######################################

# Basic syntax
bins = binsregselect('y','x','w', data=data)
print(bins)

# J ROT specified manually and require evenly-spaced binning
bins = binsregselect('y','x','w', data=data, nbinsrot=20, binspos="es")
print(bins)

# Save grid for prediction purpose
bins = binsregselect('y','x','w', data=data, simsgrid=30, savegrid = True)
grid = bins.data_grid

# Extrapolating the optimal number of bins to the full sample

bins = binsregselect('y','x','w', data=data, useeffn=1000, subset=(data.t==0))
print(bins)
bins.summary()

# Use a random subsample to select the number of bins for the full sample
bins = binsregselect('y','x','w', data=data, randcut=0.3)
print(bins)
bins.summary()

# Select the degree p and smoothness s
bins = binsregselect('y', 'x', 'w', data=data, nbins=20, pselect=(1,2,3,4))
bins.summary()

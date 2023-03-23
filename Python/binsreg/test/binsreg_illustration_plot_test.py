################################################################################
# Binsreg: illustration file for plot
# Authors: M. D. Cattaneo, R. Crump, M. Farrell, Y. Feng and Ricardo Masini
# Last update: March  17, 2023
################################################################################

import pandas as pd
import statsmodels.formula.api as smf
from plotnine import *

########## Run from the source ##################################################
# (in terminal) pip uninstall binsreg
import sys
sys.path.insert(0, '/Users/rmasini/Dropbox/binsreg/Python/binsreg/src/binsreg')
from binsregselect import binsregselect
from binsreg import binsreg
from binsglm import binsglm
from binsqreg import binsqreg
from binstest import binstest
from binspwc import binspwc
from funs import *
###################################################################################

######################################
###### Read the same data #############
####### used for STATA ################
#######################################

data = pd.read_csv("/Users/rmasini/Dropbox/binsreg/Python/binsreg/test/binsreg_sim.csv")
data.describe().T


########################################################
###### EXTERNAL PLOT USING BINSREG OUTPUT ##############
########################################################

# Run binsreg 
est = binsreg('y', 'x', 'w', data=data,line = (3,3), ci=(3,3), cb=(3,3), polyreg=4)

# Extract the plotting information
result = est.data_plot[0]

# Create the figure to plot
fig = ggplot() + labs(x='X',y ='Y')

# Add the dots
fig += geom_point(data=result.dots, mapping=aes(x='x', y='fit'), color="blue", size=2, shape='o')

# Add the line
fig += geom_line(data=result.line, mapping=aes(x='x', y='fit'), color="blue", size=0.5)

# Add the CI
fig += geom_errorbar(data=result.ci, mapping=aes(x='x', ymin='ci_l', ymax='ci_r'), color="blue", size=0.5, width = 0.02, linetype='solid')

# Add the CB
fig += geom_ribbon(data=result.cb, mapping=aes(x='x', ymin='cb_l', ymax='cb_r'), fill="blue", alpha=0.2)

# Add the polyreg
fig += geom_line(data=result.poly, mapping=aes(x='x', y='fit'), color="red", size=0.5)

# Display the plot
print(fig)
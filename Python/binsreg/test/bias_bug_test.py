import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from plotnine import *
from binsreg import *
from binsregselect_bias import binsregselect_bias
import sys
sys.path.append('/Users/rmasini/Dropbox/binsreg/Python/binsreg/src/binsreg/')
from funs import *

data = pd.read_csv("binsreg_sim.csv")
y = np.array(data.y).reshape(-1,1)
x = np.array(data.x).reshape(-1,1)
w = data.w
t = data.t
id = data.id
d = data.d

J = 4
knot = genKnot_qs(x, J)




# Test bias function
# bias(x, p, s, deriv, knot)
#p = 3
#s = 3
#deriv = 0
#b = bias(x, p, s, deriv, knot)

# Test genB function

# knot = genKnot_qs(x, J=4)

# out = genB(y, x, w = None, p=3, s=3, deriv=0, knot = knot, weights=None)



# Basic syntax

bins = binsregselect_bias('y', 'x', data=data, bins=(1,1))
print(bins)



import numpy as np
from numpy.linalg import solve
from scipy.stats import norm
from patsy import dmatrix
import warnings
import statsmodels.api as sm



# Replicate R quantile fucntion with "type =2" option.
def quantile_R(x,J):
    # m=0 for "type =2" in R
    a = np.sort(x)
    p = np.linspace(0,1,J+1)
    n = len(a)

    j = np.floor(n*p)
    g = (j==n*p)  # discontinuity index
    
    out = p.copy()
    out(g) = a(g)
    
    gamma = {False: 1, True: 0.5}.get(aleph==k)

    # Deal with case where it should be smallest value 
    if aleph < 1:
        return x[k-1]  # x[0]
    else:
        return (1.-gamma)*x[k-1] + gamma*x[k]


# Generate Design for x from the R version
def binsreg_spdes(eval, p, s, knot, deriv):
    n = len(eval)
    J = len(knot)-1
    if s == 0:
        # mimic STATA irecode
        pos = FindInterval(eval, knot)
        polyx = nanmat(n,p+1)
        h = np.diff(knot)
        eval_cen = (eval - (knot[:-1])[pos])/h[pos]
        for j in range(deriv,p+1):
            polyx[:,j] = (eval_cen**(j-deriv) * factorial(j)/factorial(j-deriv) / h[pos]**deriv).reshape(-1)
        P = nanmat(n,(p+1)*J)
        for j in range(J):
            P[:,j*(p+1):(j+1)*(p+1)] = polyx*(pos==j).reshape(-1,1)
    else:
        if len(knot) >= 3:
            ext_knot = knot
            # ext_knot = np.concatenate([np.repeat(knot[0],p+1), np.repeat(knot[1:J], p-s+1), np.repeat(knot[-1], p+1)])
        else:
            ext_knot = knot
            # ext_knot = np.concatenate([np.repeat(knot[0], p+1), np.repeat(knot[-1], p+1)])
        
        ext_knot[ext_knot<np.min(eval)] = np.min(eval)  # This line does not exist in R
        ext_knot[ext_knot>np.max(eval)] = np.max(eval) # This line does not exist in R
        P = dmatrix("bs(x, knots = ext_knot, degree=p, include_intercept=True) - 1", {"x": eval})   # The deriv does not appear (see the R code below)
        #P = splineDesign(knots = ext_knot, eval, ord = p+1, derivs = deriv)
    return P

# Using the patsy package
def binsreg_spdes_PATSY(x, p, s, knot, deriv = None):
    # the derivative is not used yet
    # n = len(x) and k = len(knot)
    # The output is a matrix with nrow = n and ncol = (k-2)*(p-s+1) + p+1
    return np.asarray(dmatrix("bs(x, knots=np.repeat(knot[1:-1], p-s+1), degree=p, include_intercept=True)-1"))
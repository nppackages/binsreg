#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Thu Mar 16 17:39:50 2023
# @author: Ricardo Masini

import numpy as np
import warnings
from binsreg.binsregselect import binsregselect
from binsreg.funs import *

def binstest(y, x, w=None, data=None, estmethod="reg", dist= None, link=None,
                quantile=None, deriv=0, at=None, nolink=False,
                testmodel=None, testmodelparfit=None, testmodelpoly=None,
                testshape=None, testshapel=None,
                testshaper=None, testshape2=None, lp=np.inf,
                bins=None, nbins=None, 
                pselect=None, sselect=None,
                binspos="qs", binsmethod="dpi", nbinsrot=None, randcut=None,
                nsims=500, simsgrid=20, simsseed=None,
                vce=None, cluster=None, asyvar=False,
                dfcheck=(20,30), masspoints="on", weights=None, subset=None,
                numdist=None, numclust=None, estmethodopt=None, **optmize):

    '''
    Data-Driven Nonparametric Shape Restriction and Parametric Model Specification Testing using Binscatter.

    Description
    -----------
    binstest implements binscatter-based hypothesis testing procedures for parametric functional
    forms of and nonparametric shape restrictions on the regression function of interest, following the results
    in \href{https://arxiv.org/abs/1902.09608}{Cattaneo, Crump, Farrell and Feng (2022a)}.
    If the binning scheme is not set by the user, the companion function binsregselect is used to implement binscatter in a
    data-driven way and inference procedures are based on robust bias correction.
    Binned scatter plots based on different methods can be constructed using the companion functions binsreg,
    binsqreg or binsglm.

    Parameters
    ----------
    y : array or str
        A vector of the outcome variable or a string if data is provided.

    x : array or str 
        A vector of the independent variable or a string if data is provided.

    w : array or str
        A vector or a matrix of control variables or a (list of) strings if data is provided.
    
    data : data frame
        Optional  pandas data frame containing variables used in the model.
    
    estmethod : str
        Estimation method. The default is estmethod="reg" for tests based on binscatter least squares regression. 
        Other options are "qreg" for quantile regression and "glm" for generalized linear regression. 
        If estmethod="glm", the option dist must be specified.

    dist : str
        Distribution of the error term. Possible values are "Gaussian" (default),"Binomial","Gamma","Poisson".
        For furher options refer to statsmodels package documentation.

    link : str
        Link function to be used together with dist. Possible values are "identity","logit","log","probit".
        The default is link = None, in that case the default link dor the specified dist is used. 
        For furher options and details refer to statsmodels package documentation.
    
    quantile : float
        The quantile to be estimated. A number strictly between 0 and 1.

    nolink : bool
        If true, the function within the inverse link function is reported instead of the conditional mean function for the outcome.

    deriv : int 
        Derivative order of the regression function for estimation, testing and plotting.
        The default is deriv=0, which corresponds to the function itself. 
        If nolink=True, deriv cannot be greater than 1.
    
    at: str
        Value of w at which the estimated function is evaluated.  The default is at="mean", which corresponds to
        the mean of w. Other options are: at="median" for the median of w, at="zero" for a vector of zeros.
        at can also be a vector of the same length as the number of columns of w (if w is a matrix) or a data frame
        containing the same variables as specified in w (when data is specified). Note that when at="mean" or at="median",
        all factor variables (if specified) are excluded from the evaluation (set as zero).
    
    testmodel : tuple or bool
        It sets the degree of polynomial and the number of smoothness constraints for parametric model specification testing.
        If testmodel=(p,s) is specified, a piecewise polynomial of degree p with s smoothness constraints is used.
        If testmodel=True or testmodel=None (default) is specified, testmodel=(1,1) is used unless the degree p and the smoothness s
        selection is requested via the option pselect (see more details in the explanation of pselect).
    
    testmodelparfit: data frame or matrix 
        Contains the evaluation grid and fitted values of the model(s) to be
        tested against.  The column contains a series of evaluation points
        at which the binscatter model and the parametric model of interest are compared with
        each other.  Each parametric model is represented by other columns, which must
        contain the fitted values at the corresponding evaluation points.
    
    testmodelpoly : int
        Degree of a global polynomial model to be tested against.
    
    testshape : tuple or bool
        It sets the degree of polynomial and the number of smoothness constraints for parametric shape restriction testing.
        If testshape=(p,s) is specified, a piecewise polynomial of degree p with s smoothness constraints is used.
        If testshape=True or testshape=None (default) is specified, testmshape=(1,1) is used unless the degree p and the smoothness s
        selection is requested via the option pselect (see more details in the explanation of pselect).timate of the regression
        function of interest for one-sided or two-sided testing.
    
    testshapel : array
        A vector of null boundary values for hypothesis testing. Each number a in the vector
        corresponds to one boundary of a one-sided hypothesis test to the left of the form
        H0: sup_x mu(x)<=a.
    
    testshaper : array
        A vector of null boundary values for hypothesis testing. Each number a in the vector
        corresponds to one boundary of a one-sided hypothesis test to the right of the form
        H0: inf_x mu(x)>=a.
    
    testshape2 : array
        A vector of null boundary values for hypothesis testing. Each number a in the vector
        corresponds to one boundary of a two-sided hypothesis test ofthe form
        H0: sup_x |mu(x)-a|=0.
    
    lp : float or numpy.inf
        A Lp metric used for (two-sided) parametric model specification testing and/or shape restriction testing. 
        The default is lp=numpy.inf, which corresponds to the sup-norm of the t-statistic. 
        Other options are lp= p for a positive real p.

    bins : tuple
        If bins=(p,s), it sets the piecewise polynomial of degree p with s smoothness constraints
        for data-driven (IMSE-optimal) selection of the partitioning/binning scheme.
        The default is bins=(0,0), which corresponds to the piecewise constant.
    
    nbins : int or bool or array
        Number of bins for partitioning/binning of x.  If nbins=True or nbins=None (default) is specified, the number
        of bins is selected via the companion command binsregselect in a data-driven, optimal way whenever possible.
        If a vector with more than one number is specified, the number of bins is selected within this vector via the companion command binsregselect.
    
    pselect : array
        vector of numbers within which the degree of polynomial p for point estimation is selected.
        Piecewise polynomials of the selected optimal degree p are used to construct dots or line if dots=True or line=True is specified,
        whereas piecewise polynomials of degree p+1 are used to construct confidence intervals or confidence band if ci=True or cb=True is specified.
        Note: To implement the degree or smoothness selection, in addition to pselect or sselect, nbins=# must be specified.
    
    sselect : array
        vector of numbers within which the number of constrains s for point estimation is selected.
        Piecewise polynomials of the selected optimal s smoothness constrains are used to construct dots or line if dots=True or line=True is specified,
        whereas piecewise polynomials of s+1 smoothness constrains are used to construct confidence intervals or confidence band if ci=True or cb=True is specified.
        If not specified, for each value p supplied in the option pselect, only the piecewise polynomial
        with the maximum smoothness is considered, i.e., s=p.
    
    binspos : array
        Position of binning knots. The default is binspos="qs", which corresponds to quantile-spaced
        binning (canonical binscatter). The other options is "es" for evenly-spaced binning.
    
    binsmethod : str
        Method for data-driven selection of the number of bins. The default is binsmethod="dpi",
        which corresponds to the IMSE-optimal direct plug-in rule.  The other option is: "rot"
        for rule of thumb implementation.

    nbinsrot : int
        Initial number of bins value used to construct the DPI number of bins selector.
        If not specified, the data-driven ROT selector is used instead.
        
    randcut : float
        Upper bound on a uniformly distributed variable used to draw a subsample for bins/degree/smoothness selection.
        Observations for which numpy.random.uniform()<=# are used. # must be between 0 and 1.
        By default, max(5,000, 0.01n) observations are used if the samples size n>5,000.

    nsims : int
        Number of random draws for constructing confidence bands. The default is nsims=500,
        which corresponds to 500 draws from a standard Gaussian random vector of size
        [(p+1)*J - (J-1)*s]. A larger number of draws is recommended to obtain the final results.
    
    simsgrid : int
        Number of evaluation points of an evenly-spaced grid within each bin used for evaluation of
        the supremum operation needed to construct confidence bands and hypothesis testing
        procedures. The default is simsgrid=20, which corresponds to 20 evenly-spaced
        evaluation points within each bin for approximating the supremum operator.
        A larger number of draws is recommended to obtain the final results.
    
    simsseed : int 
        Simulation seed.
    
    vce : str
        Procedure to compute the variance-covariance matrix estimator. Options are
            * "const" : homoskedastic variance estimator.
            * "HC0"   : heteroskedasticity-robust plug-in residuals variance estimator
                        without weights.
            * "HC1"   : heteroskedasticity-robust plug-in residuals variance estimator
                        with hc1 weights. Default.
            * "HC2"   : heteroskedasticity-robust plug-in residuals variance estimator
                        with hc2 weights.
            * "HC3"   : heteroskedasticity-robust plug-in residuals variance estimator
                        with hc3 weights.
    
    cluster: array 
        Cluster ID. Used for compute cluster-robust standard errors.
    
    asyvar : bool
        If true, the standard error of the nonparametric component is computed and the uncertainty related to control
        variables is omitted. Default is asyvar=FALSE, that is, the uncertainty related to control variables is taken into account.
    
    dfcheck : tuple
        Adjustments for minimum effective sample size checks, which take into account number of unique
        values of x (i.e., number of mass points), number of clusters, and degrees of freedom of
        the different statistical models considered. The default is dfcheck=(20, 30).
        See \href{https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2022_Stata.pdf}{Cattaneo, Crump, Farrell and Feng (2021b)} for more details.
    
    masspoints: str
        How mass points in x are handled. Available options:
            * "on"           : all mass point and degrees of freedom checks are implemented. Default.
            * "noadjust"     :  mass point checks and the corresponding effective sample size adjustments are omitted.
            * "nolocalcheck" :  within-bin mass point and degrees of freedom checks are omitted.
            * "off"          : "noadjust" and "nolocalcheck" are set simultaneously.
            * "veryfew"      : forces the function to proceed as if x has only a few number of mass points (i.e., distinct values).
                               In other words, forces the function to proceed as if the mass point and degrees of freedom checks were failed.
    
    weights : array
        An optional vector of weights to be used in the fitting process. Should be None or
        a numeric vector.
    
    subset : array
        Optional rule specifying a subset of observations to be used.

    numdist : int
        Number of distinct for selection. Used to speed up computation.
    
    numclust : int
        Number of clusters for selection. Used to speed up computation.

    **optimize : 
        Optional arguments to the GLM or QuantReg statsmodels.api optimizer. 
        For futher details refer to binsqreg and binsglm help files.

    Returns
    -------
    testshapeL : Results for testshapel, including: val , null boundary values;
                 stat, test statistics; and pval, p-value.

    testshapeR : Results for testshaper, including: val , null boundary values;
                 stat, test statistics; and pval, p-value.

    testshape2 : Results for testshape2, including: val , null boundary values;
                 stat, test statistics; and pval, p-value.

    testpoly : Results for testmodelpoly, including: val , null boundary values;
               stat, test statistics; and pval, p-value.

    testmodel : Results for testmodelparfir, including: val , null boundary values;
                stat, test statistics; and pval, p-value.

    imse_v_rot : Variance constant in IMSE, ROT selection.

    imse_b_rot : Bias constant in IMSE, ROT selection.

    imse_v_dpi : Variance constant in IMSE, DPI selection.

    imse_b_dpi : Bias constant in IMSE, DPI selection.

    options : A list containing options passed to the function, as well as total sample size (n),
              number of distinct values (Ndist) in x, number of clusters (Nclust), 
              and number of bins (nbins).

    See Also
    --------
    binsregselect, binsreg, binsqreg, binsglm, binspwc.
    
    Example
    ------- 
    >>> x = numpy.random.uniform(size = 500)
    >>> y = numpy.sin(x)+numpy.random.normal(size = 500)
    >>> out = binstest(y,x, testmodelpoly=1)
    >>> print(out)
    '''

    # param for internal use
    qrot = 2

    ####################
    ### prepare data ###
    ####################

    # extract x, y, w, by, cluster, weights and subset from data (if data is supplied)
    xname = 'x'
    yname = 'y'
    if data is not None:
        if isinstance(x, str):
            if x in data.columns:
                xname = x
                x = data[xname]
            else:
                raise Exception("x column not found in data.")
        if isinstance(y, str):
            if y in data.columns:
                yname = y
                y = data[yname]
            else:
                raise Exception("y column not found in data.")
        if w is not None:
            w = list(w)
            if all(isinstance(i, str) for i in w):
                if set(w).issubset(data.columns): w = data[w]
                else: raise Exception("w columns not found in data.")
        if cluster is not None:
            if isinstance(cluster, str):
                if cluster in data.columns: cluster = data[cluster]
                else: raise Exception("cluster column not found in data.")
        if weights is not None:
            if isinstance(weights, str):
                if weights in data.columns: weights = data[weights]
                else: raise Exception("weights column not found in data.")
        if subset is not None:
            if isinstance(subset, str):
                if subset in data.columns: subset = data[subset]
                else: raise Exception("subset column not found in data.")

    # Reshaping x, y, w, by, cluster and weights
    x = np.array(x).reshape(-1,1)
    y = np.array(y).reshape(-1,1)
    if w is not None:
        w = np.array(w).reshape(len(w),-1)
        nwvar = ncol(w)
    else: nwvar = 0
    if cluster is not None:
        cluster = np.array(cluster).reshape(len(cluster),-1)
    if weights is not None:
        weights = np.array(weights).reshape(len(weights),-1)

    # extract subset
    if subset is not None:
        subset = np.array(subset).reshape(-1)
        x = x[subset]
        y = y[subset]
        if w is not None: w = w[subset,:]
        if cluster is not None: cluster = cluster[subset]
        if weights is not None: weights = weights[subset]
        
    # remove missing values
    na_ok = complete_cases(x) & complete_cases(y)
    if w is not None:  na_ok = na_ok & complete_cases(w)
    if weights is not None: na_ok = na_ok & complete_cases(weights)
    if cluster is not None: na_ok = na_ok & complete_cases(cluster)

    y  = y[na_ok]
    x  = x[na_ok]
    if w is not None: w  = w[na_ok,:]
    if weights is not None: weights = weights[na_ok]
    if cluster is not None: cluster = cluster[na_ok]

    xmin = np.min(x)
    xmax = np.max(x)

    # evaluation point of w
    if at is None: at = "mean"

    # change method name if needed
    if dist is not None:
        if dist!="Gaussian" or link!="identity":
            estmethod = "glm"

    ### extract family obj using dist and link
    family_str = None
    if estmethod == "glm":
        if link is None: 
            family_str = f'sm.families.{dist}()'
        else: 
            family_str = f'sm.families.{dist}(sm.families.links.{link})'
        
        family = eval(family_str)
        linkinv = family.link.inverse
        linkinv_1 = family.link.inverse_deriv
        linkinv_2 = family.link.inverse_deriv2

    
    # define the estimation method
    is_qreg = False
    if estmethod=="reg":
        estmethod_name = "least squares regression"
    elif estmethod=="qreg":
        is_qreg = True
        estmethod_name = "quantile regression"
        if quantile is None: quantile = 0.5
    elif estmethod=="glm":
        estmethod_name = "generalized linear model"
    else:
        raise Exception("estmethod incorrectly specified.")
    if vce is None: vce = "HC1"
    vce_select = vce

    # analyze bins- and degree-related options
    tshaON = tmodON = False
    if ((testshapel is not None) or (testshaper is not None) or (testshape2 is not None)):
        tshaON = True
    if ((testmodelparfit is not None) or (testmodelpoly is not None)):
        tmodON = True
    if (isinstance(testshape, bool) and not testshape):
        testshape = None
    if (isinstance(testmodel, bool) and not testmodel):
        testmodel = None
    if (isinstance(testshape, bool)  and ((not tshaON) or (not isinstance(binspos, str)))):
        testshape = None
    if (isinstance(testmodel, bool)  and ((not tmodON) or (not isinstance(binspos, str)))):
        testmodel = None
    
    if not isinstance(binspos, str):
        if ((nbins is not None) or (pselect is not None) or (sselect is not None)):
            raise Exception("nbins, pselect or sselect incorrectly specified.")
    
    # 4 cases: select J, select p, user specify both, or an error

    try: 
        len_nbins = len(nbins)
    except:
        if isinstance(nbins,int): len_nbins = 1
        else: len_nbins = 0
    
    if bins is not None:
        if isinstance(nbins,int): binsp = binss = bins
        elif len(bins)==2: binsp, binsp = bins
        else: raise Exception("Bins not correctly specified.")
        plist = binsp
        slist = binss
        if  ((len_nbins==1) and (nbins!=0)):
            raise Exception("Bins or nbins is incorrectly specified.")
    else:
        plist = pselect
        slist = sselect
        binsp = binss = None

    try: 
        len_p = len(plist)
    except: 
        if isinstance(plist,int): len_p = 1
        else: len_p = 0

    try: 
        len_s = len(slist)
    except: 
        if isinstance(slist,int): len_s = 1
        else: len_s = 0

    if (len_p==1 and len_s==0):
        slist = plist
        len_s = 1
    if (len_p==0 and len_s==1):
        plist = slist
        len_p = 1

    # 1st case: select J
    selection = ""
    if isinstance(binspos,str):
        if (isinstance(nbins,bool) and nbins): selection = "J"
        if (len_nbins==1 and nbins==0): selection = "J"
        if nbins is None or len_nbins>1 or bins is not None: selection = "J"
    
    if selection=="J":
        if (len_p>1 or len_s>1):
            if nbins is None:
                raise Exception("nbins must be specified for degree/smoothness selection.")
            else:
                raise Exception("Only one p and one s are allowed to select # of bins.")
        if plist is None: plist = deriv
        if slist is None: slist = plist
        if bins is None:
            binsp = plist
            binss = slist
        len_p = 1
        len_s = 1
        if testshape is None: testshape = conc(plist+1, slist+1)
        if (isinstance(testshape,bool) and testshape):
            testshape = conc(plist+1, slist+1)
        if testmodel is None: testmodel = conc(plist+1, slist+1)
        if (isinstance(testmodel,bool) and testmodel):
            testmodel = conc(plist+1, slist+1)

    # 2nd case: select p (at least for one object)
    pselectOK = False
    if (selection!="J" and (isinstance(testshape,bool) or testshape is None) and tshaON):
        pselectOK  = True
    if (selection!="J" and (isinstance(testmodel,bool) or testmodel is None) and tmodON):
        pselectOK  = True
    if (pselectOK and len_nbins==1 and (len_p>1 or len_s>1)):
        selection = "P"

    # 3rd case: user specified
    if ((len_p<=1 and len_s<=1) and selection!="J"):
        selection = "U"
        if testshape is None:
            if (len_p==1 and len_s==1): testshape = conc(plist+1, slist+1)
            else: testshape = conc(deriv+1, deriv+1)
        if testmodel is None:
            if (len_p==1 and len_s==1): testmodel = conc(plist+1, slist+1)
            else: testmodel = conc(deriv+1, deriv+1)

    if (selection==""):
        raise Exception("Degree, smoothness, or # of bins not correctly specified.")

    ##################################################
    ######## Error Checking ##########################
    ##################################################

    if not isinstance(binspos, str):
        if np.min(binspos)<=xmin or np.max(binspos)>=xmax:
            raise Exception("Knots out of allowed range")
    else:
        if binspos!="es" and binspos!="qs":
            raise Exception("binspos incorrectly specified.")

    if bins is not None and len(bins)==2:
        if bins[0] < bins[1]:
            raise Exception("p<s not allowed.")

    if len(testshape)==2:
        if testshape[0] < testshape[1]:
            raise Exception("p<s not allowed.")

    if len(testmodel)==2:
        if testmodel[0] < testmodel[1]:
            raise Exception("p<s not allowed.")

    if len(testshape)>=1 and not isinstance(testshape, bool) and testshape[0]<deriv:
        raise Exception("p<deriv not allowed.")
    
    if len(testmodel)>=1 and not isinstance(testmodel, bool) and testmodel[0]<deriv:
        raise Exception("p<deriv not allowed.")
    if deriv < 0:
        raise Exception("derivative incorrectly specified.")
            
    if binsmethod!="dpi" and binsmethod!="rot":
        raise Exception("Bin selection method incorrectly specified.")
    
    if masspoints == "veryfew":
        raise Exception("veryfew not allowed for testing.")

    if (testshape is not None and len(testshape)>=1) and (bins is not None and len(bins)>=1):
        if not isinstance(testshape, bool) and testshape[0]<=bins[0]:
            warnings.warn("p for testing <= p for bin selection not suggested.")
  
    if (testmodel is not None and len(testmodel)>=1)  and (bins is not None and len(bins)>=1):
        if not isinstance(testmodel, bool) and testmodel[0]<=bins[0]:
            warnings.warn("p for testing <= p for bin selection not suggested.")
        
    if (nsims<2000 or simsgrid<50):
        print("Note: Setting at least nsims=2000 and simsgrid=50 is recommended to obtain the final results.")
    
    #################################################
    localcheck = massadj = True
    if masspoints=="on":
        localcheck = True
        massadj = True
    elif masspoints=="off":
        localcheck = False
        massadj = False
    elif masspoints=="noadjust":
        localcheck = True
        massadj = False
    elif masspoints=="nolocalcheck":
        localcheck = False
        massadj = True

    # effective size
    eN = N = len(x)
    Ndist = np.nan
    if massadj:
        if numdist is not None: Ndist = numdist
        else: Ndist = len(np.unique(x))
        eN = min(eN, Ndist)

    Nclust = np.nan
    if cluster is not None:
        if numclust is not None: Nclust = numclust
        else: Nclust = len(np.unique(cluster))
        eN = min(eN, Nclust)
    
    # prepare params
    tsha_p, tsha_s = testshape
    if isinstance(tsha_p,bool): tsha_p = None
    if tsha_s is not None and np.isnan(tsha_s): tsha_s = tsha_p

    tmod_p, tmod_s = testmodel
    if isinstance(tmod_p,bool): tmod_p = None
    if tmod_s is not None and np.isnan(tmod_s): tmod_s = tmod_p

    if selection=="J":
        if ((tsha_p is not None) and (tsha_p<=plist)):
            tsha_p = binsp+1 
            tsha_s = tsha_p
            warnings.warn("Degree for testshape has been changed. It must be greater than the degree for bin selection.")

        if ((tmod_p is not None) and (tmod_p<=plist)):
            tmod_p = binsp+1
            tmod_s = tmod_p
            warnings.warn("Degree for testmodel has been changed. It must be greater than the degree for bin selection.")
    
    if selection=="U":
        warnings.warn("Testing procedures are valid when nbins is much larger than the IMSE-optimal choice. Compare your choice with the IMSE-optimal one obtained by binsregselect().")

    nL = nR = nT = 0
    if testshapel is not None:
        if isinstance(testshapel,int): 
            testshapel = (testshapel,)
        nL = len(np.array(testshapel))
    if testshaper is not None:
        if isinstance(testshaper,int): 
            testshaper = (testshaper,)
        nR = len(np.array(testshaper))
    if testshape2 is not None:
        if isinstance(testshape2,int): 
            testshape2 = (testshape2,)
        nT = len(np.array(testshape2))
    ntestshape = nL + nR + nT

    if selection=="U":  selectmethod = "User-specified"
    else:
        if binsmethod=="dpi": selectmethod = "IMSE direct plug-in" 
        else: selectmethod = "IMSE rule-of-thumb"
        if selection=="J": selectmethod = selectmethod + " (select # of bins)"
        if selection=="P": selectmethod = selectmethod + " (select degree and smoothness)"


    knot = None
    if not isinstance(binspos, str):
        nbins = len(binspos)+1
        knot = np.concatenate([[xmin], binspos,[xmax]])
        position = "User-specified"
        es = False
    else:
        if binspos == "es":
            es = True
            position = "Evenly-spaced"
        else:
            es = False
            position = "Quantile-spaced"

    ### Bin selection if needed ###################
    imse_v_rot = imse_v_dpi = imse_b_rot = imse_b_dpi = np.nan
    if selection!="U":
        # check if rot can be implemented
        if binsp is None: binspcheck = 6
        else: binspcheck = binsp
        if nbinsrot is None:
            if eN <= dfcheck[0] + binspcheck + 1 + qrot:
                raise Exception("Too small effective sample size for bin selection.")  
        if np.isnan(Ndist): Ndist_sel = None
        else: Ndist_sel = Ndist
        if np.isnan(Nclust): Nclust_sel = None
        else: Ndist_sel = Nclust
        randcut1k = randcut
        if (randcut is None and  N>5000):
                randcut1k = max(5000/N, 0.01)
                warnings.warn("To speed up computation, bin/degree selection uses a subsample of roughly max(5,000, 0.01n) observations if the sample size n>5,000. To use the full sample, set randcut=1.")
        if (selection=="J"):
                binselect = binsregselect(y, x, w, deriv=deriv,
                                        bins=conc(binsp,binss), binspos=binspos, nbins=nbins,
                                        binsmethod=binsmethod, nbinsrot=nbinsrot,
                                        vce=vce_select, cluster=cluster, randcut=randcut1k,
                                        dfcheck=dfcheck, masspoints=masspoints, weights=weights,
                                        numdist=Ndist_sel, numclust=Nclust_sel)
                if np.isnan(binselect.nbinsrot_regul):
                    raise Exception("Bin selection fails.")
                if (binsmethod == "rot"):
                    nbins = binselect.nbinsrot_regul
                    imse_v_rot = binselect.imse_v_rot
                    imse_b_rot  = binselect.imse_b_rot
                elif (binsmethod == "dpi"):
                    nbins = binselect.nbinsdpi
                    imse_v_dpi = binselect.imse_v_dpi
                    imse_b_dpi = binselect.imse_b_dpi
                    if np.isnan(nbins):
                        warnings.warn("DPI selection fails. ROT choice used.")
                        nbins = binselect.nbinsrot_regul
                        imse_v_rot = binselect.imse_v_rot
                        imse_b_rot = binselect.imse_b_rot
        elif (selection=="P"):
            binselect = binsregselect(y, x, w, deriv=deriv,
                                        binspos=binspos, nbins=nbins,
                                        pselect=plist, sselect=slist,
                                        binsmethod=binsmethod, nbinsrot=nbinsrot,
                                        vce=vce_select, cluster=cluster, randcut=randcut1k,
                                        dfcheck=dfcheck, masspoints=masspoints, weights=weights,
                                        numdist=Ndist_sel, numclust=Nclust_sel)
            if np.isnan(binselect.prot_regul):
                raise Exception("Bin selection fails.")  
            if (binsmethod == "rot"):
                binsp = binselect.prot_regul
                binss =  binselect.srot_regul
                imse_v_rot = binselect.imse_v_rot
                imse_b_rot = binselect.imse_b_rot
            elif binsmethod == "dpi": 
                binsp = binselect.pdpi
                binss =  binselect.sdpi
                imse_v_dpi = binselect.imse_v_dpi
                imse_b_dpi = binselect.imse_b_dpi
                if np.isnan(binsp):
                    warnings.warn("DPI selection fails. ROT choice used.")
                    binsp  =  binselect.prot_regul
                    binss = binselect.srot_regul
                    imse_v_rot = binselect.imse_v.rot
                    imse_b_rot = binselect.imse_b_rot     
          
            if (isinstance(testshape,bool) or (testshape is None)):
                tsha_p = binsp+1
                tsha_s = binss+1
            elif tsha_p<=binsp:
                tsha_p = binsp+1
                tsha_s = tsha_p
                warnings.warn("Degree for testshape has been changed. It must be greater than the IMSE-optimal degree.")

            if (isinstance(testmodel,bool) or (testmodel is None)):
                tmod_p = binsp+1
                tmod_s = binss+1
            elif tmod_p<=binsp:
                tmod_p = binsp+1
                tmod_s = tmod_p
                warnings.warn("Degree for testmodel has been changed. It must be greater than the IMSE-optimal degree.")        

    # check eff size for testing
    tsha_fewobs = False
    if ntestshape!=0:
        if (nbins-1)*(tsha_p-tsha_s+1)+tsha_p+1+dfcheck[1]>=eN:
            tsha_fewobs = True
            warnings.warn("Too small eff. sample size for testing shape.")
    
    tmod_fewobs = False
    if testmodelparfit is not None or testmodelpoly is not None:
        if (nbins-1)*(tmod_p-tmod_s+1)+tmod_p+1+dfcheck[1]>=eN:
            tmod_fewobs = True
            warnings.warn("Too small eff. sample size for testing models.")

    # Generate knot
    if knot is None:
        if es: knot = genKnot_es(xmin, xmax, nbins)
        else: knot = genKnot_qs(x, nbins)

    # check local mass points
    knot = np.append(knot[0],np.unique(knot[1:]))
    if nbins != len(knot)-1:
        warnings.warn("Repeated knots. Some bins dropped.")
        nbins = len(knot)-1

    if localcheck:
        uniqmin = binsreg_checklocalmass(x, nbins, es, knot=knot) # mimic STATA
        if ntestshape != 0:
            if uniqmin < tsha_p+1:
                tsha_fewobs = True
                warnings.warn("Some bins have too few distinct values of x for testing shape.")
        if testmodelparfit is not None or testmodelpoly is not None:
            if uniqmin < tmod_p+1:
                tmod_fewobs = True
                warnings.warn("Some bins have too few distinct values of x for testing models.")
    
    # adjust w variables
        if w is not None: 
            if isinstance(at, str):
                if at=="mean":
                    eval_w = colWeightedMeans(x=w, w=weights)
                elif at=="median":
                    eval_w = colWeightedMedians(x=w, w=weights)
                elif at=="zero": eval_w = np.zeros(nwvar)
            else: eval_w = np.array(at).reshape(-1,1)
        else: eval_w = None

    # seed
    if simsseed is not None: np.random.seed(simsseed)
    # prepare grid for uniform inference
    x_grid = binsreg_grid(knot, simsgrid).eval

    #####################################
    ####### Shape restriction test ######
    #####################################
    stat_shapeL = pval_shapeL = np.nan
    stat_shapeR = pval_shapeR = np.nan
    stat_shape2 = pval_shape2 = np.nan
    if nL>0:
        stat_shapeL = nanmat(nL, 2)
        pval_shapeL = []
    if nR>0:
        stat_shapeR = nanmat(nR, 2)
        pval_shapeR = []
    if nT>0:
        stat_shape2 = nanmat(nT, 2)
        pval_shape2 = []
    
    if ntestshape != 0 and not tsha_fewobs:
        B = binsreg_spdes(x=x, p=tsha_p, s=tsha_s, knot=knot, deriv=0)
        k = ncol(B)
        P = cbind(B, w)
        model = binsreg_fit(y=y, x=P, weights=weights, family=family_str, is_qreg=is_qreg,
                                 quantile = quantile, cov_type = vce_select, cluster = cluster, **optmize)            
        beta = model.params[:k]
        basis_sha = binsreg_spdes(x=x_grid, p=tsha_p, s=tsha_s, knot=knot, deriv=deriv)

        if estmethod=="glm" and not nolink:
            fit_sha, se_sha = binsreg_pred(X=basis_sha, model=model, type="all",
                                            deriv=deriv, wvec=eval_w, avar=asyvar)
            basis_0 = binsreg_spdes(x=x_grid, p=tsha_p, s=tsha_s, knot=knot, deriv=0)
            fit_0 = binsreg_pred(basis_0, model, type = "xb", deriv=0, wvec=eval_w)[0]
            pred_sha_0  = linkinv_1(fit_0)

            if asyvar or deriv==0:
                se_sha  = pred_sha_0 * se_sha
                if deriv == 0: fit_sha = linkinv(fit_sha)
                if deriv == 1: fit_sha = pred_sha_0 * fit_sha
            else:
                basis_sha_1 = basis_sha
                if eval_w is not None:
                    basis_sha_0 = np.column_stack((basis_0, np.outer(np.ones(nrow(basis_0)), eval_w)))
                    basis_sha_1 = np.column_stack((basis_sha_1, np.outer(np.ones(nrow(basis_sha_1)), np.zeros(nwvar))))
                term1 = linkinv_2(fit_0).reshape(-1,1)*fit_sha.reshape(-1,1)*basis_sha_0
                term2 = pred_sha_0.reshape(-1,1)*basis_sha_1
                basis_all = term1 + term2
                fit_sha = pred_sha_0 * fit_sha
                se_sha  = binsreg_pred(basis_all, model=model, type="se", avar=True)[1]
        else:
            fit_sha, se_sha  = binsreg_pred(basis_sha, model, type = "all", deriv=deriv,
                                                 wvec=eval_w, avar=asyvar)

        pos = np.invert(np.isnan(beta))
        k_new = np.sum(pos)
        vcv_sha = model.cov_params()[:k_new,:k_new]

        for j in range(ntestshape):
            if j < nL:
                stat_shapeL[j,1] = 1
                stat_shapeL[j,0] = np.max((fit_sha - testshapel[j]) / se_sha)
            elif j < nL+nR:
                stat_shapeR[j-nL,1] = 2
                stat_shapeR[j-nL,0] = np.min((fit_sha - testshaper[j-nL]) / se_sha)
            else:
                stat_shape2[j-nL-nR,1] = 3
                if np.isfinite(lp):
                    stat_shape2[j-nL-nR,0] = np.max(np.abs((fit_sha - testshape2[j-nL-nR]) / se_sha))
                else:
                    stat_shape2[j-nL-nR,0] = np.mean(np.abs((fit_sha - testshape2[j-nL-nR]) / se_sha)**lp)**(1/lp)
    
        stat_shape = None
        if nL > 0: stat_shape = rbind(stat_shape, stat_shapeL)
        if nR > 0: stat_shape = rbind(stat_shape, stat_shapeR)
        if nT > 0: stat_shape = rbind(stat_shape, stat_shape2)
        
        # Compute p-val
        Sigma_root = lssqrtm(vcv_sha)
        num = np.matmul(basis_sha[:,pos], Sigma_root)
        denom_sha  = np.sqrt(np.sum(np.matmul(basis_sha[:,pos], vcv_sha) * basis_sha[:,pos],1))
        pval_shape = binsreg_pval(num, denom_sha, nsims, tstat=stat_shape, side=None, lp=lp)[0]
        if nL!=0:
            stat_shapeL = stat_shapeL[:,0]
            pval_shapeL = pval_shape[:nL]
        if nR!=0:
            stat_shapeR = stat_shapeR[:,0]
            pval_shapeR = pval_shape[nL:nL+nR]
        if nT!=0:
            stat_shape2 = stat_shape2[:,0]
            pval_shape2 = pval_shape[nL+nR:]
    
    ############################
    #### Specification Test ####
    ############################
    stat_poly = pval_poly = np.nan
    stat_mod = pval_mod = np.nan
    if testmodelpoly is not None:
        stat_poly = nanmat(1,2)
        pval_poly = []
    if testmodelparfit is not None:
        stat_mod = nanmat(ncol(testmodelparfit)-1,2)
        pval_mod = []

    if (testmodelparfit is not None or testmodelpoly is not None) and not tmod_fewobs:
        if tmod_p==tsha_p and tmod_s==tsha_s and "vcv_sha" in locals():
            exist_mod = True
            vcv_mod  = vcv_sha
            fit_mod  = fit_sha
            se_mod   = se_sha
            basis_mod = basis_sha
            denom_mod = denom_sha
        else:
            exist_mod = False
            B = binsreg_spdes(x=x, p=tmod_p, s=tmod_s, knot=knot, deriv=0)
            k = ncol(B)
            P = cbind(B, w)
            model = binsreg_fit(y=y, x=P, weights=weights, family=family_str, is_qreg=is_qreg,
                                 quantile=quantile,cov_type=vce_select,cluster=cluster, **optmize) 
            beta = model.params[:k]
            pos = np.invert(np.isnan(beta))
            k_new = np.sum(pos)
            if estmethod=="qreg": is_qreg = True
            else: is_qreg = False
            vcv_mod = model.cov_params()[:k_new,:k_new]
           
        ######################
        # Test poly reg
        if testmodelpoly is not None:
            if not exist_mod:
                basis_mod = binsreg_spdes(x=x_grid, p=tmod_p, s=tmod_s, knot=knot, deriv=deriv)

                if estmethod=="glm" and not nolink:
                    pred_mod_fit, pred_mod_se = binsreg_pred(X=basis_mod, model=model, type="all",
                                                            deriv=deriv, wvec=eval_w, avar=asyvar)
                    basis_0 = binsreg_spdes(x=x_grid, p=tmod_p, s=tmod_s, knot=knot, deriv=0)
                    fit_0 = binsreg_pred(basis_0, model, type = "xb", deriv=0, wvec=eval_w)[0]
                    pred_mod_0  = linkinv_1(fit_0)

                    if asyvar or deriv==0:
                        pred_mod_se  = pred_mod_0 * pred_mod_se
                        if deriv == 0: pred_mod_fit = linkinv(pred_mod_fit)
                        if deriv == 1: pred_mod_fit = pred_mod_0 * pred_mod_fit
                    else:
                        basis_mod_1 = basis_mod.copy()
                        if eval_w is not None:
                            basis_mod_0 = np.column_stack((basis_0, np.outer(np.ones((basis_0)), eval_w)))
                            basis_mod_1 = np.column_stack((basis_mod_1, np.outer(np.ones(nrow(basis_mod_1)), np.zeros(nwvar))))
                        basis_all = linkinv_2(fit_0)*pred_mod_fit*basis_mod_0 + pred_mod_0*basis_mod_1
                        pred_mod_fit = pred_mod_0 * pred_mod_fit
                        pred_mod_se  = binsreg_pred(basis_all, model=model, type="se", avar=True)[1]
                else:
                    fit_mod, se_mod = binsreg_pred(basis_mod, model, type = "all",  deriv=deriv, wvec=eval_w,
                                                avar=asyvar)

                denom_mod = np.sqrt(np.sum(np.matmul(basis_mod[:,pos],vcv_mod) * basis_mod[:,pos],1))

            # Run a poly reg
            x_p = nanmat(N, testmodelpoly+1)
            for j in range(testmodelpoly+1):  x_p[:,j] =(x**j).reshape(-1)
            P_poly = cbind(x_p, w)
            model_poly = binsreg_fit(y=y, x=P_poly, weights=weights, family=family_str, is_qreg=is_qreg,
                                 quantile = quantile, cov_type = vce_select, cluster = cluster, **optmize) 
            beta_poly = model_poly.params
            poly_fit = 0
            for j in range(deriv,testmodelpoly+1):
                poly_fit += x_grid**(j-deriv)*beta_poly[j]*factorial(j)/factorial(j-deriv)
            if eval_w is not None and deriv==0:
                 poly_fit += np.sum(eval_w * beta_poly[testmodelpoly+1:])

            stat_poly[0,1] = 3
            if np.isfinite(lp): stat_poly[0,0] = np.mean(np.abs((fit_mod - poly_fit) / se_mod)**lp)**(1/lp)
            else: stat_poly[0,0] = np.max(np.abs((fit_mod - poly_fit) / se_mod))

            Sigma_root   = lssqrtm(vcv_mod)
            num          = np.matmul(basis_mod[:,pos],Sigma_root)
            pval_poly    = binsreg_pval(num, denom_mod, nsims, tstat=stat_poly, side=None, lp=lp)[0]
            stat_poly    = stat_poly[:,0]

        ##################################
        # Test model in another data frame
        if testmodelparfit is not None:
            x_grid = testmodelparfit[:,0]
            basis_mod = binsreg_spdes(x=x_grid, p=tmod_p, s=tmod_s, knot=knot, deriv=deriv)

            if estmethod=="glm" and not nolink:
                pred_mod_fit, pred_mod_se = binsreg_pred(X=basis_mod, model=model, type="all",
                                                            deriv=deriv, wvec=eval_w, avar=asyvar)
                basis_0 = binsreg_spdes(x=x_grid, p=tmod_p, s=tmod_s, knot=knot, deriv=0)
                fit_0 = binsreg_pred(basis_0, model, type = "xb", deriv=0, wvec=eval_w)[0]
                pred_mod_0 = linkinv_1(fit_0)

                if asyvar or deriv==0:
                    pred_mod_se  = pred_mod_0 * pred_mod_se
                    if deriv == 0: pred_mod_fit = linkinv(pred_mod_fit)
                    if deriv == 1: pred_mod_fit = pred_mod_0 * pred_mod_fit
                else:
                    basis_mod_1 = basis_mod.copy()
                    if eval_w is not None:
                        basis_mod_0 = np.column_stack((basis_0, np.outer(np.ones(nrow(basis_0)),eval_w)))
                        basis_mod_1 = np.column_stack((basis_mod_1, np.outer(np.ones(nrow(basis_mod_1)), np.zeros(nwvar))))
                    basis_all = linkinv_2(fit_0)*pred_mod_fit*basis_mod_0 + pred_mod_0*basis_mod_1
                    pred_mod_fit = pred_mod_0 * pred_mod_fit
                    pred_mod_se  = binsreg_pred(basis_all, model=model, type="se",
                                                avar=True)[1]
            else:
                fit_mod, se_mod  = binsreg_pred(basis_mod, model, type = "all", deriv=deriv, 
                                            wvec=eval_w,avar=asyvar)
            
            denom_mod = np.sqrt(np.sum(np.matmul(basis_mod[:,pos],vcv_mod) * basis_mod[:,pos],1))

            for j in range(1,ncol(testmodelparfit)):
                stat_mod[j-1,1] = 3
                if np.isfinite(lp): stat_mod[j-1,0] = np.max(np.abs((fit_mod - testmodelparfit[:,j]) / se_mod))
                else: stat_mod[j-1,0] = np.mean(np.abs((fit_mod - testmodelparfit[:,j]) / se_mod)**lp)**(1/lp)

            Sigma_root = lssqrtm(vcv_mod)
            num = np.matmul(basis_mod[:,pos],Sigma_root)
            pval_mod = binsreg_pval(num, denom_mod, nsims, tstat=stat_mod, side=None, lp=lp)[0]
            stat_mod = stat_mod[:,0]
            
    ######################
    #######output#########
    ######################

    testshapeL = test_str(testshapel,stat_shapeL,pval_shapeL)
    testshapeR = test_str(testshaper,stat_shapeR,pval_shapeR)
    testshape2 = test_str(testshape2,stat_shape2,pval_shape2)
    testPoly = test_str(testmodelpoly,stat_poly,pval_poly)
    testModel = test_str(np.nan,stat_mod,pval_mod)
    
    opt = options_test(binsp=binsp, binss=binss, deriv=deriv,
                        testshape=(tsha_p, tsha_s), testmodel=(tmod_p, tmod_s),
                        binspos=position, binsmethod=selectmethod,
                        n=N, Ndist=Ndist, Nclust=Nclust,
                        nbins=nbins, knot=knot, lp=lp,
                        estmethod=estmethod_name, quantile=quantile, 
                        dist=dist,link=link)
    out = test_output(testshapeL,testshapeR,testshape2,testPoly,testModel,
                        imse_v_dpi=imse_v_dpi, imse_b_dpi=imse_b_dpi,
                        imse_v_rot=imse_v_rot, imse_b_rot=imse_b_rot,
                        options = opt)
    return out
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Sat Sep  4 17:39:50 2021
# @author: Ricardo Masini

import numpy as np
import warnings
from .binsregselect import binsregselect
# from binsregselect import binsregselect
from .funs import *
# from funs import *

def binspwc(y, x, w=None,data=None, estmethod="reg", dist=None, link=None,
                    quantile=None, deriv=0, at=None, nolink=False, by=None,
                    pwc=(3,3), testtype="two-sided", lp=np.inf,
                    bins=(2,2), bynbins=None, binspos="qs",
                    binsmethod="dpi", nbinsrot=None, samebinsby=False, randcut=None,
                    nsims=500, simsgrid=20, simsseed=None,
                    vce=None, cluster=None, asyvar=False,
                    dfcheck=(20,30), masspoints="on", weights=None, subset=None):

    '''
    Data-Driven Pairwise Group Comparison using Binscatter Methods.

    Description
    -----------
    binspwc implements hypothesis testing procedures for pairwise group comparison of binscatter estimators, following the
    results in \href{https://arxiv.org/abs/1902.09608}{Cattaneo, Crump, Farrell and Feng (2021a)}.
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

    by : array
        A vector containing the group indicator for subgroup analysis; both numeric and string variables
        are supported. When  by is specified, binsreg implements estimation and inference for each subgroup
        separately, but produces a common binned scatter plot. By default, the binning structure is selected for each
        subgroup separately, but see the option samebinsby below for imposing a common binning structure across subgroups.

    pwc : tuple
        pwc=(p,s) sets a piecewise polynomial of degree p with s smoothness constraints for testing the difference between groups. 
        The default is pwc=(3,3), which corresponds to a cubic B-spline estimate of the function of interest for each group.

    testtype : str
        Type of pairwise comparison test. The default is testtype="two-sided", which corresponds to a two-sided test of the form 
        H0: mu_1(x)=mu_2(x). Other options are: testtype="left" for the one-sided test form H0: mu_1(x)<=mu_2(x) and
        testtype="right" for the one-sided test of the form H0: mu_1(x)>=mu_2(x).
    
    lp : float or numpy.inf
        A Lp metric used for (two-sided) parametric model specification testing and/or shape restriction testing. 
        The default is lp=numpy.inf, which corresponds to the sup-norm of the t-statistic. 
        Other options are lp= p for a positive real p.

    bins : tuple
        Degree and smoothness for bin selection. The default is bins=c(2,2), which corresponds to a quadratic spline estimate.
    
    bynbins : array of int 
        A vector of the number of bins for partitioning/binning of x}, which is applied to the binscatter estimation for each group.
        If not specified, the number of bins is selected via the companion function binsregselect in a data-driven way whenever possible.
    
    binspos : array
        Position of binning knots. The default is binspos="qs", which corresponds to quantile-spaced
        binning (canonical binscatter). The other options is "es" for evenly-spaced binning, or
        a vector for manual specification of the positions of inner knots (which must be within the range of x.
    
    binsmethod : str
        Method for data-driven selection of the number of bins. The default is binsmethod="dpi",
        which corresponds to the IMSE-optimal direct plug-in rule.  The other option is: "rot"
        for rule of thumb implementation.

    nbinsrot : int
        Initial number of bins value used to construct the DPI number of bins selector.
        If not specified, the data-driven ROT selector is used instead.

    samebinsby : bool
        If true, a common partitioning/binning structure across all subgroups specified by the option by is forced.
        The knots positions are selected according to the option binspos and using the full sample. If \code{nbins}
        is not specified, then the number of bins is selected via the companion command binsregselect and
        using the full sample.
        
    randcut : float
        Upper bound on a uniformly distributed variable used to draw a subsample for bins selection.
        Observations for which numpy.random.uniform()<=# are used. # must be between 0 and 1.
    
    nsims : int
        Number of random draws for constructing confidence bands. The default is nsims=500,
        which corresponds to 500 draws from a standard Gaussian random vector of size
        [(p+1)*J - (J-1)*s].
    
    simsgrid : int
        Number of evaluation points of an evenly-spaced grid within each bin used for evaluation of
        the supremum (infimum or Lp metric) operation needed to construct confidence bands and hypothesis testing
        procedures. The default is simsgrid=20, which corresponds to 20 evenly-spaced
        evaluation points within each bin for approximating the supremum (infimum or Lp metric) operator.
    
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
        See \href{https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2021_Stata.pdf}{Cattaneo, Crump, Farrell and Feng (2021b)} for more details.
    
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
    
    Returns
    -------
    tstat : A matrix where each row corresponds to the comparison between two groups. The first column is the test statistic. 
            The second and third columns give the corresponding group numbers. The null hypothesis is mu_i(x)<=mu_j(x),
            mu_i(x)=mu_j(x), or mu_i(x)>=mu_j(x) for group i (given in the second column) and group j (given in the third column).
            The group number corresponds to the list of group names given by options_byvals.
    
    pval : A vector of p-values for all pairwise group comparisons.
    
    options : A list containing options passed to the function, as well as N_by (total sample size for each group),
              Ndist_by (number of distinct values in x for each group), Nclust_by (number of clusters for each group),
              nbins_by (number of bins for each group), and byvals (distinct values in by).

    See Also
    --------
    binsregselect, binsreg, binsqreg, binsglm, binstest.
    
    Example
    ------- 
    >>> x = numpy.random.uniform(size = 500)
    >>> y = numpy.sin(x)+numpy.random.normal(size = 500)
    >>> t  = 1*(numpy.random.uniform(size = 500)>0.5)
    >>> est = binstest(y,x, by=t)
    >>> print(est)
    '''

    # param for internal use
    qrot = 2

    ####################
    ### prepare data ###
    ####################

    # extract x, y, w, by, cluster, weights and subset from data (if data is supplied)
    xname = 'x'
    yname = 'y'
    if by is not None: byname = 'group'
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
        if by is not None:
            if isinstance(by, str):
                if by in data.columns:
                    byname = by
                    by = data[byname]
                else:
                    raise Exception("by column not found in data.")
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
    if by is not None:
        by = np.array(by).reshape(len(by),-1)
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
        if by is not None: by = by[subset]
        if cluster is not None: cluster = cluster[subset]
        if weights is not None: weights = weights[subset]
        
    # remove missing values
    na_ok = complete_cases(x) & complete_cases(y)
    if w is not None:  na_ok = na_ok & complete_cases(w)
    if by is not None: na_ok = na_ok & complete_cases(by)
    if weights is not None: na_ok = na_ok & complete_cases(weights)
    if cluster is not None: na_ok = na_ok & complete_cases(cluster)

    y  = y[na_ok]
    x  = x[na_ok]
    if w is not None: w  = w[na_ok,:]
    if by is not None: by = by[na_ok]
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
            link = 'default'
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
        if vce is None: vce = "HC1"
        vce_select = vce
    elif estmethod=="qreg":
        is_qreg = True
        estmethod_name = "quantile regression"
        if vce is None: vce = "nid"
        if vce=="iid": vce_select = "const"
        else: vce_select = "HC1"
        if quantile is None: quantile = 0.5
    elif estmethod=="glm":
        estmethod_name = "generalized linear model"
        if vce is None: vce = "HC1"
        vce_select = vce
    else:
        raise Exception("estmethod incorrectly specified.")

    ############################ Error Checking
    
    if not isinstance(binspos, str):
        if np.min(binspos)<=xmin or np.max(binspos)>=xmax:
            raise Exception("knots out of allowed range")
    else:
        if binspos!="es" and binspos!="qs":
            raise Exception("binspos incorrectly specified.")
    
    if len(bins)==2:
        if bins[0] < bins[1]:
            raise Exception("p<s not allowed.")

    if len(pwc)==2:
        if pwc[0] < pwc[1]:
            raise Exception("p<s not allowed.")
    
    if len(pwc)>=1:
        if pwc[0]<deriv:
            raise Exception("p for test cannot be smaller than deriv.")

    if (len(pwc)>=1) and (len(bins)>=1):
        if pwc[0]<=bins[0]:
            warnings.warn("p for testing > p for bin selection is suggested.")
 
    ##################################################
    # Prepare options
    bins_p, bins_s = bins
    tsha_p, tsha_s = pwc

    #########################################
    if binsmethod=="dpi": selectmethod = "IMSE direct plug-in" 
    else: selectmethod = "IMSE rule-of-thumb"
    
    nbins_all = None
    if np.isscalar(bynbins): nbins_all = bynbins    # "nbins" is reserved for use within loop
    if bynbins is not None: selectmethod = "User-specified"
    
    ###############################################
    localcheck = massadj = True
    fewmasspoints = False
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
    elif masspoints=="veryfew":
        raise Exception("veryfew not allowed for testing.")

    #############################
    # Extract byvals in by ######
    if by is not None:
        byvals = np.unique(by)
        ngroup = len(byvals)
    else:
        byvals = ["Full Sample"]
        ngroup = 1

    # extract the range for each group
    xminmat = nanmat(ngroup) 
    xmaxmat = nanmat(ngroup)
    if by is None:
        xminmat[0] = xmin
        xmaxmat[0] = xmax
    else:
        for i in range(ngroup):
            x_i = x[by == byvals[i]]
            xminmat[i] = np.min(x_i)
            xmaxmat[i] = np.max(x_i)

    knot = None
    knotlistON = False
    if not isinstance(binspos, str):
        nbins_all = len(binspos)+1
        knot = np.concatenate([[xmin], np.sort(binspos),[xmax]])
        position = "User-specified"
        es = False
        selectmethod = "User-specified"
        knotlistON = True
    else:
        if binspos == "es":
            es = True
            position = "Evenly-spaced"
        else:
            es = False
            position = "Quantile-spaced"

    Ndist = Nclust = np.nan

    ################################################
    ### Bin selection using full sample if needed ##
    selectfullON = False
    if nbins_all is None and samebinsby:
         selectfullON = True
    if selectfullON:
        # effective size
        eN = N = len(x)
        if massadj:
            Ndist = len(np.unique(x))
            eN = min(eN, Ndist)
        if cluster is not None:
            Nclust = len(np.unique(cluster))
            eN = min(eN, Nclust)

        # check if rot can be implemented
        if nbinsrot is None:
            if eN <= dfcheck[0]+bins_p+1+qrot:
                warnings.warn("too small effective sample size for bin selection.")
        if np.isnan(Ndist): Ndist_sel = None
        else: Ndist_sel = Ndist
        if np.isnan(Nclust): Nclust_sel = None
        else: Nclust_sel = Nclust
        binselect = binsregselect(y, x, w, deriv=deriv,
                                    bins=bins, binspos=binspos,
                                    binsmethod=binsmethod, nbinsrot=nbinsrot,
                                    vce=vce_select, cluster=cluster, randcut=randcut,
                                    dfcheck=dfcheck, masspoints=masspoints, weights=weights,
                                    numdist=Ndist_sel, numclust=Nclust_sel)

        if np.isnan(binselect.nbinsrot_regul):
                raise Exception("bin selection fails.")
        if binsmethod == "rot": nbins_all = binselect.nbinsrot_regul
        elif binsmethod == "dpi": nbins_all = binselect.nbinsdpi
        if np.isnan(nbins_all):
            warnings.warn("DPI selection fails. ROT choice used.")
            nbins_all = binselect.nbinsrot_regul

    # Generate knot using the full sample if needed

    if (selectfullON or (nbins_all is not None and samebinsby)) and knot is None:
        knotlistON = True
        if es: knot = genKnot_es(xmin, xmax, nbins_all)
        else: knot = genKnot_qs(x, nbins_all)

    knot_all = None
    if knotlistON: knot_all = knot 

    if selectfullON or (nbins_all is not None and samebinsby) and knot is None:
        knotlistON = True
        if es: knot = genKnot_es(xmin, xmax, nbins_all)
        else: knot = genKnot_qs(x, nbins_all)

    knot_all = None
    if knotlistON: knot_all = knot    # universal knot sequence

    ##########################################
    if simsseed is not None: np.random.seed(simsseed)
    # common grid
    uni_grid = np.linspace(np.max(xminmat), np.min(xmaxmat), num=simsgrid+2)[1:-1]

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

    ##################################################################
    ##################### ENTER the loop #############################
    ##################################################################
    # create empty lists to save results in the loop
    N_by = []
    Ndist_by = []
    Nclust_by = []
    nbins_by = []
    fit_sha = []
    se_sha = []
    nummat = [] 
    denom = []
    fit_sha = []
    se_sha = []
    nummat = []
    denom = []
    
    tstat = nanmat(int(ngroup*(ngroup-1)/2), 3); 
    pval = nanmat(int(ngroup*(ngroup-1)/2), 1)
    counter = 0

    for i in range(ngroup):
        # Take subsample
        if ngroup == 1:
            sub = np.repeat(True,len(x))
        else:
            sub = (by == byvals[i]).reshape(-1)
        y_sub = y[sub]
        x_sub = x[sub]
        if w is not None: w_sub = w[sub,:]
        else: w_sub = None
        if cluster is not None: cluster_sub = cluster[sub]
        else: cluster_sub = None
        if weights is not None: weights_sub = weights[sub]
        else: weights_sub = None

        # Effective size
        xmin = np.min(x_sub)
        xmax = np.max(x_sub)
        eN = N = len(x_sub)
        N_by +=[N]

        Ndist = np.nan
        if massadj:
            Ndist = len(np.unique(x_sub))
            eN = min(eN, Ndist)
        Ndist_by += [Ndist]

        Nclust = np.nan
        if cluster_sub is not None:
            Nclust = len(np.unique(cluster_sub))
            eN = min(eN, Nclust)
        Nclust_by += [Nclust]

        #########################################################
        ############### Bin selection within loop ###############
        nbins = knot = None                  # initialize again
        if nbins_all is not None: nbins = nbins_all
        if bynbins is not None and not np.isscalar(bynbins): nbins = bynbins[i]

        if nbins is None and not knotlistON:
            # check if rot can be implemented
            if nbinsrot is None:
                if eN <= dfcheck[0]+bins_p+1+qrot:
                    raise Exception("too small effective sample size for bin selection.")
  
            if np.isnan(Ndist): Ndist_sel = None
            else: Ndist_sel = Ndist
        
            if np.isnan(Nclust): Nclust_sel = None
            else: Nclust_sel = Nclust
        
            binselect = binsregselect(y_sub, x_sub, w_sub, deriv=deriv,
                                        bins=bins, binspos=binspos,
                                        binsmethod=binsmethod, nbinsrot=nbinsrot,
                                        vce=vce_select, cluster=cluster_sub, randcut=randcut,
                                        dfcheck=dfcheck, masspoints=masspoints, weights=weights_sub,
                                        numdist=Ndist_sel, numclust=Nclust_sel)
            
            if np.isnan(binselect.nbinsrot_regul):
                    raise Exception("bin selection fails.")
            if binsmethod == "rot": nbins = binselect.nbinsrot_regul
            elif binsmethod == "dpi":
                nbins = binselect.nbinsdpi
                if np.isnan(nbins):
                    warnings.warn("DPI selection fails. ROT choice used.")
                    nbins = binselect.nbinsrot_regul
        
        if knotlistON:
            nbins = len(knot_all)-1
            knot  = knot_all

        ###########################################
        # Checking for each case
        if (nbins-1)*(tsha_p-tsha_s+1)+tsha_p+1+dfcheck[1]>=eN:
            warnings.warn("too small effective sample size for testing shape.")
        

        ####################################
        ########### Generate knot ##########
        ####################################
        if knot is None:
            if es: knot = genKnot_es(xmin, xmax, nbins)
            else: knot = genKnot_qs(x_sub, nbins)

        # knot for few mass points
        knot = np.append(knot[0],np.unique(knot[1:]))
        if nbins != len(knot)-1:
            warnings.warn("repeated knots. Some bins dropped.")
            nbins = len(knot)-1

        # NOW, save nbins
        nbins_by += [nbins]

        # check local mass points
        if localcheck:
            uniqmin = binsreg_checklocalmass(x_sub, nbins, es, knot=knot) # mimic STATA
            if uniqmin < tsha_p+1:
                warnings.warn("some bins have too few distinct values of x for testing.")
        
        #######################################
        ###### Estimation #####################
        #######################################

        B = binsreg_spdes(x=x_sub, p=tsha_p, s=tsha_s, knot=knot, deriv=0)
        k = ncol(B)
        P = cbind(B, w_sub)
        model = binsreg_fit(y=y_sub, x=P, weights=weights_sub, family=family_str, is_qreg=is_qreg,
                                 quantile = quantile, cov_type = vce_select, cluster = cluster_sub)            
        beta = model.params[:k]
        basis_sha = binsreg_spdes(x=uni_grid, p=tsha_p, s=tsha_s, knot=knot, deriv=deriv)

        if estmethod=="glm" and not nolink:
            fit_sha_i, se_sha_i = binsreg_pred(X=basis_sha, model=model, type="all",
                                            deriv=deriv, wvec=eval_w, avar=asyvar)
            basis_0 = binsreg_spdes(x=uni_grid, p=tsha_p, s=tsha_s, knot=knot, deriv=0)
            fit_0 = binsreg_pred(basis_0, model, type = "xb", deriv=0, wvec=eval_w)[0]
            pred_sha_0  = linkinv_1(fit_0)

            if asyvar or deriv==0:
                se_sha_i  = pred_sha_0 * se_sha_i
                if deriv == 0: fit_sha_i = linkinv(fit_sha_i)
                if deriv == 1: fit_sha_i = pred_sha_0 * fit_sha_i
            else:
                basis_sha_1 = basis_sha
                if eval_w is not None:
                    basis_sha_0 = np.column_stack((basis_0, np.outer(np.ones(nrow(basis_0)), eval_w)))
                    basis_sha_1 = np.column_stack((basis_sha_1, np.outer(np.ones(nrow(basis_sha_1)), np.zeros(nwvar))))
                basis_all = linkinv_2(fit_0)*fit_sha_i*basis_sha_0 + pred_sha_0*basis_sha_1
                fit_sha_i = pred_sha_0 * fit_sha_i
                se_sha_i  = binsreg_pred(basis_all, model=model, type="se", avar=True)[1]
        else:
            fit_sha_i, se_sha_i  = binsreg_pred(basis_sha, model, type = "all", deriv=deriv,
                                                 wvec=eval_w, avar=asyvar)
        
        fit_sha += [fit_sha_i]
        se_sha += [se_sha_i]

        pos = np.invert(np.isnan(beta))
        k_new = np.sum(pos)
        vcv_sha = model.cov_params()[:k_new,:k_new]
        Sigma_root = lssqrtm(vcv_sha)
        
        nummat += [np.matmul((basis_sha[:,pos]).reshape(-1,k_new), Sigma_root)]
        denom += [np.sqrt(np.sum(np.matmul(basis_sha[:,pos], vcv_sha) * basis_sha[:,pos],1))]

        # second loop over 1:(i-1)
        if i>0:
            for j in range(i):
                if testtype=="left":
                    tstat[counter,:] = np.array((np.max((fit_sha[i]-fit_sha[j]) / np.sqrt(se_sha[i]**2+se_sha[j]**2)), i, j))
                elif testtype=="right":
                    tstat[counter,:] = np.array((np.min((fit_sha[i]-fit_sha[j]) / np.sqrt(se_sha[i]**2+se_sha[j]**2)), i, j))
                else:
                    if not np.isfinite(lp):
                        tstat[counter,:] = np.array((np.max(np.abs((fit_sha[i]-fit_sha[j]) / np.sqrt(se_sha[i]**2+se_sha[j]**2))), i, j))
                    else:
                        tstat[counter,:] = np.array((np.mean(((fit_sha[i]-fit_sha[j]) / np.sqrt(se_sha[i]**2+se_sha[j]**2))**lp)**(1/lp), i, j))               
                pval[counter,0] = binspwc_pval(nummat[i], nummat[j], denom[i], denom[j], nsims, tstat=tstat[counter,0], testtype=testtype, lp=lp)
                counter += 1
        

    ######################################
    ########### Output ###################
    ######################################

    opt = options_pwc(bins_p=bins_p, bins_s=bins_s, deriv=deriv, pwc_p=tsha_p, pwc_s=tsha_s, byname=None,
                        testtype=testtype, binspos=position, binsmethod=selectmethod, N_by=N_by,
                        Ndist_by=Ndist_by, Nclust_by=Nclust_by, nbins_by=nbins_by, byvals=byvals,
                        lp=lp, estmethod=estmethod_name, quantile=quantile, dist=dist, link=link)
    out = pwc_output(tstat, pval, opt)
    return out 
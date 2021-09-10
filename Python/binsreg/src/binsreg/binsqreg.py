#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Sat Sep  4 17:39:50 2021
# @author: Ricardo Masini

import numpy as np
import pandas as pd
from scipy.stats import norm
import warnings
from plotnine import (ggplot,theme_bw,aes,
                      geom_point,geom_line,geom_errorbar,geom_ribbon,
                      scale_color_manual,guide_legend,theme,labs)
from .binsregselect import binsregselect
# from binsregselect import binsregselect
from .funs import *
# from funs import *

def binsqreg(y, x, w=None, data=None, at=None, quantile = 0.5, deriv=0,
            dots=(0,0), dotsgrid=0, dotsgridmean=True, line=None, linegrid=20,
            ci=None, cigrid=0, cigridmean=True, cb=None, cbgrid=20,
            polyreg=None, polyreggrid=20, polyregcigrid=0,
            by=None, bycolors=None, bysymbols=None, bylpatterns=None,
            legendTitle=None, legendoff=False,
            nbins=None, binspos="qs", binsmethod="dpi", nbinsrot=None, samebinsby=False, randcut=None,
            nsims=500, simsgrid=20, simsseed=None,
            vce="HC1",cluster=None, asyvar=False, level=95,
            noplot=False, dfcheck=(20,30), masspoints="on", weights=None, subset=None, plotxrange=None, plotyrange=None):

    '''
    Data-Driven Binscatter Quantile Regression with Robust Inference Procedures and Plots
    
    Description
    ----------- 
    binsqreg implements binscatter quantile regression with robust inference procedures and plots, following the
    results in \href{https://arxiv.org/abs/1902.09608}{Cattaneo, Crump, Farrell and Feng (2021a)}.
    Binscatter provides a flexible way to describe the quantile relationship between two variables, after
    possibly adjusting for other covariates, based on partitioning/binning of the independent variable of interest.
    The main purpose of this function is to generate binned scatter plots with curve estimation with robust pointwise confidence intervals and
    uniform confidence band.  If the binning scheme is not set by the user, the companion function
    binsregselect is used to implement binscatter in a data-driven way. Hypothesis testing about the function of interest
    can be conducted via the companion function binstest.
    
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
    
    at: str
        Value of w at which the estimated function is evaluated.  The default is at="mean", which corresponds to
        the mean of w. Other options are: at="median" for the median of w, at="zero" for a vector of zeros.
        at can also be a vector of the same length as the number of columns of w (if w is a matrix) or a data frame
        containing the same variables as specified in w (when data is specified). Note that when at="mean" or at="median",
        all factor variables (if specified) are excluded from the evaluation (set as zero).
    
    quantile : float
        The quantile to be estimated. A number strictly between 0 and 1.

    deriv : int 
        Derivative order of the regression function for estimation, testing and plotting.
        The default is deriv=0, which corresponds to the function itself.

    dots: tuple 
        dots=(p,s) sets a piecewise polynomial of degree p with s smoothness constraints for
        point estimation and plotting as "dots". The default is dots=(0,0), which corresponds to
        piecewise constant (canonical binscatter)

    dotsgrid : int
        Number of dots within each bin to be plotted. Given the choice, these dots are point estimates
        evaluated over an evenly-spaced grid within each bin. The default is dotsgrid=0, and only
        the point estimates at the mean of x within each bin are presented.

    dotsgridmean : bool
        If true, the dots corresponding to the point estimates evaluated at the mean of x within each bin
        are presented. By default, they are presented, i.e., dotsgridmean=True.

    line : tuple 
        line=(p,s) sets a piecewise polynomial of degree  p with  s smoothness constraints for plotting as a "line".
        By default, the line is not included in the plot unless explicitly specified.  Recommended specification is
        line=(3,3), which adds a cubic B-spline estimate of the regression function of interest to the binned scatter plot.

    linegrid : int
        Number of evaluation points of an evenly-spaced grid within each bin used for evaluation of
        the point estimate set by the line=c(p,s) option. The default is linegrid=20}
        which corresponds to 20 evenly-spaced evaluation points within each bin for fitting/plotting the line.
    
    ci : tuple 
        ci=c(p,s) sets a piecewise polynomial of degree  p with  s smoothness constraints used for
        constructing confidence intervals. By default, the confidence intervals are not included in the plot
        unless explicitly specified.  Recommended specification is ci=(3,3), which adds confidence
        intervals based on cubic B-spline estimate of the regression function of interest to the binned scatter plot.
    
    cigrid : int
        Number of evaluation points of an evenly-spaced grid within each bin used for evaluation of the point
        estimate set by the ci=(p,s) option. The default is cigrid=1, which corresponds to 1
        evenly-spaced evaluation point within each bin for confidence interval construction.
    
    cigridmean : bool
        If true, the confidence intervals corresponding to the point estimates evaluated at the mean of x within each bin
        are presented. The default is cigridmean=True.
    
    cb : tuple 
        cb=(p,s) sets a the piecewise polynomial of degree  p with  s smoothness constraints used for
        constructing the confidence band. By default, the confidence band is not included in the plot unless
        explicitly specified.  Recommended specification is cb=c(3,3), which adds a confidence band
        based on cubic B-spline estimate of the regression function of interest to the binned scatter plot.
    
    cbgrid : int 
        Number of evaluation points of an evenly-spaced grid within each bin used for evaluation of the point
        estimate set by the cb=(p,s) option. The default is cbgrid=20, which corresponds
        to 20 evenly-spaced evaluation points within each bin for confidence interval construction.
    
    polyreg : int
        Degree of a global polynomial regression model for plotting. By default, this fit is not included
        in the plot unless explicitly specified. Recommended specification is polyreg=3, which
        adds a cubic (global) polynomial fit of the regression function of interest to the binned scatter plot.
    
    polyreggrid : int
        number of evaluation points of an evenly-spaced grid within each bin used for evaluation of
        the point estimate set by the polyreg=p option. The default is polyreggrid=20,
        which corresponds to 20 evenly-spaced evaluation points within each bin for confidence
        interval construction.

    polyregcigrid : int
        Number of evaluation points of an evenly-spaced grid within each bin used for constructing
        confidence intervals based on polynomial regression set by the polyreg=p option.
        The default is polyregcigrid=0, which corresponds to not plotting confidence
        intervals for the global polynomial regression approximation.
    
    by : array
        A vector containing the group indicator for subgroup analysis; both numeric and string variables
        are supported. When  by is specified, binsreg implements estimation and inference for each subgroup
        separately, but produces a common binned scatter plot. By default, the binning structure is selected for each
        subgroup separately, but see the option samebinsby below for imposing a common binning structure across subgroups.
    
    bycolors : list or tuple of srt
        An ordered list of colors for plotting each subgroup series defined by the option by.
        Refer to plotnine package documentation.
    
    bysymbols : list of tuple of str
        An ordered list of symbols for plotting each subgroup series defined by the option by.
        Refer to plotnine package documentation.
    
    bylpatterns : list of tuple of str
        An ordered list of line patterns for plotting each subgroup series defined by the option by.
        The options are 'solid','dashed','dashdot', and 'dotted'
    
    legendTitle : str
        Title of legend.
     
    legendoff : bool 
        If true, no legend is added.
    
    nbins : int 
        Number of bins for partitioning/binning of x.  If not specified, the number of bins is
        selected via the companion function binsregselect in a data-driven, optimal way whenever possible.
    
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
        Dummy variable to select the variance-covariance matrix estimator.
        It is currently NOT implemented in statsmodel.api for quantile regression. 
            
    cluster: array 
        Cluster ID. Used for compute cluster-robust standard errors.
    
    asyvar : bool
        If true, the standard error of the nonparametric component is computed and the uncertainty related to control
        variables is omitted. Default is asyvar=FALSE, that is, the uncertainty related to control variables is taken into account.
    
    level : int 
        Nominal confidence level for confidence interval and confidence band estimation. Default is level=95.
    
    noplot : bool
        If true, no plot produced.
    
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
    
    plotxrange : tuple
        plotxrange=(min, max) specifies a range of the x-axis for plotting. Observations outside the range are dropped in the plot.
    
    plotyrange : tuple
        plotyrange=(min, max) specifies a range of the y-axis for plotting. Observations outside the range are dropped in the plot.
    
    Returns
    -------
    bins_plot : A ggplot object for binscatter plot.
    
    data_plot : A list containing data for plotting. Each item is a sublist of data frames for each group. 
                Each sublist may contain the following data frames:
   
                    data_dots : Data for dots. It contains: x, evaluation points; bin, the indicator of bins;
                                isknot, indicator of inner knots; mid, midpoint of each bin; and fit, fitted values.
                    data.line : Data for line. It contains: x, evaluation points;  bin, the indicator of bins;
                                isknot, indicator of inner knots; mid, midpoint of each bin; and fit, fitted values.
                    data.ci : Data for CI. It contains: x, evaluation points;  bin, the indicator of bins;
                              isknot, indicator of inner knots; mid, midpoint of each bin;
                              ci_l and ci_r, left and right boundaries of each confidence intervals.
                    data.cb : Data for CB. It contains: x, evaluation points;  bin, the indicator of bins;
                              isknot, indicator of inner knots; mid, midpoint of each bin;
                              cb_l and cb_r, left and right boundaries of the confidence band.
                    data.poly : Data for polynomial regression. It contains: x, evaluation points; bin, the indicator of bins;
                                isknot, indicator of inner knots; mid, midpoint of each bin; and fit, fitted values.
                    data.polyci : Data for confidence intervals based on polynomial regression. It contains:  x, evaluation points;
                                  bin, the indicator of bins; isknot, indicator of inner knots; mid, midpoint of each bin;
                                  polyci_l and polyci_r, left and right boundaries of each confidence intervals.
                    cval_by : A vector of critical values for constructing confidence band for each group.
                    options : A list containing options passed to the function, as well as N_by (total sample size for each group),
                              Ndist_by (number of distinct values in  x for each group), Nclust_by (number of clusters for each group),
                              and nbins_by (number of bins for each group), and byvals (number of distinct values in by).
    
    See Also
    --------
    binsregselect, binsreg, binsglm, binstest, binspwc.
    
    Example
    ------- 
    >>> x = numpy.random.uniform(size = 500)
    >>> y = numpy.sin(x) + numpy.random.normal(size = 500)
    >>> est = binsqreg(y,x)
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

    ##### Error Checking ##########################

    if deriv < 0:
        raise Exception("derivative incorrectly specified.")

    if dotsgrid<0 or linegrid<0 or cigrid<0 or cbgrid<0 or polyreggrid<0 or polyregcigrid<0:
        raise Exception("# of evaluation points incorrectly specified.")

    if nbins is not None:
        if nbins < 0:
            raise Exception("# of bins incorrectly specified.")
        
    if not isinstance(binspos, str):
        if np.min(binspos)<=xmin or np.max(binspos)>=xmax:
            raise Exception("knots out of allowed range")
    else:
        if binspos!="es" and binspos!="qs":
            raise Exception("binspos incorrectly specified.")
        
    if len(dots)==2:
        if dots[0] < dots[1]:
            raise Exception("p<s not allowed.")

    if line is not None:
        if len(line)==2:
            if line[0]<line[1]:
                raise Exception("p<s not allowed.")

    if ci is not None:
        if len(ci)==2:
            if ci[0]<ci[1]:
                raise Exception("p<s not allowed.")

    if cb is not None:
        if len(cb)==2:
            if cb[0]<cb[1]:
                raise Exception("p<s not allowed.")

    if dots[0] < deriv:
        raise Exception("p<deriv not allowed.")

    if line is not None:
        if line[0] < deriv:
            raise Exception("p<deriv not allowed.")

    if ci is not None:
        if ci[0] < deriv:
            raise Exception("p<deriv not allowed.")

    if cb is not None:
        if cb[0] < deriv:
            raise Exception("p<deriv not allowed.")

    if binsmethod!="dpi" and binsmethod!="rot":
        raise Exception("bin selection method incorrectly specified.")


    if w is not None:
        if isinstance(at,str):
            if at!="mean" and at!="median" and at!="zero":
                raise Exception("at incorrectly specified.")
        else:
            if not np.isscalar(at):
                if len(at)!=nwvar:
                    raise Exception("length of at not equal to # of w variables.")
            
    ##################################################
    # Prepare options
    dots_p, dots_s = dots
    dotsmean = 0
    if dotsgridmean: dotsmean = 1

    if line is None: linegrid = 0
    elif np.isscalar(line): line_p = line_s =  line
    elif len(line)==2: line_p, line_s = line

    cimean = 0
    if cigridmean: cimean = 1
    if ci is None: cigrid  = cimean = 0
    elif np.isscalar(ci): ci_p = ci_s = ci
    elif len(ci)==2: ci_p, ci_s = ci

    if cb is None: cbgrid = 0
    elif np.isscalar(cb): cb_p = cb_s = cb
    elif len(cb)==2: cb_p, cb_s = cb

    if polyreg is None: polyreggrid  = polyregcigrid = 0

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
        fewmasspoints = True

    #############################
    # Extract byvals in by ######
    if by is not None:
        byvals = np.unique(by)
        ngroup = len(byvals)
    else:
        byvals = ["Full Sample"]
        ngroup = 1

    ########################################
    # Default plotting options
    if bycolors is None:
        bycolors = ["navy", "maroon", "forestgreen", "darkorange", "lavenderblush3",
                    "khaki", "sienna", "steelblue", "brown", "gold", "gray45"]
    bycolors = np.resize(bycolors, ngroup)

    if bysymbols is None:
        bysymbols = np.array(['o',  # circle
                    '^',  # triangle up
                    's',  # square
                    'D',  # Diamond
                    'v',  # triangle down
                    '*',  # star
                    'p',  # pentagon
                    '8',  # octagon
                    '<',  # triangle left
                    'h',  # hexagon1    
                    '>',  # triangle right
                    'H',  # hexagon1
                    'd'])  # thin diamond
    bysymbols = np.resize(bysymbols, ngroup)

    if bylpatterns is None: bylpatterns = np.array(['solid'])
    bylpatterns = np.resize(bylpatterns, ngroup)

    # Legend
    if ngroup==1: legendoff = True    # turn off legend when only one group is available
    else:
        if legendTitle is None: legendTitle = byname + " ="

    #########################################
    if binsmethod=="dpi": selectmethod = "IMSE direct plug-in" 
    else: selectmethod = "IMSE rule-of-thumb"

    nbins_all = nbins         # "nbins" is reserved for use within loop
    if nbins is not None: selectmethod = "User-specified"

    knot = None
    knotlistON = False
    if not isinstance(binspos, str):
        nbins = len(binspos)+1
        knot = np.concatenate([[xmin], np.sort(x),[xmax]])
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
    fullfewobs = fewobs = selectfullON = False
    if fewmasspoints: fullfewobs = fewobs = True
    if (not fullfewobs and nbins is None and by is None) or (by is not None and samebinsby): 
        selectfullON = True
    if selectfullON:
        eN = N = len(x) # effective size
        if massadj:
            Ndist = len(np.unique(x))
            eN = min(eN, Ndist)
        if cluster is not None:
            Nclust = len(np.unique(cluster))
            eN = min(eN, Nclust)

        # check if rot can be implemented
        if nbinsrot is None:
            if eN <= dfcheck[0]+dots_p+1+qrot:
                warnings.warn("too small effective sample size for bin selection. # of mass of points or clusters used and by option ignored.")
                fewobs = fullfewobs = True
                byvals = "Full Sample"
                es = False
                position = "Quantile-spaced"
        if not fewobs:
            if np.isnan(Ndist): Ndist_sel = None
            else: Ndist_sel = Ndist
            if np.isnan(Nclust): Nclust_sel = None
            else: Nclust_sel = Nclust
            binselect = binsregselect(y, x, w, deriv=deriv,
                                        bins=dots, binspos=binspos,
                                        binsmethod=binsmethod, nbinsrot=nbinsrot,
                                        vce=vce, cluster=cluster, randcut=randcut,
                                        dfcheck=dfcheck, masspoints=masspoints, weights=weights,
                                        numdist=Ndist_sel, numclust=Nclust_sel)
            if np.isnan(binselect.nbinsrot_regul):
                raise Exception("bin selection fails.")
            if binsmethod == "rot": nbins = binselect.nbinsrot_regul
            elif binsmethod == "dpi": nbins = binselect.nbinsdpi
            if np.isnan(nbins):
                warnings.warn("DPI selection fails. ROT choice used.")
                nbins = binselect.nbinsrot_regul

    # Generate knot using the full sample if needed

    if (selectfullON or (nbins is not None and samebinsby)) and not fullfewobs and knot is None:
        knotlistON = True
        if es: knot = genKnot_es(xmin, xmax, nbins)
        else: knot = genKnot_qs(x, nbins)

    knot_all = None
    if knotlistON: knot_all = knot    # universal knot sequence

    ##########################################
    if simsseed is not None: np.random.seed(simsseed)
    alpha = (100-(100 - level)/2)/100
   

    ##################################################################
    ##################### ENTER the loop #############################
    ##################################################################
    # save results
    N_by = []
    Ndist_by = []
    Nclust_by = []
    nbins_by = []
    cval_by = []
    data_plot = []   

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

        #######################################################
        ############### Bin selection if needed ###############
        nbins = knot = None                       # initialize again
        if (nbins_all is None) and (not knotlistON) and (not fullfewobs):
            # check if rot can be implemented
            if (nbinsrot is None) and (eN <= dfcheck[0]+dots_p+1+qrot):
                warnings.warn("too small effective sample size for bin selection. # of mass points or clusters used.")
                fewobs = True
                nbins = eN
                es = False
            if not fewobs:
                if np.isnan(Ndist): Ndist_sel = None
                else: Ndist_sel = Ndist
                if np.isnan(Nclust): Nclust_sel = None
                else: Nclust_sel = Nclust

                binselect = binsregselect(y_sub, x_sub, w_sub, deriv=deriv,
                                            bins=dots, binspos=binspos,
                                            binsmethod=binsmethod, nbinsrot=nbinsrot,
                                            vce=vce, cluster=cluster_sub, randcut=randcut,
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

        if nbins_all is not None: nbins = nbins_all
        if knotlistON:
            nbins = len(knot_all)-1
            knot  = knot_all
        if fullfewobs:
            fewobs = True
            nbins = eN

        ###########################################
        # Checking Error  for each case
        dots_fewobs = line_fewobs = ci_fewobs = cb_fewobs = polyreg_fewobs = False
        if not fewobs:
            if (nbins-1)*(dots_p-dots_s+1)+dots_p+1+dfcheck[1]>=eN: 
                fewobs = True
                nbins = eN
                es = False
                warnings.warn("too small effective sample size for dots. # of mass points or clusters used.")
            if line is not None:
                if ((nbins-1)*(line_p-line_s+1)+line_p+1+dfcheck[1]>=eN):
                    line_fewobs = True
                    warnings.warn("too small effective sample size for line.")
            if ci is not None:
                if (nbins-1)*(ci_p-ci_s+1)+ci_p+1+dfcheck[1]>=eN:
                    ci_fewobs = True
                    warnings.warn("too small effective sample size for ci.")
            if cb is not None:
                if (nbins-1)*(cb_p-cb_s+1)+cb_p+1+dfcheck[1]>=eN:
                    cb_fewobs = True
                    warnings.warn("too small effective sample size for line.")

        if polyreg is not None:
            if polyreg+1>eN:
                polyreg_fewobs = True
                warnings.warn("too small effective sample size for polynomial fit.")


        ####################################
        ########### Generate knot ##########
        ####################################
        fewmass = False
        if fewobs and (not np.isnan(Ndist)) and eN==Ndist: fewmass = True
        if fewmasspoints: fewmass = True

        if knot is None:
            if not fewmass:
                if es: knot = genKnot_es(xmin, xmax, nbins)
                else: knot = genKnot_qs(x_sub, nbins)
        else:
            if fewobs and (not np.isnan(Ndist)) and (eN!=Ndist):
                knot = genKnot_qs(x_sub, nbins)

        # knot for few mass points
        if fewmass:
            knot = np.unique(x_sub)
            if fewmasspoints:
                eN = nbins = Ndist = len(knot)
        else:
            knot = np.append(knot[0],np.unique(knot[1:]))
            if nbins != len(knot)-1:
                warnings.warn("repeated knots. Some bins dropped.")
                nbins = len(knot)-1

        # check local mass points
        if not fewobs and localcheck:
            uniqmin = binsreg_checklocalmass(x_sub, nbins, es, knot=knot) # mimic STATA
            if dots is not None: 
                if uniqmin < dots_p+1:
                    dots_fewobs = True
                    warnings.warn("some bins have too few distinct values of x for dots.")
            if line is not None:
                if uniqmin < line_p+1:
                    line_fewobs = True
                    warnings.warn("some bins have too few distinct values of x for line.")
            if ci is not None:
                if uniqmin < ci_p+1:
                    ci_fewobs = True
                    warnings.warn("some bins have too few distinct values of x for CI.")
            if cb is not None:
                if uniqmin < cb_p+1:
                    cb_fewobs = True
                    warnings.warn("some bins have too few distinct values of x for CB.")

        # NOW, save nbins
        nbins_by += [nbins]
    

        #######################################
        ###### Prepare Plots ##################
        #######################################
        data_by = data_str()    # initialize data list
        # data_dots = data.line = data.ci= data.cb = data.poly = data.polyci = None

        # adjust w variables
        if w_sub is not None: 
            if isinstance(at, str):
                if at=="mean":
                    eval_w = colWeightedMeans(x=w_sub, w=weights_sub)
                elif at=="median":
                    eval_w = colWeightedMedians(x=w_sub, w=weights_sub)
                elif at=="zero": eval_w = np.zeros(nwvar)
            else: eval_w = np.array(at).reshape(-1,1)
        else: eval_w = None

        ##################################
        # Dots and CIs for Small eN case
        ##################################

        if dotsmean+dotsgrid != 0 and fewobs:
            warnings.warn("dots=c(0,0) used.")
            k = nbins
            if not np.isnan(Ndist):
                if eN == Ndist:
                    dots_x = knot
                    # RENEW knot, each value forms a category
                    xcat_few  = FindInterval(x_sub,knot)
                else:
                    dots_x = (knot[1:]+knot[:-1])/2
                    xcat_few = FindInterval(x_sub,knot)
            else:
                dots_x = (knot[1:]+knot[:-1])/2
                xcat_few  = FindInterval(x_sub,knot)

            design = binsreg_spdes(x=x_sub, p=0, s=0, deriv=0, knot=xcat_few)
            if w_sub is not None: design = np.column_stack((design,w_sub))
            model = binsreg_fit(y=y, x=design, is_qreg = True, quantile = quantile,
                                weights = weights_sub, cov_type = vce, cluster = cluster_sub)
            beta = model.params[:k].values
            beta[np.isnan(beta)] = 0
            vcv = model.cov_params()

            dots_fit = beta.copy()
            dots_fit_0 = beta.copy()
            if eval_w is not None:
                coeff_w = model.params[k:]
                coeff_w[np.isnan(coeff_w)] = 0
                dots_fit += dots_fit + np.sum(eval_w * coeff_w)
            
            data_dots = pd.DataFrame({'group': str(byvals[i]),
                                      'x': dots_x,
                                      'fit': dots_fit})
            data_by.dots = data_dots
            if cigrid+cimean!=0:
                warnings.warn("ci=(0,0) used.")
                basis_all = np.column_stack((np.identity(len(dots_x)), np.outer(np.ones(dots_x), eval_w)))
                dots_se  = np.sqrt(np.sum(np.matmul(basis_all,vcv) * basis_all,0))
                ci_arm = norm.ppf(alpha)*dots_se
                ci_l = dots_fit - ci_arm
                ci_r = dots_fit + ci_arm
                data_ci = pd.DataFrame({'group':str(byvals[i]),
                                      'x':dots_x,
                                      'ci_l':ci_l,
                                      'ci_r':ci_r})

                data_by.ci = data_ci
            
        ##########################################
        ########## Usual Case ####################
        ##########################################

        dotsON = lineON = ciON = cbON = polyON = False
        xmean = None

        ################ Dots ####################
        if dotsmean+dotsgrid !=0 and not dots_fewobs and not fewobs:
            dotsON = True
        if dotsON:
            dots_x = []
            dots_bin = []
            dots_isknot = []
            dots_mid = []
            if dotsmean!=0:
                xcat  = FindInterval(x_sub, knot)
                xmean = [np.mean(x_sub[xcat==i]) for i in range(len(knot)-1)]
                dots_x += xmean
                dots_bin += list(range(nbins))
                dots_isknot += [0]*nbins
                dots_mid += [0]*nbins
            if dotsgrid!=0:
                grid = binsreg_grid(knot, dotsgrid, addmore=True)
                dots_x = np.concatenate((dots_x, grid.eval))
                dots_bin = np.concatenate((dots_bin, grid.bin))
                dots_isknot = np.concatenate((dots_isknot, grid.isknot))
                dots_mid = np.concatenate((dots_mid, grid.mid))
            B = binsreg_spdes(x=x_sub, p=dots_p, s=dots_s, deriv=0, knot=knot)
            if w_sub is not None: design = np.column_stack((B, w_sub))
            else: design = B
            model_dots = binsreg_fit(y=y_sub, x=design, is_qreg=True, quantile = quantile, weights = weights_sub)
            check_drop(model_dots.params, ncol(B))

            basis = binsreg_spdes(x=dots_x, p=dots_p, s=dots_s, knot=knot, deriv=deriv)
            dots_fit, dots_se  = binsreg_pred(basis, model_dots, type = "xb", deriv=deriv, wvec=eval_w)
            
            dots_fit[dots_isknot==1] = np.nan
            data_dots = pd.DataFrame({'group':str(byvals[i]),
                                      'x':dots_x,
                                      'bin':dots_bin,
                                      'isknot':dots_isknot,
                                      'mid':dots_mid,
                                      'fit':dots_fit})
            data_by.dots = data_dots    

        ################ Line ####################
        if linegrid !=0 and not line_fewobs and not fewobs:
            lineON = True
        if lineON:
            grid = binsreg_grid(knot, linegrid, addmore=True)
            line_x = grid.eval
            line_bin = grid.bin
            line_isknot = grid.isknot
            line_mid = grid.mid

            line_reg_ON = True
            if dotsON:
                if line_p==dots_p & line_s==dots_s:
                    model_line = model_dots
                    line_reg_ON = False
            if line_reg_ON:
                B  = binsreg_spdes(x=x_sub, p=line_p, s=line_s, deriv=0, knot=knot)
                if w_sub is not None: design = np.column_stack((B, w_sub))
                else: design = B
                model_line = binsreg_fit(y=y_sub, x=design, is_qreg=True, quantile=quantile, weights=weights_sub)
                check_drop(model_line.params, ncol(B))

            basis = binsreg_spdes(x=line_x, p=line_p, s=line_s, knot=knot, deriv=deriv)
            line_fit, line_se = binsreg_pred(basis, model_line, type = "xb", deriv=deriv, wvec=eval_w)

            if line_s == 0 or line_s-deriv <= 0: line_fit[line_isknot==1] = np.nan
            data_line = pd.DataFrame({'group': str(byvals[i]),
                                      'x': line_x,
                                      'bin': line_bin,
                                      'isknot': line_isknot,
                                      'mid': line_mid, 
                                      'fit': line_fit})
            data_by.line = data_line

        ############### Poly fit #########################
        if polyreggrid!=0 and not polyreg_fewobs:
            polyON = True
        if polyON:
            grid = binsreg_grid(knot, polyreggrid, addmore=True)
            poly_x = grid.eval
            poly_bin = grid.bin
            poly_isknot = grid.isknot
            poly_mid = grid.mid

            # Run a poly reg
            x_p = nanmat(N, polyreg+1)
            for j in range(polyreg+1):  x_p[:,j] = (x_sub**j).reshape(-1)
            if w_sub is not None: design = np.column_stack((x_p, w_sub))
            else: design = x_p
            model_poly = binsreg_fit(y=y_sub, x=design, is_qreg=True, quantile=quantile, weights=weights_sub)
            beta_poly = model_poly.params
            beta_poly[np.isnan(beta_poly)] = 0
            poly_fit  = 0
            for j in range(deriv,polyreg+1):
                poly_fit +=  poly_x**(j-deriv)*beta_poly[j]*factorial(j)/factorial(j-deriv)
            if eval_w is not None and deriv==0:
                poly_fit += np.sum(beta_poly[polyreg+1:]*eval_w)
    
            data_poly = pd.DataFrame({'group': str(byvals[i]),
                                      'x': poly_x,
                                      'bin': poly_bin,
                                      'isknot': poly_isknot,
                                      'mid': poly_mid, 
                                      'fit': poly_fit})
            data_by.poly = data_poly

            # add CI?
            if polyregcigrid!=0:
                grid = binsreg_grid(knot, polyregcigrid, addmore=True)
                polyci_x = grid.eval
                polyci_bin = grid.bin
                polyci_isknot = grid.isknot
                polyci_mid = grid.mid

                npolyci_x = len(polyci_x)
                basis_polyci = nanmat(npolyci_x, polyreg+1)
                for j in range(polyreg+1):
                    if j>=deriv:
                        basis_polyci[:,j] = polyci_x^(j-deriv)*factorial(j)/factorial(j-deriv)
                    else:
                        basis_polyci[:,j] = np.zeros(npolyci_x)
                
                if eval_w is not None:
                    if deriv==0:
                        basis_polyci = np.column_stack((basis_polyci, np.outer(np.ones(basis_polyci.shape[0]), eval_w)))
                    else:          
                        basis_polyci = np.column_stack((basis_polyci, np.outer(np.ones(basis_polyci.shape[0]), np.zeros(nwvar))))

                polyci_fit, polyci_se = binsreg_pred(basis_polyci, model=model_poly, type="all", avar=True)

                polyci_arm = norm.ppf(alpha)*polyci_se
                polyci_l = polyci_fit - polyci_arm
                polyci_r = polyci_fit + polyci_arm

                data_polyci = pd.DataFrame({'group': str(byvals[i]),
                                      'x': polyci_x,
                                      'bin': polyci_bin,
                                      'isknot': polyci_isknot,
                                      'mid': polyci_mid, 
                                      'polyci_l': polyci_l,
                                      'polyci_r': polyci_r})
                data_by.polyci = data_polyci
            
        ################ CI ####################
        if cimean+cigrid !=0 and not ci_fewobs and not fewobs:
            ciON = True
        if ciON:
            ci_x = []
            ci_bin = []
            ci_isknot = [] 
            ci_mid = []
            if cimean!=0:
                if xmean is not None:
                    ci_x += xmean
                else:
                    xcat  = FindInterval(x_sub, knot)
                    ci_x += [np.mean(x_sub[xcat==i]) for i in range(len(knot)-1)]
                ci_bin += list(range(nbins))
                ci_isknot += [0]*nbins
                ci_mid += [0]*nbins
            if cigrid!=0:
                grid = binsreg_grid(knot, cigrid, addmore=True)
                ci_x = np.concatenate((ci_x, grid.eval))
                ci_bin = np.concatenate((ci_bin, grid.bin))
                ci_isknot = np.concatenate((ci_isknot, grid.isknot))
                ci_mid = np.concatenate((ci_mid, grid.mid))
            
            ci_reg_ON = True
            if lineON: 
                if ci_p==line_p and ci_s==line_s:
                    model_ci = model_line 
                    ci_reg_ON = False
            if ci_reg_ON:
                if dotsON:
                    if ci_p==dots_p and ci_s==dots_s:
                        model_ci = model_dots
                        ci_reg_ON = False
            if ci_reg_ON:
                B  = binsreg_spdes(x=x_sub, p=ci_p, s=ci_s, deriv=0, knot=knot)
                if w_sub is not None: design = np.column_stack((B, w_sub))
                else: design = B
                model_ci = binsreg_fit(y=y_sub, x=design, is_qreg=True, quantile=quantile, weights=weights_sub)
                check_drop(model_ci.params, ncol(B))
            
            basis = binsreg_spdes(x=ci_x, p=ci_p, s=ci_s, knot=knot, deriv=deriv)
            ci_fit, ci_se = binsreg_pred(X=basis, model=model_ci, type="all", deriv=deriv, wvec=eval_w, avar=asyvar)
            
            ci_arm = norm.ppf(alpha)*ci_se
            ci_l = ci_fit - ci_arm
            ci_r = ci_fit + ci_arm
            ci_l[ci_isknot==1] = np.nan
            ci_r[ci_isknot==1] = np.nan
            data_ci = pd.DataFrame({'group': str(byvals[i]),
                                      'x': ci_x,
                                      'bin': ci_bin,
                                      'isknot': ci_isknot,
                                      'mid': ci_mid, 
                                      'ci_l': ci_l,
                                      'ci_r': ci_r})
            data_by.ci = data_ci

        ################ CB ###############################
        cval = np.nan
        if cbgrid !=0 and not cb_fewobs and not fewobs:
            cbON = True
        if cbON:
            grid = binsreg_grid(knot, cbgrid, addmore=True)
            cb_x = grid.eval
            cb_bin = grid.bin
            cb_isknot = grid.isknot
            cb_mid = grid.mid
            cb_reg_ON = True
            if ciON:
                if cb_p==ci_p & cb_s==ci_s:
                    model_cb = model_ci
                    cb_reg_ON = False
            if cb_reg_ON:
                if lineON:
                    if cb_p==line_p & cb_s==line_s:
                        model_cb = model_line
                        cb_reg_ON = False
            if cb_reg_ON:
                if dotsON:
                    if cb_p==dots_p & cb_s==dots_s:
                        model_cb = model_dots
                        cb_reg_ON = False
            if cb_reg_ON:
                B = binsreg_spdes(x=x_sub, p=cb_p, s=cb_s, deriv=0, knot=knot)
                if w_sub is not None: design = np.column_stack((B, w_sub))
                else: design = B
                model_cb = binsreg_fit(y=y_sub, x=design, is_qreg=True, quantile=quantile, weights=weights_sub)
                check_drop(model_cb.params, ncol(B))
            
            basis = binsreg_spdes(x=cb_x, p=cb_p, s=cb_s, knot=knot, deriv=deriv)
            pos = np.invert(np.isnan(model_cb.params[:ncol(basis)]))
            k_new = np.sum(pos)
            cb_fit, cb_se = binsreg_pred(X=basis, model=model_cb, type="all", deriv=deriv, wvec=eval_w, avar=asyvar)

            ### Compute cval ####
            x_grid = binsreg_grid(knot, simsgrid).eval
            basis_sim = binsreg_spdes(x=x_grid, p=cb_p, s=cb_s, knot=knot, deriv=deriv)
            sim_fit, sim_se = binsreg_pred(X=basis_sim, model=model_cb, type="all", avar=True)
            vcv = model_cb.cov_params()[:k_new,:k_new]
            Sigma_root = lssqrtm(vcv)
            num = np.matmul(basis_sim[:,pos], Sigma_root)
            
            
            pval, cval = binsreg_pval(num, sim_se, rep=nsims, tstat=None, side="two", alpha=level)
            cb_arm = cval*cb_se
            cb_l = cb_fit - cb_arm
            cb_r = cb_fit + cb_arm
            if (cb_s == 0 | cb_s - deriv <=0):
                cb_l[cb_isknot==1] = np.nan
                cb_r[cb_isknot==1] = np.nan
            data_cb = pd.DataFrame({'group': str(byvals[i]),
                                      'x': cb_x,
                                      'bin': cb_bin,
                                      'isknot': cb_isknot,
                                      'mid': cb_mid, 
                                      'cb_l': cb_l,
                                      'cb_r': cb_r})
            data_by.cb = data_cb
        cval_by += [cval]

        # Save all data for each group
        data_plot += [data_by]
        # names(data_plot)[i] = paste("Group", byvals[i], sep=" ")  # Uncommented in R
    
    #############
    # END Loop
    #############

    ########################################
    ############# Plotting ? ################
    binsplot = None
    if not noplot:
        binsplot = ggplot() + theme_bw()
        x = fit = ci_l = ci_r = cb_l = cb_r = polyci_l = polyci_r = group = None
        xr_min = yr_min = -np.inf
        xr_max = yr_max = np.inf
        if plotxrange is not None:
            xr_min = plotxrange[0]
            if len(plotxrange)==2: xr_max = plotxrange[1]
        if plotyrange is not None:
            yr_min = plotyrange[0]
            if len(plotyrange)==2: yr_max = plotyrange[1]
        
        for i in range(ngroup):
            data_by = data_plot[i]
            if data_by.dots is not None:
                index = (complete_cases(data_by.dots[['x','fit']]) & (data_by.dots["x"]>=xr_min) & (data_by.dots["x"]<=xr_max) &
                        (data_by.dots["fit"]>=yr_min) & (data_by.dots["fit"]<=yr_max))
                if not legendoff:
                    binsplot = binsplot + geom_point(aes(x='x', y='fit', colour='group'), data_by.dots.loc[index], shape=bysymbols[i], size=2)
                else:
                    binsplot = binsplot + geom_point(aes(x='x', y='fit'), data_by.dots.loc[index], colour=bycolors[i], shape=bysymbols[i], size=2)
            
            if data_by.line is not None:
                index = ((data_by.line["x"]>=xr_min) & (data_by.line["x"]<=xr_max) &
                        (((data_by.line["fit"]>=yr_min) & (data_by.line["fit"]<=yr_max)) | np.isnan(data_by.line["fit"])))
                if not legendoff:
                    binsplot = binsplot + geom_line(aes(x='x', y='fit', colour='group'), data_by.line.loc[index], linetype=bylpatterns[i])
                else:
                    binsplot = binsplot + geom_line(aes(x='x', y='fit'), data_by.line.loc[index], color=bycolors[i], linetype=bylpatterns[i])
            if data_by.poly is not None:
                index = ((data_by.poly["x"]>=xr_min) & (data_by.poly["x"]<=xr_max) &
                         (data_by.poly["fit"]>=yr_min) & (data_by.poly["fit"]<=yr_max))
                binsplot = binsplot + geom_line(aes(x='x', y='fit'), data_by.poly.loc[index], color=bycolors[i], inetype=bylpatterns[i])         
            
            if data_by.polyci is not None:
                index = ((data_by.polyci["x"]>=xr_min) & (data_by.polyci["x"]<=xr_max) &
                         (data_by.polyci["polyci_l"]>=yr_min) & (data_by.polyci["polyci_r"]<=yr_max))
                binsplot = binsplot + geom_errorbar(aes(x='x', ymin='polyci_l', ymax='polyci_r'), data_by.polyci.loc[index], alpha=1, color=bycolors[i], linetype=bylpatterns[i])
            
            if data_by.ci is not None:
                index = (complete_cases(data_by.ci[["x", "ci_l", "ci_r"]]) & (data_by.ci["x"]>=xr_min) &
                         (data_by.ci["x"]<=xr_max) & (data_by.ci["ci_l"]>=yr_min) & (data_by.ci["ci_r"]<=yr_max))
                binsplot = binsplot + geom_errorbar(aes(x='x', ymin='ci_l', ymax='ci_r'), data_by.ci.loc[index],color=bycolors[i], alpha =1, width = 0.05, linetype=bylpatterns[i])
            
            if data_by.cb is not None:
                index = ((data_by.cb["x"]>=xr_min) & (data_by.cb["x"]<=xr_max) &
                        (((data_by.cb["cb_l"]>=yr_min) & (data_by.cb["cb_r"]<=yr_max)) | np.isnan(data_by.cb["cb_l"])))
                binsplot = binsplot + geom_ribbon(aes(x='x', ymin='cb_l', ymax='cb_r'), data_by.cb.loc[index],alpha=0.2, fill=bycolors[i])
        
        # Add legend ?
        if not legendoff:
            binsplot = binsplot + scale_color_manual(name=legendTitle, values = bycolors, 
                                                    guide = guide_legend(override_aes ={"linetype":bylpatterns,"shape":bysymbols}))
        else:
            binsplot = binsplot + theme(legend_position="none")
        binsplot = binsplot + labs(x=xname, y=yname)

    ######################################
    ########### Output ###################
    ######################################

    opt = options_qreg(dots=dots, line=line, ci=ci, cb=cb, polyreg=polyreg, deriv=deriv, 
                        quantile=quantile,binspos=position, binsmethod=selectmethod,
                        N_by=N_by, Ndist_by=Ndist_by, Nclust_by=Nclust_by,
                        nbins_by=nbins_by, byvals=byvals)

    out = binsqreg_output(bins_plot=binsplot, data_plot=data_plot, cval_by=cval_by, options=opt)
    return out

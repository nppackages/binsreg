#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Thu Mar 16 17:39:50 2023
# @author: Ricardo Masini

import numpy as np
import pandas as pd
from scipy.stats import norm
import warnings
from plotnine import (ggplot,theme_bw,aes,
                      geom_point,geom_line,geom_errorbar,geom_ribbon,
                      scale_color_manual,guide_legend,theme,labs,xlim)
from binsreg.binsregselect import binsregselect
from binsreg.funs import *

def binsqreg(y, x, w=None, data=None, at=None, quantile=0.5, deriv=0,
            dots=None, dotsgrid=0, dotsgridmean=True, line=None, linegrid=20,
            ci=None, cigrid=0, cigridmean=True, cb=None, cbgrid=20,
            polyreg=None, polyreggrid=20, polyregcigrid=0,
            by=None, bycolors=None, bysymbols=None, bylpatterns=None,
            legendTitle=None, legendoff=False,
            nbins=None, binspos="qs", binsmethod="dpi", nbinsrot=None, samebinsby=False, randcut=None,
            pselect=None, sselect=None,
            nsims=500, simsgrid=20, simsseed=None,
            vce="nid",cluster=None, asyvar=False, level=95,
            noplot=False, dfcheck=(20,30), masspoints="on", weights=None, subset=None, 
            plotxrange=None, plotyrange=None, **optimize):

    '''
    Data-Driven Binscatter Quantile Regression with Robust Inference Procedures and Plots.

    Description
    -----------
    binsreg implements binscatter quantile regression with robust inference procedures and plots, following the
    results in \href{https://arxiv.org/abs/1902.09608}{Cattaneo, Crump, Farrell and Feng (2022a)}.
    Binscatter provides a flexible way to describe the quantile relationship between two variables, after
    possibly adjusting for other covariates, based on partitioning/binning of the independent variable of interest.
    The main purpose of this function is to generate binned scatter plots with curve estimation with robust pointwise confidence intervals and
    uniform confidence band.  If the binning scheme is not set by the user, the companion function
    binsregselect is used to implement binscatter in a data-driven way.
    Hypothesis testing about the regression function of interest can be conducted via the companion binstest.

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

    dots: tuple or bool
        If dots=(p,s), a piecewise polynomial of degree p with s smoothness constraints is used for
        point estimation and plotting as "dots". The default is dots=(0,0), which corresponds to
        piecewise constant (canonical binscatter). If dots=True, the default dots=(0,0) is used unless
        the degree p and smoothness s selection is requested via the option pselect (see more details in the explanation of pselect).
        If dots=False is specified, the dots are not included in the plot.

    dotsgrid : int
        Number of dots within each bin to be plotted. Given the choice, these dots are point estimates
        evaluated over an evenly-spaced grid within each bin. The default is dotsgrid=0, and only
        the point estimates at the mean of x within each bin are presented.

    dotsgridmean : bool
        If true, the dots corresponding to the point estimates evaluated at the mean of x within each bin
        are presented. By default, they are presented, i.e., dotsgridmean=True.

    line : tuple  or bool
        If line=(p,s), a piecewise polynomial of degree  p with s smoothness constraints is used for plotting as a "line".
        If line=True is specified, line=(0,0) is used unless the degree p and smoothness selection is requested via the option
        pselect (see more details in the explanation of pselect and sselect). If line=False or line=None (default) is specified, the line is not included in the plot.
        The default is line=None.

    linegrid : int
        Number of evaluation points of an evenly-spaced grid within each bin used for evaluation of
        the point estimate set by the line=c(p,s) option. The default is linegrid=20}
        which corresponds to 20 evenly-spaced evaluation points within each bin for fitting/plotting the line.
    
    ci : tuple or bool
        If ci=(p,s), a piecewise polynomial of degree  p with  s smoothness constraints is used for
        constructing confidence intervals. If ci=True is specified, ci=(1,1) is used unless the degree p and smoothness s
        selection is requested via the option pselect (see more details in the explanation of pselect).
        If ci=False or ci=None (default) is specified, the confidence intervals are not included in the plot.
    
    cigrid : int
        Number of evaluation points of an evenly-spaced grid within each bin used for evaluation of the point
        estimate set by the ci=(p,s) option. The default is cigrid=1, which corresponds to 1
        evenly-spaced evaluation point within each bin for confidence interval construction.
    
    cigridmean : bool
        If true, the confidence intervals corresponding to the point estimates evaluated at the mean of x within each bin
        are presented. The default is cigridmean=True.
    
    cb : tuple or bool
        If cb=(p,s), a the piecewise polynomial of degree  p with  s smoothness constraints is used for
        constructing the confidence band. If the option cb=True is specified, cb=(1,1) is used unless the degree p and smoothness s          
        selection is requested via the option pselect (see more details in the explanation of pselect).
        If cb=False or cb=None (default) is specified, the confidence band is not included in the plot.
    
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
    
    nbins : int or bool or array
        Number of bins for partitioning/binning of x.  If nbins=True or nbins=None (default) is specified, the number
        of bins is selected via the companion command binsregselect in a data-driven, optimal way whenever possible.
        If a vector with more than one number is specified, the number of bins is selected within this vector via the companion command binsregselect.
    
    binspos : array
        Position of binning knots. The default is binspos="qs", which corresponds to quantile-spaced
        binning (canonical binscatter). The other options are "es" for evenly-spaced binning, or
        a vector for manual specification of the positions of inner knots (which must be within the range of x).
    
    binsmethod : str
        Method for data-driven selection of the number of bins. The default is binsmethod="dpi",
        which corresponds to the IMSE-optimal direct plug-in rule.  The other option is: "rot"
        for rule of thumb implementation.

    nbinsrot : int
        Initial number of bins value used to construct the DPI number of bins selector.
        If not specified, the data-driven ROT selector is used instead.

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
    
    samebinsby : bool
        If true, a common partitioning/binning structure across all subgroups specified by the option by is forced.
        The knots positions are selected according to the option binspos and using the full sample. If \code{nbins}
        is not specified, then the number of bins is selected via the companion command binsregselect and
        using the full sample.
        
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
    
    level : int 
        Nominal confidence level for confidence interval and confidence band estimation. Default is level=95.
    
    noplot : bool
        If true, no plot is produced.
    
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
    
    plotxrange : tuple
        plotxrange=(min, max) specifies a range of the x-axis for plotting. Observations outside the range are dropped in the plot.
    
    plotyrange : tuple
        plotyrange=(min, max) specifies a range of the y-axis for plotting. Observations outside the range are dropped in the plot.
    
    **optimize : 
        Optional arguments to the QuantReg statsmodels.api optimizer. 
        For futher details \href{https://www.statsmodels.org/dev/generated/statsmodels.regression.quantile_regression.QuantReg.fit.html}{sm.QuantReg.fit()}.
    
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
                    imse_v_rot : Variance constant in IMSE, ROT selection.
                    imse_b_rot : Bias constant in IMSE, ROT selection.
                    imse_v_dpi : Variance constant in IMSE, DPI selection.
                    imse_b_dpi : Bias constant in IMSE, DPI selection.
                    cval_by : A vector of critical values for constructing confidence band for each group.
                    options : A list containing options passed to the function, as well as N_by (total sample size for each group),
                              Ndist_by (number of distinct values in  x for each group), Nclust_by (number of clusters for each group),
                              and nbins_by (number of bins for each group), and byvals (number of distinct values in by).
                              The degree and smoothness of polynomials for dots, line, confidence intervals and confidence band for each group are saved
                              in dots, line, ci, and cb.
    
    See Also
    --------
    binsregselect, binsreg, binsglm, binstest, binspwc.
    
    Example
    ------- 
    >>> x = numpy.random.uniform(size = 500)
    >>> y = numpy.sin(x) + numpy.random.normal(size = 500)
    >>> out = binsqreg(y,x)
    >>> print(out)
    >>> out.summary()
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
    # define the support of plot
    if plotxrange is not None and len(plotxrange)==2:
        xsc_min, xsc_max = plotxrange
    else:
        xsc_min = xmin
        xsc_max = xmax

    # evaluation point of w
    if at is None: at = "mean"

    if (vce=="iid"): vce_select = "const"
    else: vce_select = "HC1"


    ##########################################
    # analyze bins- and degree-related options
    if (isinstance(dots, bool) and not dots):
            dots = None
            dotsgrid = 0
    if (isinstance(line, bool) and not line): line = None
    if (isinstance(ci, bool) and not ci): ci = None
    if (isinstance(cb, bool) and not cb): cb = None

    # 4 cases: select J, select p, user specify both, or an error
    try: 
        len_nbins = len(nbins)
    except:
        if isinstance(nbins,(bool,int)): len_nbins = 1
        else: len_nbins = 0
    plist = pselect
    slist = sselect
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
        slist = plist.copy()
        len_s = 1
    if (len_p==0 and len_s==1):
        plist = slist.copy()
        len_p = 1

    if not isinstance(binspos,str):
        if (nbins is not None or pselect is not None or sselect is not None):
            raise Exception("binspos not correctly specified")
        
    try: 
        len_dots = len(dots)
    except:
        if isinstance(dots,(bool,int)): len_dots = 1
        else: len_dots = 0
        
    # 1st case: select J
    selection = ""
    if isinstance(binspos,str):
        if isinstance(nbins,bool):
            if nbins: selection = "J"
        else:
            if (len_nbins==1 and nbins==0): selection = "J"
            if nbins is None: selection = "J"
        if (len_nbins>1): selection = "J"
    
    if selection=="J":
        if (len_p>1 or len_s>1):
            if nbins is None:
                raise Exception("nbins must be specified for degree/smoothness selection.")
            else:
                raise Exception("Only one p and one s are allowed to select # of bins.")
        if plist is None: plist = deriv
        if slist is None: slist = plist
        if (not isinstance(dots, bool) and dots is not None):
            plist, slist = dots
            if np.isnan(slist): slist = plist.copy()
        
        if dots is None: dots = conc(plist, slist)
        if (isinstance(dots,bool) and dots): dots = conc(plist, slist)
        if (isinstance(line,bool) and line): line = conc(plist, slist)
        if (isinstance(ci,bool) and ci): ci = conc(plist+1, slist+1)
        if (isinstance(cb,bool) and cb): cb = conc(plist+1, slist+1)
        len_p = 1
        len_s = 1
    

    # 2nd case: select p (at least for one object)
    pselectOK = False
    if selection!="J":
        if dots is None: pselectOK  = True
        if (isinstance(dots,bool) and dots): pselectOK  = True
        if (isinstance(line,bool) and line): pselectOK  = True
        if (isinstance(ci,bool) and ci): pselectOK  = True
        if (isinstance(cb,bool) and cb): pselectOK  = True
    
    if (pselectOK and len_nbins==1 and (len_p>1 or len_s>1)):
        selection = "P"

    # 3rd case: user specified
    if ((len_p<=1 and len_s<=1) and selection!="J"):
        selection = "U"
        if dots is None:
            if (len_p==1 and len_s==1): dots = conc(plist, slist)
            else: dots = conc(deriv, deriv)
        
        if (isinstance(dots,bool) and dots):
            if (len_p==1 and len_s==1): dots = conc(plist, slist)
            else: dots = conc(deriv, deriv)
        
        if (np.isnan(dots[1])): dots[1] = dots[0].copy()

        if (isinstance(line,bool) and line):
            if (len_p==1 and len_s==1): line = conc(plist, slist)
            else: line = dots.copy()
        
        if (isinstance(ci,bool) and ci):
            if (len_p==1 and len_s==1): ci = conc(plist+1, slist+1)
            else: ci = conc(dots[0]+1, dots[1]+1)
        
        if (isinstance(cb,bool) and cb):
            if (len_p==1 and len_s==1): cb = conc(plist+1, slist+1)
            else: cb = conc(dots[0]+1, dots[1]+1)
        
    if (selection==""):
        raise Exception("Degree, smoothness, or # of bins not correctly specified")

    ##################################################
    ######## Error Checking ##########################
    ##################################################

    if deriv < 0:
        raise Exception("derivative incorrectly specified.")

    if dotsgrid<0 or linegrid<0 or cigrid<0 or cbgrid<0 or polyreggrid<0 or polyregcigrid<0:
        raise Exception("# of evaluation points incorrectly specified.")
        
    if not isinstance(binspos, str):
        if np.min(binspos)<=xmin or np.max(binspos)>=xmax:
            raise Exception("Knots out of allowed range.")
    else:
        if binspos!="es" and binspos!="qs":
            raise Exception("binspos incorrectly specified.")
        
    if dots is not None and len(dots)==2:
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

    # if dots is not None and dots[0] < deriv:
    #     raise Exception("p<deriv not allowed.")

    # if line is not None:
    #     if line[0] < deriv:
    #         raise Exception("p<deriv not allowed.")

    # if ci is not None:
    #     if ci[0] < deriv:
    #         raise Exception("p<deriv not allowed.")

    # if cb is not None:
    #     if cb[0] < deriv:
    #         raise Exception("p<deriv not allowed.")

    if binsmethod!="dpi" and binsmethod!="rot":
        raise Exception("Bin selection method incorrectly specified.")

    if w is not None:
        if isinstance(at,str):
            if at!="mean" and at!="median" and at!="zero":
                raise Exception("at incorrectly specified.")
        else:
            if not np.isscalar(at):
                if len(at)!=nwvar:
                    raise Exception("Length of at not equal to # of w variables.")
            
    ##################################################
    ####### Prepare options ##########################
    ##################################################
    
    if dots is None:
        dots_p = None
        dots_s = None
    else:
        dots_p, dots_s = dots
        if isinstance(dots_p, bool): dots_p = None
        if np.isnan(dots_s): dots_s = dots_p
    dotsmean = 0
    if dotsgridmean: dotsmean = 1

    if line is None: 
        linegrid = 0
        line_p = None
        line_s = None
    else:
        line_p, line_s = line
        if isinstance(line_p, bool): line_p = None
        if np.isnan(line_s): line_s = line_p

    cimean = 0
    if cigridmean: cimean = 1
    if ci is None:
        cigrid  = 0
        cimean = 0
        ci_p = None
        ci_s = None
    else:
        ci_p, ci_s = ci
        if isinstance(ci_p, bool): ci_p = None
        if np.isnan(ci_s): ci_s = ci_p

    if cb is None: 
        cbgrid = 0
        cb_p = None
        cb_s = None
    else:
        cb_p, cb_s = cb
        if isinstance(cb_p, bool): cb_p = None
        if np.isnan(cb_s): cb_s = cb_p 

    if polyreg is None:
        polyreggrid  = 0
        polyregcigrid = 0

     # Add a warning about degrees for estimation and inference
    if selection=="J":
        if ((ci_p is not None) and (ci_p<=dots_p)):
            ci_p = dots_p+1 
            ci_s = ci_p
            warnings.warn("Degree for ci has been changed. It must be greater than the degree for dots.")

        if ((cb_p is not None) and (cb_p<=dots_p)):
            cb_p = dots.p+1
            cb_s = cb.p
            warnings.warn("Degree for cb has been changed. It must be greater than the degree for dots.")
        
    if selection=="U":
        if ((ci is not None) or (cb is not None)):
            warnings.warn("Confidence intervals/bands are valid when nbins is much larger than the IMSE-optimal choice. Compare your choice with the IMSE-optimal one obtained by binsregselect().")

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

    # Extract byvals in by
    if by is not None:
        byvals = np.unique(by)
        ngroup = len(byvals)
    else:
        byvals = ["Full Sample"]
        ngroup = 1

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
    nbins_all = nbins         # "nbins" is reserved for use within loop

    if selection=="U":  selectmethod = "User-specified"
    else:
        if binsmethod=="dpi": selectmethod = "IMSE direct plug-in" 
        else: selectmethod = "IMSE rule-of-thumb"
        if selection=="J": selectmethod = selectmethod + " (select # of bins)"
        if selection=="P": selectmethod = selectmethod + " (select degree and smoothness)"

    knot = None
    knotlistON = False
    if not isinstance(binspos, str):
        nbins = len(binspos)+1
        knot = np.concatenate([[xmin], np.sort(binspos),[xmax]])
        position = "User-specified"
        es = False
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
    imse_v_rot = np.repeat(np.nan, ngroup)
    imse_v_dpi = np.repeat(np.nan, ngroup)
    imse_b_rot = np.repeat(np.nan, ngroup)
    imse_b_dpi = np.repeat(np.nan, ngroup)
    fullfewobs = selectfullON = False
    if fewmasspoints: fullfewobs = True
    if (not fullfewobs and selection!="U" and by is None) or (by is not None and samebinsby): 
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
            if dots_p is None: dotspcheck = 6
            else: dotspcheck = dots_p
            if eN <= dfcheck[0]+dotspcheck+1+qrot:
                warnings.warn("Too small effective sample size for bin selection. # of mass of points or clusters used and by option ignored.")
                fullfewobs = True
                byvals = "Full Sample"
                es = False
                position = "Quantile-spaced"
        if not fullfewobs:
            if np.isnan(Ndist): Ndist_sel = None
            else: Ndist_sel = Ndist
            if np.isnan(Nclust): Nclust_sel = None
            else: Nclust_sel = Nclust

            randcut1k = randcut
            if (randcut is None and  N>5000):
                randcut1k = max(5000/N, 0.01)
                warnings.warn("To speed up computation, bin/degree selection uses a subsample of roughly max(5,000, 0.01n) observations if the sample size n>5,000. To use the full sample, set randcut=1.")

            if (selection=="J"):
                binselect = binsregselect(y, x, w, deriv=deriv,
                                        bins=dots, binspos=binspos, nbins=nbins_all,
                                        binsmethod=binsmethod, nbinsrot=nbinsrot,
                                        vce=vce_select, cluster=cluster, randcut=randcut1k,
                                        dfcheck=dfcheck, masspoints=masspoints, weights=weights,
                                        numdist=Ndist_sel, numclust=Nclust_sel)
                if np.isnan(binselect.nbinsrot_regul):
                    raise Exception("Bin selection fails.")
                if (binsmethod == "rot"):
                    nbins = binselect.nbinsrot_regul
                    imse_v_rot = np.repeat(binselect.imse_v_rot, ngroup)
                    imse_b_rot  = np.repeat(binselect.imse_b_rot, ngroup)
                elif (binsmethod == "dpi"):
                    nbins = binselect.nbinsdpi
                    imse_v_dpi = np.repeat(binselect.imse_v_dpi, ngroup)
                    imse_b_dpi = np.repeat(binselect.imse_b_dpi, ngroup)
                    if np.isnan(nbins):
                        warnings.warn("DPI selection fails. ROT choice used.")
                        nbins = binselect.nbinsrot_regul
                        imse_v_rot = np.repeat(binselect.imse_v_rot, ngroup)
                        imse_b_rot = np.repeat(binselect.imse_b_rot, ngroup)
            elif (selection=="P"):
                binselect = binsregselect(y, x, w, deriv=deriv,
                                            binspos=binspos, nbins=nbins_all,
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
                    imse_v_rot = np.repeat(binselect.imse_v_rot, ngroup)
                    imse_b_rot = np.repeat(binselect.imse_b_rot, ngroup)
                elif binsmethod == "dpi": 
                    binsp = binselect.pdpi
                    binss =  binselect.sdpi
                    imse_v_dpi= np.repeat(binselect.imse_v_dpi, ngroup)
                    imse_b_dpi = np.repeat(binselect.imse_b_dpi, ngroup)
                    if np.isnan(binsp):
                        warnings.warn("DPI selection fails. ROT choice used.")
                        binsp  =  binselect.prot_regul
                        binss = binselect.srot_regul
                        imse_v_rot = np.repeat(binselect.imse_v.rot, ngroup)
                        imse_b_rot = np.repeat(binselect.imse_b_rot, ngroup)       
                if (isinstance(dots,bool) or dots is None):
                    dots_p = binsp
                    dots_s = binss
                if isinstance(line,bool) :
                    line_p = binsp
                    line_s = binss
                if ((ci is not None)  and (not isinstance(ci,bool)) and (ci_p<=binsp)):
                    ci_p = binsp+1
                    ci_s = ci.p
                    warnings.warn("Degree for ci has been changed. It must be greater than the IMSE-optimal degree.")
                if isinstance(ci,bool):
                    ci_p = binsp+1
                    ci_s = binss+1
                if ((cb is not None) and (not isinstance(cb,bool)) and (cb_p<=binsp)):
                    cb_p = binsp+1
                    cb_s = cb_p
                    warnings.warn("Degree for cb has been changed. It must be greater than the IMSE-optimal degree.")
                if isinstance(cb,bool):
                    cb_p = binsp+1
                    cb_s = binss+1

    # Generate knot using the full sample if needed

    if (selectfullON or (selection=="U" and samebinsby)) and not fullfewobs and knot is None:
        knotlistON = True
        nbins_full = nbins
        if es: knot = genKnot_es(xmin, xmax, nbins)
        else: knot = genKnot_qs(x, nbins)

    knot_all = None
    if knotlistON: knot_all = knot.copy()    # universal knot sequence

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
    dots_by = nanmat(ngroup,2)
    line_by = nanmat(ngroup,2)
    ci_by = nanmat(ngroup,2)
    cb_by = nanmat(ngroup,2)
    data_plot = []  # list storing graph data

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
        fewobs = False
        if (selection!="U") and (not knotlistON) and (not fullfewobs):
            # check if rot can be implemented
            if dots_p is None: dotspcheck = 6
            else: dotspcheck = dots_p
            if (nbinsrot is None) and (eN <= dfcheck[0]+dotspcheck+1+qrot):
                warnings.warn("Too small effective sample size for bin selection. # of mass points or clusters used.")
                fewobs = True
                nbins = eN
                es = False
            if not fewobs:
                if np.isnan(Ndist): Ndist_sel = None
                else: Ndist_sel = Ndist
                if np.isnan(Nclust): Nclust_sel = None
                else: Nclust_sel = Nclust
                randcut1k = randcut
                if randcut is None and N>5000:
                    randcut1k = max(5000/N, 0.01)
                    warnings.warn("To speed up computation, bin/degree selection uses a subsample of roughly max(5,000, 0.01n) observations if the sample size n>5,000. To use the full sample, set randcut=1.")
                if selection=="J":
                    binselect = binsregselect(y_sub, x_sub, w_sub, deriv=deriv,
                                            bins=dots, binspos=binspos, nbins=nbins_all,
                                            binsmethod=binsmethod, nbinsrot=nbinsrot,
                                            vce=vce_select, cluster=cluster_sub, randcut=randcut1k,
                                            dfcheck=dfcheck, masspoints=masspoints, weights=weights_sub,
                                            numdist=Ndist_sel, numclust=Nclust_sel)
                    if np.isnan(binselect.nbinsrot_regul):
                        raise Exception("Bin selection fails.")
                    if binsmethod == "rot": 
                        nbins = binselect.nbinsrot_regul
                        imse_v_rot[i] = binselect.imse_v_rot
                        imse_b_rot[i] = binselect.imse_b_rot
                    elif binsmethod == "dpi":
                        nbins = binselect.nbinsdpi
                        imse_v_dpi[i] = binselect.imse_v_dpi
                        imse_b_dpi[i] = binselect.imse_b_dpi
                        if np.isnan(nbins):
                            warnings.warn("DPI selection fails. ROT choice used.")
                            nbins = binselect.nbinsrot_regul
                            imse_v_rot[i] = binselect.imse_v_rot
                            imse_b_rot[i] = binselect.imse_b_rot

                elif selection=="P":
                    binselect = binsregselect(y_sub, x_sub, w_sub, deriv=deriv,
                                     binspos=binspos, nbins=nbins_all,
                                     pselect=plist, sselect=slist,
                                     binsmethod=binsmethod, nbinsrot=nbinsrot,
                                     vce=vce_select, cluster=cluster_sub, randcut=randcut1k,
                                     dfcheck=dfcheck, masspoints=masspoints, weights=weights_sub,
                                     numdist=Ndist_sel, numclust=Nclust_sel)
                    if np.isnan(binselect.prot_regul):
                        raise Exception("Bin selection fails.")
                    binsp = binss = np.nan
                    if (binsmethod == "rot"):
                        binsp = binselect.prot_regul
                        binss = binselect.srot_regul
                        imse_v_rot[i] = binselect.imse_v_rot
                        imse_b_rot[i] = binselect.imse_b_rot
                    elif (binsmethod == "dpi"):
                        binsp = binselect.pdpi
                        binss = binselect.sdpi
                        imse_v_dpi[i] = binselect.imse_v_dpi
                        imse_b_dpi[i] = binselect.imse_b_dpi
                        if np.isnan(binsp):
                            warnings.warn("DPI selection fails. ROT choice used.")
                            binsp = binselect.prot_regul
                            binss = binselect.srot_regul
                            imse_v_rot[i] = binselect.imse_v_rot
                            imse_b_rot[i] = binselect.imse_b_rot
                    if (isinstance(dots, bool) or (dots is None)):
                        dots_p = binsp
                        dots_s = binss
                    if isinstance(line, bool):
                        line_p = binsp
                        line_s = binss
                    if ((ci is not None) and (not isinstance(ci, bool)) and (ci_p<=binsp)):
                        ci_p = binsp+1
                        ci_s = ci_p.copy()
                        warnings.warn("Degree for ci has been changed. It must be greater than the IMSE-optimal degree.")
                    if isinstance(ci, bool):
                        ci_p = binsp+1
                        ci_s = binss+1
                    if ((cb is not None) and (not isinstance(cb, bool)) and (cb_p<=binsp)):
                        cb_p = binsp+1
                        cb_s = cb_p
                        warnings.warn("Degree for ci has been changed. It must be greater than the IMSE-optimal degree.")
                    if isinstance(cb, bool):
                        cb_p = binsp+1
                        cb_s = binss+1
                    nbins = nbins_all

        if (selection=="U"): nbins = nbins_all
        if knotlistON:
            nbins = len(knot_all)-1
            knot  = knot_all.copy()
        if fullfewobs:
            fewobs = True
            nbins = eN
        
        if (dotsmean+dotsgrid !=0): dots_by[i,:] = conc(dots_p, dots_s)
        if line is not None: line_by[i,:] = conc(line_p, line_s)
        if ci is not None: ci_by[i,:] = conc(ci_p, ci_s)
        if cb is not None: cb_by[i,:] = conc(cb_p, cb_s)


        ###########################################
        # Checking Error for each case
        dots_fewobs = line_fewobs = ci_fewobs = cb_fewobs = polyreg_fewobs = False
        if not fewobs:
            if (nbins-1)*(dots_p-dots_s+1)+dots_p+1+dfcheck[1]>=eN: 
                fewobs = True
                nbins = eN
                es = False
                warnings.warn("Too small effective sample size for dots. # of mass points or clusters used.")
            if line is not None:
                if ((nbins-1)*(line_p-line_s+1)+line_p+1+dfcheck[1]>=eN):
                    line_fewobs = True
                    warnings.warn("Too small effective sample size for line.")
            if ci is not None:
                if (nbins-1)*(ci_p-ci_s+1)+ci_p+1+dfcheck[1]>=eN:
                    ci_fewobs = True
                    warnings.warn("Too small effective sample size for ci.")
            if cb is not None:
                if (nbins-1)*(cb_p-cb_s+1)+cb_p+1+dfcheck[1]>=eN:
                    cb_fewobs = True
                    warnings.warn("Too small effective sample size for line.")

        if polyreg is not None:
            if polyreg+1>eN:
                polyreg_fewobs = True
                warnings.warn("Too small effective sample size for polynomial fit.")


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
                warnings.warn("Repeated knots. Some bins dropped.")
                nbins = len(knot)-1

        # check local mass points
        if not fewobs and localcheck:
            uniqmin = binsreg_checklocalmass(x_sub, nbins, es, knot=knot) # mimic STATA
            if dots is not None: 
                if uniqmin < dots_p+1:
                    dots_fewobs = True
                    warnings.warn("Some bins have too few distinct values of x for dots.")
            if line is not None:
                if uniqmin < line_p+1:
                    line_fewobs = True
                    warnings.warn("Some bins have too few distinct values of x for line.")
            if ci is not None:
                if uniqmin < ci_p+1:
                    ci_fewobs = True
                    warnings.warn("Some bins have too few distinct values of x for CI.")
            if cb is not None:
                if uniqmin < cb_p+1:
                    cb_fewobs = True
                    warnings.warn("Some bins have too few distinct values of x for CB.")

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
        # Dots and CI for Small eN case
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
                                weights = weights_sub, cov_type = vce,
                                cluster = cluster_sub, **optimize)
            beta = model.params[:k].values
            beta[np.isnan(beta)] = 0
            vcv = model.cov_params()

            dots_fit = beta.copy()
            if eval_w is not None:
                coeff_w = model.params[k:]
                coeff_w[np.isnan(coeff_w)] = 0
                dots_fit = dots_fit + np.sum(eval_w * coeff_w)
            data_dots = pd.DataFrame({'group': str(byvals[i]),
                                      'x': dots_x,
                                      'fit': dots_fit})
            data_by.dots = data_dots
            if cigrid+cimean!=0:
                warnings.warn("ci=(0,0) used.")
                basis_all = np.identity(len(dots_x))
                if eval_w is not None:
                    basis_all = np.column_stack((basis_all, np.outer(np.repeat(1, len(dots_x)), eval_w)))
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
            model_dots = binsreg_fit(y=y_sub, x=design, is_qreg=True, quantile = quantile, weights = weights_sub, **optimize)
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
        if linegrid !=0  and not line_fewobs and not fewobs:
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
                model_line = binsreg_fit(y=y_sub, x=design, is_qreg=True, quantile=quantile, weights=weights_sub, **optimize)
                check_drop(model_line.params, ncol(B))

            basis = binsreg_spdes(x=line_x, p=line_p, s=line_s, knot=knot, deriv=deriv)
            line_fit, line_se = binsreg_pred(basis, model_line, type = "xb",deriv=deriv, wvec=eval_w)
            if line_s == 0 or line_s-deriv <= 0: line_fit[line_isknot==1] = np.nan
            data_line = pd.DataFrame({'group': str(byvals[i]),
                                      'x': line_x,
                                      'bin': line_bin,
                                      'isknot': line_isknot,
                                      'mid': line_mid, 
                                      'fit': line_fit})
            data_by.line = data_line

        ############### Poly fit #########################
        if ((polyreggrid!=0) and (not noplot) and (not polyreg_fewobs)):
            polyON = True
        if polyON:
            if w_sub is not None:
                print("Note: When additional covariates w are included, the polynomial fit may not always be close to the binscatter fit.")
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
            model_poly = binsreg_fit(y=y_sub, x=design, is_qreg=True, quantile=quantile, weights=weights_sub, **optimize)
            beta_poly = model_poly.params
            beta_poly[np.isnan(beta_poly)] = 0
            poly_fit  = 0
            for j in range(deriv,polyreg+1):
                poly_fit = poly_fit + poly_x**(j-deriv)*beta_poly[j]*factorial(j)/factorial(j-deriv)
            if eval_w is not None and deriv==0:
                poly_fit = poly_fit + np.sum(beta_poly[polyreg+1:]*eval_w)
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
                        basis_polyci[:,j] = np.repeat(0, npolyci_x)
                
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
                model_ci = binsreg_fit(y=y_sub, x=design, is_qreg=True, quantile=quantile, weights=weights_sub, **optimize)
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
            if (nsims<2000 or simsgrid<50):
                print("Note: Setting at least nsims=2000 and simsgrid=50 is recommended to obtain the final results.")
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
                model_cb = binsreg_fit(y=y_sub, x=design, is_qreg=True, quantile=quantile, weights=weights_sub, **optimize)
                check_drop(model_cb.params, ncol(B))
            basis = binsreg_spdes(x=cb_x, p=cb_p, s=cb_s, knot=knot, deriv=deriv)
            pos = np.invert(np.isnan(model_cb.params[:ncol(basis)]))
            k_new = np.sum(pos)
            cb_fit, cb_se = binsreg_pred(X=basis, model=model_cb, type="all", deriv=deriv, wvec=eval_w, avar=asyvar)

            ### Compute cval ####
            x_grid = binsreg_grid(knot, simsgrid).eval
            basis_sim = binsreg_spdes(x=x_grid, p=cb_p, s=cb_s, knot=knot, deriv=deriv)
            sim_fit,sim_se = binsreg_pred(X=basis_sim, model=model_cb, type="all", avar=True)
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

        # save bin information
        if nbins==len(knot):
            data_by.data_bin  = pd.DataFrame({'group' : str(byvals[i]),
                                              'bin_id' : np.arange(1,nbins+1),
                                              'left_endpoint' : knot,
                                              'right.endpoint': knot})
        else:
            data_by.data_bin =  pd.DataFrame({'group' : str(byvals[i]),
                                              'bin_id' : np.arange(1,nbins+1),
                                              'left_endpoint' : knot[:-1],
                                              'right.endpoint': knot[1:]})

        # Save all data for each group
        data_plot += [data_by]
        # names(data_plot)[i] = paste("Group", byvals[i], sep=" ")  # Uncommented in R
    
    #############
    # END Loop
    #############

    ########################################
    ############# Plotting ? ###############
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
                binsplot = binsplot + geom_line(aes(x='x', y='fit'), data_by.poly.loc[index], color=bycolors[i], linetype=bylpatterns[i])         
            
            if data_by.polyci is not None:
                index = ((data_by.polyci["x"]>=xr_min) & (data_by.polyci["x"]<=xr_max) &
                         (data_by.polyci["polyci_l"]>=yr_min) & (data_by.polyci["polyci_r"]<=yr_max))
                binsplot = binsplot + geom_errorbar(aes(x='x', ymin='polyci_l', ymax='polyci_r'), data_by.polyci.loc[index], alpha=1, width = 0.05, color=bycolors[i], linetype=bylpatterns[i])
            
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
        binsplot = binsplot + labs(x=xname, y=yname) + xlim(xsc_min, xsc_max)
        print(binsplot)

    ######################################
    ########### Output ###################
    ######################################

    opt = options_qreg(dots=dots_by, line=line_by, ci=ci_by, cb=cb_by,
                        polyreg=polyreg, deriv=deriv, quantile=quantile,
                        binspos=position, binsmethod=selectmethod,
                        N_by=N_by, Ndist_by=Ndist_by, Nclust_by=Nclust_by,
                        nbins_by=nbins_by, byvals=byvals)

    out = binsqreg_output(bins_plot=binsplot, data_plot=data_plot, cval_by=cval_by,
                            imse_v_dpi=imse_v_dpi, imse_b_dpi=imse_b_dpi,
                            imse_v_rot=imse_v_rot, imse_b_rot=imse_b_rot,
                            options = opt)
    return out
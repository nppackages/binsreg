#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Thu Mar 16 10:20:24 2023
# @author: Ricardo Masini

import numpy as np
import pandas as pd
import warnings
from binsreg.funs import *

def binsregselect(y, x, w=None, data=None, deriv=0, bins=None, pselect=None, sselect=None,
                    binspos="qs", nbins=None,binsmethod="dpi", nbinsrot=None, simsgrid=20, savegrid=False,
                    vce="HC1", useeffn=None, randcut=None, cluster=None,
                    dfcheck=(20,30), masspoints="on", weights=None, subset=None,
                    norotnorm=False, numdist=None, numclust=None):

    '''
    Data-Driven IMSE-Optimal Partitioning/Binning Selection for Binscatter.
    
    Description
    -----------
    binsregselect implements data-driven procedures for selecting the number of bins for binscatter
    estimation. The selected number is optimal in minimizing integrated mean squared error (IMSE).

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
    
    deriv : int 
        Derivative order of the regression function for estimation, testing and plotting.
        The default is deriv=0, which corresponds to the function itself.

    bins : tuple
        bins=(p,s) set a piecewise polynomial of degree p with s smoothness constraints
        for data-driven (IMSE-optimal) selection of the partitioning/binning scheme. The default is
        bins=(0,0), which corresponds to piecewise constant (canonical binscatter).
    
    pselect: array 
        Vector of integers within which the degree of polynomial p for point estimation is selected.
        Note: To implement the degree or smoothness selection, in addition to pselect or sselect
        nbins=# must be specified.

    sselect: array 
        Vector of integers within which the number of smoothness constraints s for point estimation is selected.
        If not specified, for each value p supplied in the option pselect, only the
        piecewise polynomial with the maximum smoothness is considered, i.e., s=p.

    binspos : array
        Position of binning knots. The default is binspos="qs", which corresponds to quantile-spaced
        binning (canonical binscatter). The other options is binspos="es" for evenly-spaced binning.
    
    nbins: int
        Number of bins for degree/smoothness selection. If nbins=True or nbins=None (default) is specified,
        the function selects the number of bins instead, given the specified degree and smoothness.
        If a vector with more than one number is specified, the command selects the number of bins within this vector.
    
    binsmethod : str
        Method for data-driven selection of the number of bins. The default is binsmethod="dpi",
        which corresponds to the IMSE-optimal direct plug-in rule.  The other option is: "rot"
        for rule of thumb implementation.

    nbinsrot : int
        Initial number of bins value used to construct the DPI number of bins selector.
        If not specified, the data-driven ROT selector is used instead.

    simsgrid : int
        Number of evaluation points of an evenly-spaced grid within each bin used for evaluation of
        the supremum (infimum or Lp metric) operation needed to construct confidence bands and hypothesis testing
        procedures. The default is simsgrid=20, which corresponds to 20 evenly-spaced
        evaluation points within each bin for approximating the supremum (infimum or Lp metric) operator.
    
    savegrid : bool
        If true, a data frame produced containing grid.
    
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
     
    useeffn : int
        Effective sample size to be used when computing the (IMSE-optimal) number of bins. This option
        is useful for extrapolating the optimal number of bins to larger (or smaller) datasets than
        the one used to compute it.
    
    randcut : float
        Upper bound on a uniformly distributed variable used to draw a subsample for bins/degree/smoothness selection.
        Observations for which numpy.random.uniform()<=# are used. # must be between 0 and 1.
    
    cluster: array 
        Cluster ID. Used for compute cluster-robust standard errors.
    
    dfcheck :tuple
        Adjustments for minimum effective sample size checks, which take into account number of unique
        values of x (i.e., number of mass points), number of clusters, and degrees of freedom of
        the different statistical models considered. The default is dfcheck=(20, 30).
        See Cattaneo, Crump, Farrell and Feng (2024c) for more details.
    
    masspoints: str
        How mass points in x are handled. Available options:
            * "on"           : all mass point and degrees of freedom checks are implemented. Default.
            * "noadjust"     :  mass point checks and the corresponding effective sample size adjustments are omitted.
            * "nolocalcheck" :  within-bin mass point and degrees of freedom checks are omitted.
            * "off"          : "noadjust" and "nolocalcheck" are set simultaneously.
            * "veryfew"      : forces the function to proceed as if x has only a few number of mass points (i.e., distinct values).
                               In other words, forces the function to proceed as if the mass point and degrees of freedom checks were failed.

    norotnorm : bool
        If true, a uniform density rather than normal density used for ROT selection.
    
    numdist : int
        Number of distinct for selection. Used to speed up computation.
    
    numclust : int 
        Number of clusters for selection. Used to speed up computation.
    
    weights : array
        An optional vector of weights to be used in the fitting process. Should be None or
        a numeric vector.
    
    subset : array
        Optional rule specifying a subset of observations to be used.
    
    Returns
    -------
    nbinsrot_poly : ROT number of bins, unregularized.
    
    nbinsrot_regul : ROT number of bins, regularized.
    
    nbinsrot_uknot : ROT number of bins, unique knots.
    
    nbinsdpi : DPI number of bins.
    
    nbinsdpi_uknot : DPI number of bins, unique knots.
    
    prot_poly : ROT degree of polynomials, unregularized.

    prot_regul : ROT degree of polynomials, regularized.
    
    prot_uknot : ROT degree of polynomials, unique knots.
    
    pdpi : DPI degree of polynomials.
    
    pdpi_uknot : DPI degree of polynomials, unique knots.

    srot_poly : ROT number of smoothness constraints, unregularized.

    srot_regul : ROT number of smoothness constraints, regularized.

    srot_uknot : ROT number of smoothness constraints, unique knots.

    sdpi : DPI number of smoothness constraints.

    sdpi_uknot : DPI number of smoothness constraints, unique knots.

    imse_v_rot : Variance constant in IMSE expansion, ROT selection.

    imse_b_rot : Bias constant in IMSE expansion, ROT selection.

    imse_v_dpi : Variance constant in IMSE expansion, DPI selection.

    imse_b_dpi : Bias constant in IMSE expansion, DPI selection.

    int.result : Intermediate results, including a matrix of degree and smoothness (deg_mat),
                 the selected numbers of bins (vec_nbinsrot_poly, vec_nbinsrot_regul,vec_nbinsrot_uknot,
                 vec_nbinsdpi, vec_nbinsdpi.uknot), and the bias and variance constants in 
                 IMSE (vec_imse_b_rot, vec_imse_v_rot, vec_imse_b_dpi, vec_imse_v_dpi) 
                 under each rule (ROT or DPI), corresponding to each pair of degree and smoothness
                 each row in deg_mat.
    
    opt :  A list containing options passed to the function, as well as total sample size (n),
           number of distinct values (Ndist) in x, and number of clusters (Nclust).
    
    data_grid : A data frame containing grid.


    See Also
    --------
    binsreg, binsglm, binsqreg, binstest, binspwc.
    
    Example
    ------- 
    >>> x = numpy.random.uniform(size = 500)
    >>> y = numpy.sin(x) + numpy.random.normal(size = 500)
    >>> out = binsregselect(y,x)
    >>> print(out)
    >>> out.summary()
    '''

    # param for internal use
    rot_lb = 1
    qrot = 2

    # extract x, y, w, cluster, weights and subset from data (if data is supplied)
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

    # Reshaping x, y, w, cluster and weights
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

    ## analyze bins- and degree-related options
    if bins is not None:
        bins_p = bins[1]
        if (len(bins)==1): bins_s = bins_p
        elif (len(bins)==2): bins_s = bins[1]
        else: raise Exception("Bins not correctly specified.")
        plist = bins[0]; slist = bins[1]
    else:
        plist = pselect; slist = sselect

    try: len_p = len(plist)
    except: len_p=0
    try: len_s = len(slist)
    except: len_s=0

    if (len_p==1 & len_s==0):
        slist = plist
        len_s = 1
    if (len_p==0 & len_s==1):
        plist = slist
        len_p = 1

    selectJ = None
    if nbins is True: selectJ = True
    if nbins is None: selectJ = True
    else:
        try: len_nbins = len(nbins)
        except:
            if np.isscalar(nbins): len_nbins = 1
            else: len_nbins = 0
        if (len_nbins==1 and nbins==0): selectJ = True
        if (len_nbins>1 or bins is not None or (len_p==1 and len_s==1)): selectJ = True
    
    if isinstance(selectJ,bool) and selectJ:
        if (len_p>1 or len_s>1): raise Exception("Only one p and one s are allowed for J selection.")
        if (plist is None):
            plist = deriv
            len_p = 1
        if (slist is None):
            slist = plist
            len_s = 1
    
    if (selectJ is None and (len_p>1 or len_s>1)): selectJ = False

    if selectJ is None: raise Exception("Degree, smoothness, or # of bins are not correctly specified.")
  
    # find all compatible pairs
    if selectJ: deg_mat = cbind(plist, slist)
    else:
        if (len_p>0 and len_s==0): deg_mat = cbind(plist, plist)
        if (len_p==0 and len_s>0): deg_mat = cbind(slist, slist)
        if (len_p>0 and len_s>0):
            deg_mat = np.array([(p,s) for p in plist for s in slist])
            comp_ind = (deg_mat[:,0]>=deg_mat[:,1]) & (deg_mat[:,1]>deriv)   # p>=s and p>deriv
            if (sum(comp_ind)==0): 
                raise Exception("Degree and smoothness incompatible")
            deg_mat = deg_mat[comp_ind,:]
        
    # colnames(deg_mat) = = ("degree", "smoothness")

    ########################## error checking
    
    if bins is not None and len(bins)==2:
        if bins[0] < bins[1]:
            raise Exception("p<s not allowed.")
    
    if binsmethod!="dpi" and binsmethod!="rot":
        raise Exception("bin selection method incorrectly specified.")
    
    if binspos!="es" and binspos!="qs":
        raise Exception("binspos incorrectly specified.")
    
    #####################################
    rot_fewobs = dpi_fewobs = False
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
    elif masspoints=="veryfew":
        rot_fewobs = dpi_fewobs = True

    # effective size
    eN = N =len(x)
    Ndist = np.nan
    if massadj:
        if numdist is not None: Ndist =numdist
        else: Ndist = len(np.unique(x))
        eN = min(eN, Ndist)
    Nclust = np.nan
    if cluster is not None: 
        if numclust is not None: Nclust =numclust
        else: Nclust = len(np.unique(cluster))
        eN = min(eN, Nclust)
    
    # take a subsample?
    x_full = x.copy()
    eN_sub = eN
    Ndist_sub = Ndist
    Nclust_sub = Nclust
    if randcut is not None:
        subsample =(np.random.uniform(size = N)<=randcut)
        x = x[subsample]
        y = y[subsample]
        if w is not None: w = w[subsample,:]
        if weights is not None: weights = weights[subsample]
        if cluster is not None: cluster = cluster[subsample]

        eN_sub = len(x)
        if massadj:
            Ndist_sub = len(np.unique(x))
            eN_sub = min(eN_sub, Ndist_sub)
        if cluster is not None:
            Nclust_sub =len(np.unique(cluster))
            eN_sub = min(eN_sub, Nclust_sub)


    # redefine cluster if needed
    if massadj:
        if (cluster is None): cluster = x.copy()
        # note: clustered at mass point level
        else:
            if (Nclust_sub > Ndist_sub):
                cluster = x.copy()
                warnings.warn("# mass points < # clusters. Clustered at mass point level.")
    
    # Prepare params
    # p,s = bins

    if binspos == "es": es = True
    else: es = False

    # Store options
    if es: position = "Evenly-spaced"
    else: position = "Quantile-spaced"

    if binsmethod == "dpi": selectmethod = "IMSE direct plug-in"
    else: selectmethod = "IMSE rule-of-thumb"

    if selectJ: selectmethod = selectmethod + " (select # of bins)"
    else: selectmethod = selectmethod + " (select degree and smoothness)"

    nj = len(deg_mat)
    vec_J_rot_poly = nanmat(nj)
    vec_J_rot_regul = nanmat(nj)
    vec_J_rot_uniq = nanmat(nj)
    vec_J_dpi = nanmat(nj)
    vec_J_dpi_uniq = nanmat(nj)
    vec_imse_v_rot = nanmat(nj)
    vec_imse_b_rot = nanmat(nj)
    vec_imse_v_dpi = nanmat(nj)
    vec_imse_b_dpi = nanmat(nj)
    #### START loop here #########
    for j in range(nj):
        p,s = deg_mat[j,:]
        # Run rot selection
        J_rot_regul = J_rot_poly = imse_v_rot = imse_b_rot = np.nan
        if nbinsrot is not None: J_rot_regul = nbinsrot
        if np.isnan(J_rot_regul) and not rot_fewobs:
            # checking
            if eN_sub <= dfcheck[0]+p+1+qrot:
                rot_fewobs = True
                warnings.warn("Too small effective sample size for bin selection.")
            if not rot_fewobs:
                J_rot_poly, imse_b_rot, imse_v_rot = binsregselect_rot(y, x, w, p, s, deriv, eN=eN_sub, es=es,
                                                qrot=qrot, norotnorm=norotnorm, weights=weights)
            J_rot_regul = max(J_rot_poly, int(np.ceil((2*(p+1-deriv)/(1+2*deriv)*rot_lb*eN_sub)**(1/(2*p+3)))))

        # repeated knots?
        J_rot_uniq = J_rot_regul
        if not es and not np.isnan(J_rot_regul):
            J_rot_uniq = len(np.unique(genKnot_qs(x, J_rot_regul)[1:]))

        # Run dpi selection
        J_dpi = imse_v_dpi = imse_b_dpi = np.nan
        if binsmethod == "dpi" and not dpi_fewobs:
            # check if dpi can be implemented
            if not np.isnan(J_rot_uniq):
                if (p-s+1)*(J_rot_uniq-1)+p+2+dfcheck[1]>=eN_sub:
                    dpi_fewobs = True
                    warnings.warn("Too small effective sample size for DPI selection.")

                # check empty bins
                if localcheck:
                    uniqmin = binsreg_checklocalmass(x, J_rot_regul, es, knot = None) # mimic STATA
                    if uniqmin < p+2:
                        dpi_fewobs = True
                        warnings.warn("Some bins have too few distinct values of x for DPI selection.")
            else: dpi_fewobs = True

            if not dpi_fewobs:
                J_dpi, imse_v_dpi, imse_b_dpi = binsregselect_dpi(y, x, w, p, s, deriv, es=es, vce=vce, cluster=cluster, nbinsrot=J_rot_uniq, weights=weights)
                imse_v_dpi = imse_v_dpi * eN_sub
                                                                
        J_dpi_uniq = J_dpi

        if useeffn is not None or randcut is not None:
            if useeffn is not None: scaling = (useeffn/eN)**(1/(2*p+2+1))
            if randcut is not None: scaling = (eN/eN_sub)**(1/(2*p+2+1))
            if not np.isnan(J_rot_poly):  J_rot_poly  =int(np.ceil(J_rot_poly * scaling))
            if not np.isnan(J_rot_regul): J_rot_regul =int(np.ceil(J_rot_regul * scaling))
            if not np.isnan(J_rot_uniq):  J_rot_uniq  =int(np.ceil(J_rot_uniq * scaling))
            if not np.isnan(J_dpi):       J_dpi       =int(np.ceil(J_dpi * scaling))
            if not np.isnan(J_dpi_uniq):  J_dpi_uniq  =int(np.ceil(J_dpi_uniq * scaling))

        # save results
        vec_J_rot_poly[j] = J_rot_poly
        vec_J_rot_regul[j] = J_rot_regul
        vec_J_rot_uniq[j] = J_rot_uniq
        vec_J_dpi[j] = J_dpi
        vec_J_dpi_uniq[j] = J_dpi_uniq
        vec_imse_b_rot[j] = imse_b_rot
        vec_imse_v_rot[j] = imse_v_rot
        vec_imse_b_dpi[j] = imse_b_dpi
        vec_imse_v_dpi[j] = imse_v_dpi
    ##### END loop here  ########

    # extract results
    imse_b_dpi_upd = np.nan 
    imse_v_dpi_upd = np.nan 
    if selectJ:
        if  (nbins is not None) and (not isinstance(nbins, bool)) and len(nbins)>1:
            J_rot_poly  = nbins[which_min(np.abs(vec_J_rot_poly-nbins))]
            J_rot_regul = nbins[which_min(np.abs(vec_J_rot_regul-nbins))]
            J_rot_uniq  = nbins[which_min(np.abs(vec_J_rot_uniq-nbins))]
            J_dpi       = nbins[which_min(np.abs(vec_J_dpi-nbins))]
            J_dpi_uniq  = nbins[which_min(np.abs(vec_J_dpi_uniq-nbins))]
        ord_rot_poly = deg_mat.copy().reshape(-1,)
        ord_rot_regul = deg_mat.copy().reshape(-1,)
        ord_rot_uniq = deg_mat.copy().reshape(-1,)
        ord_dpi = deg_mat.copy().reshape(-1,)
        ord_dpi_uniq = deg_mat.copy().reshape(-1,)
    else:
        ind_rot_poly = which_min(np.abs(vec_J_rot_poly-nbins))
        ind_dpi      = which_min(np.abs(vec_J_dpi-nbins))

        ord_rot_poly  = deg_mat[ind_rot_poly,:]
        ord_rot_regul = deg_mat[which_min(np.abs(vec_J_rot_regul-nbins)),:]
        ord_rot_uniq  = deg_mat[which_min(np.abs(vec_J_rot_uniq-nbins)),:]
        ord_dpi       = deg_mat[ind_dpi,:]
        ord_dpi_uniq  = deg_mat[which_min(abs(vec_J_dpi_uniq-nbins)),:]
        J_rot_poly = nbins
        J_rot_regul = nbins
        J_rot_uniq = nbins
        J_dpi = nbins
        J_dpi_uniq = nbins

        imse_v_rot = vec_imse_v_rot[ind_rot_poly]
        imse_b_rot = vec_imse_b_rot[ind_rot_poly]
        imse_v_dpi = vec_imse_v_dpi[ind_dpi]
        imse_b_dpi = vec_imse_b_dpi[ind_dpi]

        if (nbins!=vec_J_dpi[ind_dpi]):
            J_dpi, imse_v_dpi, imse_b_dpi = binsregselect_dpi(y, x, w, ord_dpi[0], ord_dpi[1], deriv, es=es, vce=vce, cluster=cluster, nbinsrot=nbins, weights=weights)
            imse_b_dpi_upd = imse_b_dpi
            imse_v_dpi_upd = imse_v_dpi * eN_sub
            imse_b_dpi = imse_b_dpi_upd
            imse_v_dpi = imse_v_dpi_upd

    
    # Generate a knot vector
    if binsmethod == "rot": Jselect = J_rot_uniq
    else: Jselect = J_dpi
    
    knot = data_grid = np.nan
    if not np.isnan(Jselect) and useeffn is None:
        if es: knot = genKnot_es(np.min(x_full), np.max(x_full), Jselect)
        else: knot = genKnot_qs(x_full, Jselect)
        knot = np.append(knot[0],np.unique(knot[1:]))
        Jselect = len(knot)-1
        if binsmethod=="dpi": J_dpi_uniq =Jselect
        
        # a grid dataset
        if savegrid:
            grid = binsreg_grid(knot=knot, ngrid=simsgrid, addmore=True)
            data_grid = np.column_stack((grid.eval, grid.bin, grid.isknot))
            # data_grid = np.column_stack((data_grid, np.zeros(nrow(data_grid),len(wname))))
            data_grid = pd.DataFrame(data_grid)
            data_grid.columns = (xname, "binreg_bin", "binsreg_isknot")


    ######################
    #######output#########
    ######################

    int_result=inter_result(vec_nbinsrot_poly=vec_J_rot_poly, vec_nbinsrot_regul=vec_J_rot_regul,
                    vec_nbinsrot_uknot=vec_J_rot_uniq, vec_nbinsdpi=vec_J_dpi, vec_nbinsdpi_uknot=vec_J_dpi_uniq,
                    vec_imse_b_rot=vec_imse_b_rot, vec_imse_v_rot=vec_imse_v_rot,
                    vec_imse_b_dpi=vec_imse_b_dpi, vec_imse_v_dpi=vec_imse_v_dpi,
                    deg_mat=deg_mat)

    opt = options_select(deriv=deriv,selectJ=selectJ, binspos=position, 
                         binsmethod=selectmethod,n=N, Ndist=Ndist, Nclust=Nclust)
    
    out = binsregselect_output(nbinsrot_poly=J_rot_poly, nbinsrot_regul=J_rot_regul, nbinsrot_uknot=J_rot_uniq,
                                nbinsdpi=J_dpi, nbinsdpi_uknot=J_dpi_uniq, imse_b_rot=imse_b_rot, 
                                imse_v_rot=imse_v_rot, imse_b_dpi=imse_b_dpi, imse_v_dpi=imse_v_dpi,
                                int_result = int_result,
                                prot_poly=ord_rot_poly[0], srot_poly=ord_rot_poly[1],
                                prot_regul=ord_rot_regul[0], srot_regul=ord_rot_regul[1],
                                prot_uknot=ord_rot_uniq[0], srot_uknot=ord_rot_uniq[1],
                                pdpi=ord_dpi[0], sdpi=ord_dpi[1],
                                pdpi_uknot=ord_dpi_uniq[0], sdpi_uknot=ord_dpi_uniq[1],
                                options = opt, knot=knot, data_grid=data_grid)
    return out
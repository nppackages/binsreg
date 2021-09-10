#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Sat Sep  4 17:39:50 2021
# @author: Ricardo Masini
# binsreg package supporting classes and functions

import numpy as np
from numpy.linalg import solve
from scipy.stats import norm
import warnings
import statsmodels.api as sm

#######################################################
########## Classes Definitons #########################
#######################################################

class pwc_output:
    def __init__(self, tstat, pval, options):
        self.tstat = tstat
        self.pval = pval
        self.options = options
    
    def __repr__(self):
        print("Call: binspwc\n")
        print("Pairwise Group Comparison")
        fw = 40
        fr = 15
        print("Group Variable                     =  ", "t")
        print("Estimation Method (estmethod)      =  ", self.options.estmethod)
        if self.options.estmethod=="generalized linear model":
            print("Distribution                       =  ", self.options.dist)
            print("Link                               =  ", self.options.link)
        if self.options.estmethod=="quantile regression":
            print("Quantile                           =  ", self.options.quantile)
        print("degree (p)                         =  ", str(self.options.pwc_p).ljust(fr))
        print("smooth (s)                         =  ", str(self.options.pwc_s).ljust(fr))
        print("Derivative (deriv)                 =  ", str(self.options.deriv).ljust(fr))
        print("Bin selection:")
        print("  Method (binsmethod)              =  ", self.options.binsmethod.ljust(fr))
        print("  Placement (binspos)              =  ", self.options.binspos.ljust(fr))
        print("  degree (p)                       =  ", str(self.options.bins_p).ljust(fr))
        print("  smooth (s)                       =  ", str(self.options.bins_s).ljust(fr))
        print("")
        c = (20,12,12)
        cc = sum(c)+2
        for i in range(nrow(self.tstat)):
            g1, g2 = self.tstat[i,1:].astype(int)
            group1 = self.options.byvals[g1]
            group2 = self.options.byvals[g2]
            print(f"Group {group1} vs. Group{group2}")
            print("="*cc)
            print("Group".ljust(c[0]), f"{group1}".rjust(c[1]),f"{group2}".rjust(c[2]))
            print("-"*cc)
            print("Sample size".ljust(c[0]), f"{self.options.N_by[g1]}".rjust(c[1]),f"{self.options.N_by[g2]}".rjust(c[2]))
            print("# of distinct values".ljust(c[0]), f"{self.options.Ndist_by[g1]}".rjust(c[1]),f"{self.options.Ndist_by[g2]}".rjust(c[2]))
            print("# of clusters".ljust(c[0]), f"{self.options.Nclust_by[g1]}".rjust(c[1]),f"{self.options.Nclust_by[g2]}".rjust(c[2]))
            print("# of bins".ljust(c[0]), f"{self.options.nbins_by[g1]}".rjust(c[1]),f"{self.options.nbins_by[g2]}".rjust(c[2]))
            print("-"*cc)
            print('')
            
            print(f"diff = Group {group1} - Group {group2}")
            print("="*cc)
            if self.options.testtype=="left":
                print("H0:".ljust(c[0]), "sup T".rjust(c[1]),"p value".rjust(c[2]))
                print("-"*cc)
                print("diff<=0".ljust(c[0]), "{:3.3f}".format(self.tstat[i,0]).rjust(c[1]),"{:3.3f}".format(self.pval[i,0]).rjust(c[2]))
                print("-"*cc)
            elif self.options.testtype=="right":
                print("H0:".ljust(c[0]), "inf T".rjust(c[1]),"p value".rjust(c[2]))
                print("-"*cc)
                print("diff>=0".ljust(c[0]), "{:3.3f}".format(self.tstat[i,0]).rjust(c[1]),"{:3.3f}".format(self.pval[i,0]).rjust(c[2]))
                print("-"*cc)
            else:
                if not np.isfinite(self.options.lp):
                    print("H0:".ljust(c[0]), "sup |T|".rjust(c[1]),"p value".rjust(c[2]))
                    print("-"*cc)
                    print("diff=0".ljust(c[0]), "{:3.3f}".format(self.tstat[i,0]).rjust(c[1]),"{:3.3f}".format(self.pval[i,0]).rjust(c[2]))
                    print("-"*cc)
                else:
                    print("H0:".ljust(c[0]), f"L-{self.options.lp} of T".rjust(c[1]),"p value".rjust(c[2]))
                    print("-"*cc)
                    print("diff=0".ljust(c[0]), "{:3.3f}".format(self.tstat[i,0]).rjust(c[1]),"{:3.3f}".format(self.pval[i,0]).rjust(c[2]))
                    print("-"*cc)
        return('')

class test_output:
    def __init__(self,testshapeL,testshapeR,testshape2,
                testpoly,testmodel,options):
        self.testshapeL = testshapeL
        self.testshapeR = testshapeR
        self.testshape2 = testshape2
        self.testpoly = testpoly
        self.testmodel = testmodel
        self.options = options

    def __repr__(self):
        print("Call: binstest\n")
        fw = 40
        fr = 15
        print("Sample size (n)                    =  ", str(self.options.n).ljust(fr))
        print("# of distinct values (Ndist)       =  ", str(self.options.Ndist).ljust(fr))
        print("# of clusters (Nclust)             =  ", str(self.options.Nclust).ljust(fr))
        print("Estimation Method (estmethod)      =  ", self.options.estmethod)
        print("Derivative (deriv)                 =  ", str(self.options.deriv).ljust(fr))
        print("Bin selection:")
        print("  Method (binsmethod)              =  ", self.options.binsmethod.ljust(fr))
        print("  Placement (binspos)              =  ", self.options.binspos.ljust(fr))
        print("  degree (p)                       =  ", str(self.options.bins_p).ljust(fr))
        print("  smooth (s)                       =  ", str(self.options.bins_s).ljust(fr))
        print("  # of bins (nbins)                =  ", str(self.options.nbins).ljust(fr))
        print("")

        c = (15,12,12)
        cc = sum(c)+2
        if self.testshapeL.val is not None or self.testshapeR.val is not None or self.testshape2.val is not None:
            print("Shape Restriction Tests:")
            print("degree (p) = ", self.options.testshape[0], "; smooth (s) = ", self.options.testshape[1])

            if self.testshapeL.val is not None:
                print("="*cc)
                print("H0: sup mu <=".ljust(c[0]), "sup T".rjust(c[1]),"p value".rjust(c[2]))
                print("-"*cc)
                for i in range(len(self.testshapeL.val)):
                    print(str(self.testshapeL.val[i]).center(c[0]),
                        "{:3.3f}".format(self.testshapeL.stat[i]).rjust(c[1]),
                        "{:3.3f}".format(self.testshapeL.pval[i]).rjust(c[2]))
                print("-"*cc)
                print("\n")

            if self.testshapeR.val is not None:
                print("="*cc)
                print("H0: inf mu >=".ljust(c[0]), "inf T".rjust(c[1]),"p value".rjust(c[2]))
                print("-"*cc)
                for i in range(len(self.testshapeR.val)):
                    print(str(self.testshapeR.val[i]).center(c[0]),
                        "{:3.3f}".format(self.testshapeR.stat[i]).rjust(c[1]),
                        "{:3.3f}".format(self.testshapeR.pval[i]).rjust(c[2]))
                print("-"*cc)
                print("\n")

            if self.testshape2.val is not None:
                print("="*cc)
                if np.isfinite(self.options.lp):
                    norm = f"L-{self.options.lp} of |T|"
                else:
                    norm = f"sup |T|"
                print("H0: mu =".ljust(c[0]), norm.rjust(c[1]),"p value".rjust(c[2]))
                print("-"*cc)
                for i in range(len(self.testshape2.val)):
                    print(str(self.testshape2.val[i]).center(c[0]),
                        "{:3.3f}".format(self.testshape2.stat[i]).rjust(c[1]),
                        "{:3.3f}".format(self.testshape2.pval[i]).rjust(c[2]))
                print("-"*cc)
                print("\n")

        try: dum = all(np.isnan(self.testmodel.stat))
        except: dum = np.isnan(self.testmodel.stat)

        if self.testpoly.val is not None or not dum:
            print("Model Specification Tests:")
            print("degree (p) = ", self.options.testmodel[0], "; smooth (s) = ", self.options.testmodel[1])

            if self.testpoly.val is not None:
                print("="*cc)
                if np.isfinite(self.options.lp):
                    norm = f"L-{self.options.lp} of |T|"
                else:
                    norm = f"sup |T|"
                print("H0: mu =".ljust(c[0]), norm.rjust(c[1]),"p value".rjust(c[2]))
                print("-"*cc)
                print(("poly. p="+str(self.testpoly.val)).ljust(c[0]),
                        "{:3.3f}".format(self.testpoly.stat[0]).rjust(c[1]),
                        "{:3.3f}".format(self.testpoly.pval[0]).rjust(c[2]))
                print("-"*cc)

            if not dum:
                print("="*cc)
                if np.isfinite(self.options.lp):
                    norm = f"L-{self.options.lp} of |T|"
                else:
                    norm = f"sup |T|"
                print("H0: mu =".ljust(c[0]), norm.rjust(c[1]),"p value".rjust(c[2]))
                print("-"*cc)
                for i in range(len(self.testmodel.stat)):
                    print(f"model{i+1}".center(c[0]),
                        "{:3.3f}".format(self.testmodel.stat[0]).rjust(c[1]),
                        "{:3.3f}".format(self.testmodel.pval[0]).rjust(c[2]))
                print("-"*cc)
        return ''

class test_str:
    def __init__(self,val,stat,pval):
        self.val = val
        self.stat = stat
        self.pval = pval

class options_pwc:
    def __init__(self, bins_p, bins_s, deriv, pwc_p, pwc_s, byname,
                testtype, binspos, binsmethod, N_by, Ndist_by, Nclust_by, nbins_by, byvals,
                lp, estmethod, quantile, dist, link):
        self.bins_p = bins_p
        self.bins_s = bins_s
        self.deriv = deriv
        self.pwc_p = pwc_p
        self.pwc_s = pwc_s
        self.byname = byname
        self.testtype = testtype
        self.binspos = binspos
        self.binsmethod = binsmethod
        self.N_by = N_by
        self.Ndist_by = Ndist_by
        self.Nclust_by = Nclust_by
        self.nbins_by = nbins_by
        self.byvals = byvals
        self.lp = lp
        self.estmethod = estmethod
        self.quantile = quantile
        self.dist = dist
        self.link = link

class options_test:
    def __init__(self, bins_p, bins_s, deriv, testshape, testmodel,
                binspos, binsmethod,n, Ndist, Nclust, nbins, knot,lp,
                estmethod, quantile, dist, link):
        self.bins_p = bins_p
        self.bins_s = bins_s
        self.deriv = deriv
        self.testshape = testshape
        self.testmodel = testmodel
        self.binspos = binspos
        self.binsmethod = binsmethod
        self.n = n
        self.Ndist = Ndist
        self.Nclust = Nclust
        self.nbins = nbins
        self.knot = knot
        self.lp = lp
        self.estmethod = estmethod
        self.quantile = quantile
        self.dist = dist
        self.link = link

class options_select:
    def __init__(self, bins_p, bins_s, deriv, binspos, binsmethod,
                    n, Ndist, Nclust):
        self.bins_p = bins_p
        self.bins_s = bins_s
        self.deriv = deriv
        self.binspos = binspos
        self.binsmethod = binsmethod
        self.n = n
        self.Ndist = Ndist
        self.Nclust = Nclust

class binsregselect_output:
    def __init__(self, nbinsrot_poly, nbinsrot_regul, nbinsrot_uknot, nbinsdpi, nbinsdpi_uknot,
                    imse_b_rot, imse_v_rot, imse_b_dpi, imse_v_dpi, options, knot,data_grid):
        self.nbinsrot_poly = nbinsrot_poly
        self.nbinsrot_regul = nbinsrot_regul
        self.nbinsrot_uknot =nbinsrot_uknot
        self.nbinsdpi = nbinsdpi
        self.nbinsdpi_uknot = nbinsdpi_uknot
        self.imse_b_rot = imse_b_rot
        self.imse_v_rot = imse_v_rot
        self.imse_b_dpi = imse_b_dpi
        self.imse_v_dpi = imse_v_dpi
        self.options = options
        self.knot = knot
        self.data_grid = data_grid
    
    def __repr__(self):
        print("Call: binsregselect\n")
        fw = 40
        fr = 15
        print("Sample size (n)                    =  ", str(self.options.n).ljust(fr))
        print("# of distinct values (Ndist)       =  ", str(self.options.Ndist).ljust(fr))
        print("# of clusters (Nclust)             =  ", str(self.options.Nclust).ljust(fr))
        print("Derivative (deriv)                 =  ", str(self.options.deriv).ljust(fr))
        print("Bin selection:")
        print("  Method (binsmethod)              =  ", self.options.binsmethod.ljust(fr))
        print("  Placement (binspos)              =  ", self.options.binspos.ljust(fr))
        print("  degree (p)                       =  ", str(self.options.bins_p).ljust(fr))
        print("  smooth (s)                       =  ", str(self.options.bins_s).ljust(fr))
        print("")

        fw = 10
        fr = 10
        ft = 55
        n_dec = 3
        print('='*ft)
        print("method".ljust(fw),"# bins".rjust(fr),"df".rjust(fr),"imse-bias^2".rjust(fr), "imse-var".rjust(fr))
        print('-'*ft)

        ROT_poly_df = np.nan
        if not np.isnan(self.nbinsrot_poly):
            ROT_poly_df = self.options.bins_p+1+(self.nbinsrot_poly-1)*(self.options.bins_p-self.options.bins_s+1)
        print('ROT-POLY'.ljust(fw),str(self.nbinsrot_poly).rjust(fr),str(ROT_poly_df).rjust(fr),
                ("{:.{}f}".format(self.imse_b_rot, n_dec)).rjust(fr),("{:.{}f}".format(self.imse_v_rot, n_dec)).rjust(fr))

        ROT_regul_df = np.nan
        if not np.isnan(self.nbinsrot_regul):
            ROT_regul_df = self.options.bins_p+1+(self.nbinsrot_regul-1)*(self.options.bins_p-self.options.bins_s+1)
        print('ROT-REGUL'.ljust(fw),str(self.nbinsrot_regul).rjust(fr),str(ROT_regul_df).rjust(fr),
                str(np.nan).rjust(fr),str(np.nan).rjust(fr))

        ROT_uknot_df = np.nan
        if not np.isnan(self.nbinsrot_uknot):
            ROT_uknot_df = self.options.bins_p+1+(self.nbinsrot_uknot-1)*(self.options.bins_p-self.options.bins_s+1)
        print('ROT-UKNOT'.ljust(fw),str(self.nbinsrot_uknot).rjust(fr),str(ROT_uknot_df).rjust(fr),
                str(np.nan).rjust(fr),str(np.nan).rjust(fr))

        if self.options.binsmethod=="IMSE direct plug-in":
            DPI_df = np.nan
            if not np.isnan(self.nbinsdpi):
                DPI_df = self.options.bins_p+1+(self.nbinsdpi-1)*(self.options.bins_p-self.options.bins_s+1)
            print('DPI'.ljust(fw),str(self.nbinsdpi).rjust(fr),str(DPI_df).rjust(fr),
                    ("{:.{}f}".format(self.imse_b_dpi, n_dec)).rjust(fr),("{:.{}f}".format(self.imse_v_dpi, n_dec)).rjust(fr))

            DPI_uknot_df = np.nan
            if not np.isnan(self.nbinsdpi_uknot):
                DPI_uknot_df = self.options.bins_p+1+(self.nbinsdpi_uknot-1)*(self.options.bins_p-self.options.bins_s+1)
            print('DPI-UKNOT'.ljust(fw),str(self.nbinsdpi_uknot).rjust(fr),str(DPI_uknot_df).rjust(fr),
                    str(np.nan).rjust(fr),str(np.nan).rjust(fr))
        print('-'*ft)
        return ''

class grid_str:
    def __init__(self, eval, bin, isknot, mid):
        self.eval = eval
        self.bin = bin
        self.isknot = isknot
        self.mid = mid

class data_str:
    def __init__(self, dots = None, line = None, ci = None, cb = None, poly = None, polyci= None):
        self.dots = dots
        self.line = line
        self.ci = ci
        self.cb = cb
        self.poly = poly
        self.polyci = polyci

class options:
    def __init__(self,dots,line,ci,cb,polyreg,deriv,binspos,
                 binsmethod, N_by, Ndist_by, Nclust_by,
                 nbins_by, byvals):
        self.dots=dots
        self.line=line
        self.ci=ci
        self.cb=cb
        self.polyreg=polyreg
        self.deriv=deriv
        self.binspos=binspos
        self.binsmethod=binsmethod
        self.N_by=N_by
        self.Ndist_by=Ndist_by
        self.Nclust_by=Nclust_by
        self.nbins_by=nbins_by
        self.byvals=byvals

class options_glm:
    def __init__(self,dots,line,ci,cb,polyreg,deriv, dist, link,
                 binspos,binsmethod, N_by, Ndist_by, Nclust_by,
                 nbins_by, byvals):
        self.dots=dots
        self.line=line
        self.ci=ci
        self.cb=cb
        self.polyreg=polyreg
        self.deriv=deriv
        self.dist=dist
        self.link=link
        self.binspos=binspos
        self.binsmethod=binsmethod
        self.N_by=N_by
        self.Ndist_by=Ndist_by
        self.Nclust_by=Nclust_by
        self.nbins_by=nbins_by
        self.byvals=byvals

class options_qreg:
    def __init__(self,dots,line,ci,cb,polyreg,deriv, quantile,
                 binspos,binsmethod, N_by, Ndist_by, Nclust_by,
                 nbins_by, byvals):
        self.dots=dots
        self.line=line
        self.ci=ci
        self.cb=cb
        self.polyreg=polyreg
        self.deriv=deriv
        self.quantile=quantile
        self.binspos=binspos
        self.binsmethod=binsmethod
        self.N_by=N_by
        self.Ndist_by=Ndist_by
        self.Nclust_by=Nclust_by
        self.nbins_by=nbins_by
        self.byvals=byvals

class binsreg_output:
    def __init__(self,bins_plot, data_plot, cval_by,options):
        self.bins_plot=bins_plot
        self.data_plot=data_plot
        self.cval_by=cval_by
        self.options = options

    def __repr__(self):
        print("Call: binsreg\n")
        fw = 35
        fr = 15
        print("Binscatter Plot")
        print("Bin selection method (binsmethod) = ", self.options.binsmethod.rjust(fr))
        print("Placement (binspos)               = ", str(self.options.binspos).rjust(fr))
        print("Derivative (deriv)                = ".ljust(fw), str(self.options.deriv).rjust(fr))
        print(" ")

        for i in range(len(self.options.byvals)):
            print("Group (by)                        = ", str(self.options.byvals[i]).rjust(fr))
            print("Sample size (n)                   = ", str(self.options.N_by[i]).rjust(fr))
            print("# of distinct values (Ndist)      = ", str(self.options.Ndist_by[i]).rjust(fr))
            print("# of clusters (Nclust)            = ", str(self.options.Nclust_by[i]).rjust(fr))
            print("dots, degree (p)                  = ", str(self.options.dots[0]).rjust(fr))
            print("dots, smooth (s)                  = ", str(self.options.dots[1]).rjust(fr))
            print("# of bins (nbins)                 = ", str(self.options.nbins_by[i]).rjust(fr))
            print("\n")
            
            fw = 10
            fr = 5
            print('='*(fw+4*fr))
            print(" ".ljust(fw),'p'.rjust(fr),'s'.rjust(fr), 'df'.rjust(fr))
            print('-'*(fw+4*fr))
            dots_df = self.options.dots[0]+1+(self.options.nbins_by[i]-1)*(self.options.dots[0]-self.options.dots[1]+1)
            print('dots'.rjust(fw),str(self.options.dots[0]).rjust(fr),str(self.options.dots[1]).rjust(fr),str(dots_df).rjust(fr))
            
            if self.options.line is not None:
                line_df = self.options.line[0]+1+(self.options.nbins_by[i]-1)*(self.options.line[0]-self.options.line[1]+1)
                print('line'.rjust(fw),str(self.options.line[0]).rjust(fr),str(self.options.line[1]).rjust(fr),str(line_df).rjust(fr))
            
            if self.options.ci is not None:
                ci_df = self.options.ci[0]+1+(self.options.nbins_by[i]-1)*(self.options.ci[0]-self.options.ci[1]+1)
                print('Ci'.rjust(fw),str(self.options.ci[0]).rjust(fr),str(self.options.ci[1]).rjust(fr),str(ci_df).rjust(fr))

            if self.options.cb is not None:
                cb_df = self.options.cb[0]+1+(self.options.nbins_by[i]-1)*(self.options.cb[0]-self.options.cb[1]+1)
                print('CB'.rjust(fw),str(self.options.cb[0]).rjust(fr),str(self.options.cb[1]).rjust(fr),str(cb_df).rjust(fr))

            if self.options.polyreg is not None:
                polyreg_df = self.options.polyreg + 1
                print('polyreg'.rjust(fw),str(self.options.polyreg).rjust(fr),'nan'.rjust(fr),str(polyreg_df).rjust(fr))
            print('-'*(fw+4*fr))
            return ''

class binsglm_output:
    def __init__(self,bins_plot, data_plot, cval_by, options):
        self.bins_plot=bins_plot
        self.data_plot=data_plot
        self.cval_by=cval_by
        self.options=options

    def __repr__(self):
        print("Call: binsglm\n")
        fw = 35
        fr = 15
        print("Binscatter Plot")
        print("Bin selection method (binsmethod) = ", self.options.binsmethod.rjust(fr))
        print("Placement (binspos)               = ", str(self.options.binspos).rjust(fr))
        print("Derivative (deriv)                = ", str(self.options.deriv).rjust(fr))
        print("Distribution Family               = ", self.options.dist.rjust(fr))
        print("Link Function                     = ", self.options.link.rjust(fr))
        print(" ")

        for i in range(len(self.options.byvals)):
            print("Group (by)                        = ", str(self.options.byvals[i]).rjust(fr))
            print("Sample size (n)                   = ", str(self.options.N_by[i]).rjust(fr))
            print("# of distinct values (Ndist)      = ", str(self.options.Ndist_by[i]).rjust(fr))
            print("# of clusters (Nclust)            = ", str(self.options.Nclust_by[i]).rjust(fr))
            print("dots, degree (p)                  = ", str(self.options.dots[0]).rjust(fr))
            print("dots, smooth (s)                  = ", str(self.options.dots[1]).rjust(fr))
            print("# of bins (nbins)                 = ", str(self.options.nbins_by[i]).rjust(fr))
            print("\n")
            
            fw = 10
            fr = 5
            print('='*(fw+4*fr))
            print(" ".ljust(fw),'p'.rjust(fr),'s'.rjust(fr), 'df'.rjust(fr))
            print('-'*(fw+4*fr))
            dots_df = self.options.dots[0]+1+(self.options.nbins_by[i]-1)*(self.options.dots[0]-self.options.dots[1]+1)
            print('dots'.rjust(fw),str(self.options.dots[0]).rjust(fr),str(self.options.dots[1]).rjust(fr),str(dots_df).rjust(fr))
            
            if self.options.line is not None:
                line_df = self.options.line[0]+1+(self.options.nbins_by[i]-1)*(self.options.line[0]-self.options.line[1]+1)
                print('line'.rjust(fw),str(self.options.line[0]).rjust(fr),str(self.options.line[1]).rjust(fr),str(line_df).rjust(fr))
            
            if self.options.ci is not None:
                ci_df = self.options.ci[0]+1+(self.options.nbins_by[i]-1)*(self.options.ci[0]-self.options.ci[1]+1)
                print('Ci'.rjust(fw),str(self.options.ci[0]).rjust(fr),str(self.options.ci[1]).rjust(fr),str(ci_df).rjust(fr))

            if self.options.cb is not None:
                cb_df = self.options.cb[0]+1+(self.options.nbins_by[i]-1)*(self.options.cb[0]-self.options.cb[1]+1)
                print('CB'.rjust(fw),str(self.options.cb[0]).rjust(fr),str(self.options.cb[1]).rjust(fr),str(cb_df).rjust(fr))

            if self.options.polyreg is not None:
                polyreg_df = self.options.polyreg + 1
                print('polyreg'.rjust(fw),str(self.options.polyreg).rjust(fr),'nan'.rjust(fr),str(polyreg_df).rjust(fr))
            print('-'*(fw+4*fr))
            return ''

class binsqreg_output:
    def __init__(self,bins_plot, data_plot, cval_by, options):
        self.bins_plot=bins_plot
        self.data_plot=data_plot
        self.cval_by=cval_by
        self.options=options

    def __repr__(self):
        print("Call: binsqreg\n")
        fw = 35
        fr = 15
        print("Binscatter Plot")
        print("Bin selection method (binsmethod) = ", self.options.binsmethod.rjust(fr))
        print("Placement (binspos)               = ", str(self.options.binspos).rjust(fr))
        print("Derivative (deriv)                = ", str(self.options.deriv).rjust(fr))
        print("Quantile                          = ", str(self.options.quantile).rjust(fr))
        print(" ")

        for i in range(len(self.options.byvals)):
            print("Group (by)                        = ", str(self.options.byvals[i]).rjust(fr))
            print("Sample size (n)                   = ", str(self.options.N_by[i]).rjust(fr))
            print("# of distinct values (Ndist)      = ", str(self.options.Ndist_by[i]).rjust(fr))
            print("# of clusters (Nclust)            = ", str(self.options.Nclust_by[i]).rjust(fr))
            print("dots, degree (p)                  = ", str(self.options.dots[0]).rjust(fr))
            print("dots, smooth (s)                  = ", str(self.options.dots[1]).rjust(fr))
            print("# of bins (nbins)                 = ", str(self.options.nbins_by[i]).rjust(fr))
            print("\n")
            
            fw = 10
            fr = 5
            print('='*(fw+4*fr))
            print(" ".ljust(fw),'p'.rjust(fr),'s'.rjust(fr), 'df'.rjust(fr))
            print('-'*(fw+4*fr))
            dots_df = self.options.dots[0]+1+(self.options.nbins_by[i]-1)*(self.options.dots[0]-self.options.dots[1]+1)
            print('dots'.rjust(fw),str(self.options.dots[0]).rjust(fr),str(self.options.dots[1]).rjust(fr),str(dots_df).rjust(fr))
            
            if self.options.line is not None:
                line_df = self.options.line[0]+1+(self.options.nbins_by[i]-1)*(self.options.line[0]-self.options.line[1]+1)
                print('line'.rjust(fw),str(self.options.line[0]).rjust(fr),str(self.options.line[1]).rjust(fr),str(line_df).rjust(fr))
            
            if self.options.ci is not None:
                ci_df = self.options.ci[0]+1+(self.options.nbins_by[i]-1)*(self.options.ci[0]-self.options.ci[1]+1)
                print('Ci'.rjust(fw),str(self.options.ci[0]).rjust(fr),str(self.options.ci[1]).rjust(fr),str(ci_df).rjust(fr))

            if self.options.cb is not None:
                cb_df = self.options.cb[0]+1+(self.options.nbins_by[i]-1)*(self.options.cb[0]-self.options.cb[1]+1)
                print('CB'.rjust(fw),str(self.options.cb[0]).rjust(fr),str(self.options.cb[1]).rjust(fr),str(cb_df).rjust(fr))

            if self.options.polyreg is not None:
                polyreg_df = self.options.polyreg + 1
                print('polyreg'.rjust(fw),str(self.options.polyreg).rjust(fr),'nan'.rjust(fr),str(polyreg_df).rjust(fr))
            print('-'*(fw+4*fr))
            return ''

#######################################################
########## Functions Definitons #######################
#######################################################

def factorial(n):
    return np.prod(range(1,n+1))

# n choose k function
def comb(n,k):
    return np.math.factorial(n)//np.math.factorial(k)//np.math.factorial(n-k)

def nanmat(n, m = None):
    # Create a (n x m) matrix of NaN
    if m is None: M = np.empty((n,))
    else: M = np.empty((n,m,))    
    M.fill(np.nan)
    return M

def ncol(x):
    try: d = x.shape[1]
    except: d = 1
    return d

def nrow(x):
    try: d = x.shape[0]
    except: d = 1
    return d

def cbind(x,y):
    if x is not None:
        if y is not None: return np.column_stack((x,y))
        else: return x
    else:
        if y is not None: return y
        else: return None

def rbind(x,y):
    if x is not None:
        if y is not None: return np.concatenate((x,y))
        else: return x
    else:
        if y is not None: return y
        else: return None

# generalized square root of a  matrix
def lssqrtm(A):
  u, s, vh =  np.linalg.svd(A)
  s[s<0] = 0
  return np.matmul(np.matmul(u, np.diag(np.sqrt(s))), vh)

def complete_cases(x):
    return np.all(np.invert(np.isnan(x)), axis = 1)

def genKnot_es(xmin, xmax, J):
    return np.linspace(xmin,xmax,J+1)

# Quantile knot list (including xmin and xmax as boundaries)
def genKnot_qs(x, J):
    return np.quantile(x, np.linspace(0,1,J+1),)

def FindInterval(x,knot):
    bin_index = np.searchsorted(knot, x,side = 'right')-1
    bin_index[x==knot[-1]]= len(knot)-2   # add the x's at the last knot to the last bin
    return bin_index

# Check local mass points
def binsreg_checklocalmass(x, J, es, knot = None):
    if knot is None:
        if es: knot = genKnot_es(np.min(x), np.max(x), J)
        else: knot = genKnot_qs(x, J)
    position = FindInterval(x,knot)
    n = [len(np.unique(x[position==i])) for i in range(J)]
    return np.min(n)

def colWeightedMeans(x, w = None):
    if w is not None: x = x*w
    return x.mean(0)

def colWeightedMedians(x, w = None):
    if w is not None: x = x*w
    return np.quantile(x,0.5,0)

# grid generation
def binsreg_grid(knot, ngrid, addmore = False):
    eval  = np.cumsum(np.append(knot[0],np.repeat(np.diff(knot)/(ngrid+1), ngrid+1)))
    eval = eval[1:-1]
    bin = isknot = mid = np.nan
    if addmore:
        bin = np.repeat(np.arange(len(knot)-1), ngrid+1)
        bin = bin[:-1]
        isknot = np.tile(np.append(np.repeat(0, ngrid), 1), len(knot)-1)
        isknot = isknot[:-1]
        id_mid = np.repeat(0, ngrid+1)
        id_mid[int(np.floor(ngrid/2))] = 1
        mid = np.tile(id_mid, len(knot)-1)
        mid = mid[:-1]
    return grid_str(eval, bin, isknot, mid)

# Wrapper for statsmodel.api
def binsreg_fit(y, x, weights = None, family = None, is_qreg = False, quantile = None,
                cov_type = None, cluster = None):
    
    aux0=aux1=aux2=aux3=aux4=''
    if cov_type=='cluster': 
        aux3 = 'cov_type = cov_type, cov_kwds={\'groups\': cluster}'
    elif cov_type is not None:
        aux3 = 'cov_type = cov_type'

    if is_qreg: 
        aux0 = 'QuantReg'
        if aux3=='':
            aux4 = 'q = quantile'
        else:
            aux4 = ', q = quantile'
    else:
        if family is None:
            if weights is None:
                aux0 = 'OLS'
            else:
                aux0 = 'WLS'
                aux2 = ', weights = weights'
        else:
            aux0 = 'GLM'
            aux1 = f', family = {family}'
            if weights is not None:
                aux2 = ', weights = weights'
    model = f'sm.{aux0}(y,x{aux1}{aux2}).fit({aux3}{aux4})'
    return eval(model)

# check drop, display warning
def check_drop(beta, k): 
    if any(np.isnan(beta[:k+1])): warnings.warn("some X-based variables dropped")

# Internal pred function (model is long regression, NA handled inside)
def binsreg_pred(X, model, type="xb", deriv=0, wvec=None, avar=False):
    k = ncol(X)
    fit = np.nan
    if type == "xb" or type == "all":
        coef = model.params
        coef[np.isnan(coef)] = 0
        if wvec is not None and deriv==0:
            fit = np.matmul(X,coef[:k]) + np.sum(wvec*coef[k:])
        else:
            fit = np.matmul(X,coef[:k])
    se = np.nan
    if type == "se" or type == "all":
        vcv = model.cov_params()
        if avar: 
            pos = np.invert(np.isnan(model.params[:k]))
            k_new = np.sum(pos)
            vcv = vcv[:k_new, :k_new]
        else:
            if wvec is not None:
                if deriv==0: X = np.column_stack((X, np.outer(np.ones(nrow(X)),wvec)))
                else: X = np.column_stack((X, np.outer(np.ones(nrow(X)), np.zeros(len(wvec)))))
            pos = np.invert(np.isnan(model.params))
        se = np.sqrt(np.sum(np.matmul(X[:,pos],vcv) * X[:,pos],1))
    return fit, se

# pval, cval simulation (tstat should be 2-col matrix)
def binsreg_pval(num, denom, rep, tstat=None, side=None, alpha=95, lp=np.inf):
    tvec = nanmat(rep)
    pval = np.nan
    if tstat is not None: pval = np.zeros(nrow(tstat))
    cval = np.nan
    k = ncol(num)

    for i in range(rep):
        eps = np.random.normal(size = k).reshape(-1,1)
        tx  = np.matmul(num, eps)/denom.reshape(-1,1)

        # for critical value
        if side is not None:
            if side == "two":
                if np.isinf(lp): tvec[i] = np.max(np.abs(tx))
                else: tvec[i] = (np.mean(np.abs(tx)**lp))**(1/lp)
            elif side == "left":
                tvec[i] = np.max(tx)
            elif side == "right":
                tvec[i] = np.min(tx)

        # for p value
        if tstat is not None:
            for j in range(nrow(tstat)):
                # 1: left; 2: right; 3: two-sided
                if tstat[j,1] == 1:
                    pval[j] += (np.max(tx) >= tstat[j,0])
                elif tstat[j,1] == 2:
                    pval[j] += (np.min(tx) <= tstat[j,0])
                elif tstat[j,1] == 3:
                    if np.isinf(lp): pval[j] = pval[j] + (np.max(np.abs(tx)) >= tstat[j,0])
                    else: pval[j] = pval[j] + (np.mean(np.abs(tx)**lp)**(1/lp) >= tstat[j,0])

    if tstat is not None: pval = pval / rep
    if side is not None: cval = np.quantile(tvec, alpha/100)
    return pval, cval

# pval used only by binspwc
def binspwc_pval(nummat1, nummat2, denom1, denom2, rep, tstat=None, testtype=None, lp=np.inf):
    pval = 0
    k1 = ncol(nummat1)
    k2 = ncol(nummat2)
    nummat1 = nummat1.reshape(-1,k1)
    nummat2 = nummat2.reshape(-1,k2)
    denom1 = denom1.reshape(len(denom1),-1)
    denom2 = denom2.reshape(len(denom2),-1)
    for i in range(rep):
        eps1 = np.random.normal(size = k1).reshape(-1,1)
        eps2 = np.random.normal(size = k2).reshape(-1,1)
        tx  = (np.matmul(nummat1, eps1) - np.matmul(nummat2, eps2))/np.sqrt(denom1**2+denom2**2)

        # for p value
        if testtype == "left":
            pval += (np.max(tx) >= tstat)
        elif testtype == "right":
            pval += (np.min(tx) <= tstat)
        else:
            if not np.isfinite(lp): pval += (np.max(np.abs(tx)) >= tstat)
            else: pval += (np.mean(np.abs(tx)**lp)**(1/lp) >= tstat)
    pval = pval / rep
    return pval

# Spline by recursion
def bspline(x, knot, p = 0, s = 0, deriv = 0):
    n = len(x)
    m = len(knot)
    if m-p<2:
        raise Exception("m<2+p is not allowed")
    B = np.zeros((n,m-1)) 
    B[range(n),FindInterval(x,knot)]=1 # degree = 0
    for k in range(1,p+1):
        aux = np.zeros((n,m-k-1))
        dem = knot[k:]-knot[:-k]
        for j in range(m-k-1):
            if dem[j]>0:
                aux[:,j] = ((x-knot[j])/dem[j])*B[:,j] 
            if dem[j+1]>0:
                aux[:,j] += (1- ((x-knot[j+1])/dem[j+1]))*B[:,j+1]
        B = aux
    return B

# Copy of the MATA code
def binsreg_spdes(x, p, s, knot, deriv):
    n = len(x)
    k = len(knot)
    #The  resulting spline is C**(s-1)
    bin_ind = FindInterval(x,knot).reshape(-1)

    # Extending the knots
    if len(knot) > 2:
        ext_knot = np.concatenate([np.repeat(knot[0],p+1), np.repeat(knot[1:-1], p-s+1), np.repeat(knot[-1], p+1)])
    else:
        ext_knot = np.concatenate([np.repeat(knot[0], p+1), np.repeat(knot[-1], p+1)])
    
    # index of the left knot for each x in the ext_knot 
    ind_lk = bin_ind*(p-s+1) + p
    # knot matrix: left and right endpoints of local supports
    lk = nanmat(n,p)
    rk = nanmat(n,p)
    for i in range(p):
        lk[:,p-1-i] = ext_knot[ind_lk -i]
        rk[:,i] = ext_knot[ind_lk+i+1]

    # Recursive Formula
    bs = np.ones((n,1))
    x = np.array(x).reshape(-1,1)
    if p > 0:
        if p < deriv:
            bs = np.zeros((n,p+1))
        elif p > deriv:	# loop: up to the degree of "degree-deriv"
            for i in range(p-deriv):
                tl = lk[:,-(i+1):]
                tr = rk[:,:(i+1)]
                w = (x-tl)/(tr-tl)
                bs = np.column_stack(((1-w)*bs,np.zeros(n))) + np.column_stack((np.zeros(n),w*bs))
            if deriv > 0: # loop: derivative
                for i in range(p-deriv,p):
                    tl = lk[:,-(i+1):]
                    tr = rk[:,:(i+1)]
                    w = 1/(tr-tl)
                    bs = (np.column_stack((np.zeros(n),w*bs)) - np.column_stack((w*bs,np.zeros(n))))*(i+1)
        else: # degree=deriv, so only need to account for derivative
            for i in range(p):
                tl = lk[:,-(i+1):]
                tr = rk[:,:(i+1)]
                w = 1/(tr-tl)
                bs = (np.column_stack((np.zeros(n),w*bs)) - np.column_stack((w*bs,np.zeros(n))))*(i+1)
    
    # bs should be an n by p+1 matrix containing nonzeros of the design matrix
    # The design should be a n by (k-2)*(p-s+1) + p+1 matrix 
    design = np.zeros((n,(k-2)*(p-s+1)+p+1))
    for i in range(k-1):
        r_ind = (bin_ind==i)
        c_ind = i*(p-s+1)
        design[r_ind,c_ind:(c_ind+p+1)] = bs[r_ind,:]

    return design

# slightly modified mean function
def binsreg_summ(x, w = None, std=False):
    sig = np.nan
    if w is None:
        mu = np.mean(x)
        if std: sig = np.std(x,ddof = 1)
    else:
        mu = np.average(x, weights=w)
        if std: sig = np.sqrt(np.sum(w*(x-mu)**2)/(np.sum(w)-1))
    return mu, sig

# IMSE V constant
def imse_vcons(p, deriv):
    n = p + 1
    V = 1/np.array([np.arange(j,n+j) for j in range(1,n+1)])
    Vderiv = np.zeros((n, n))
    for r in range(deriv,n):
        for c in range(deriv,n):
            Vderiv[r,c] = (1/(r + c + 1 - 2*deriv)
                            *(np.math.factorial(r)/np.math.factorial(r-deriv))
                            *(np.math.factorial(c)/np.math.factorial(c-deriv)))
    return np.sum(np.diag(np.linalg.solve(V, Vderiv)))


# IMSE B constant
def imse_bcons(p, deriv, s = 0):
    ord = p + 1
    if (s == 0):
        bcons = 1 / (2*(ord-deriv) + 1) / np.math.factorial(ord-deriv)**2 / comb(2*(ord-deriv), ord-deriv)**2
    else:
        if p==0: bernum = 1/6
        elif p==1: bernum = 1/30
        elif p==2: bernum = 1/42
        elif p==3: bernum = 1/30
        elif p==4: bernum = 5/66
        elif p==5: bernum = 691/2730
        elif p==6: bernum = 7/6
        else: raise Exception ("p>6 not allowed.")
        bcons = 1 / np.math.factorial(2*(ord-deriv)) * bernum
    return bcons

# ROT selector
def binsregselect_rot(y, x, w, p, s, deriv, eN, es=False, norotnorm=False, 
                        qrot=2, den_alpha=0.975, weights = None):
    x = (x - np.min(x)) / (np.max(x) - np.min(x))
    ord = p+1
    N = len(x)

    x_p = nanmat(N, p+qrot+1)
    for j in range(p+qrot+1):  x_p[:,j] = (x**j).reshape(-1)
    if w is not None: P = np.column_stack((x_p, w))
    else: P = x_p

    est = binsreg_fit(y, P, weights = weights)
    beta = est.params
    est = est.fittedvalues

    # variance constant
    s2 =  binsreg_fit(y**2, P, weights = weights).fittedvalues - est**2
    if norotnorm: fz = 1
    else:
        mu, sig = binsreg_summ(x, w=weights, std=True)
        fz = norm.pdf(x, loc = mu, scale = sig)
        # trim density from below
        cutval = norm.pdf(norm.ppf(den_alpha)*sig, loc = 0, scale = sig)
        fz[fz<cutval] = cutval
    if es: s2 = s2 / fz
    else: s2 = s2 * (fz**(2*deriv))
    s2 = binsreg_summ(s2, w = weights, std=False)[0]
    vcons = imse_vcons(p, deriv)
    imse_v = vcons * s2

    # bias constant
    bcons = imse_bcons(p, deriv = deriv, s=0)
    mu_m_hat = 0
    for j in range(p,p+qrot):
        mu_m_hat +=  x**(j-p)*beta[j+1]*np.math.factorial(j+1)/np.math.factorial(j-p)
    mu_m_hat = mu_m_hat**2
    if not es:
        mu_m_hat = mu_m_hat / (fz**(2*ord-2*deriv))
    imse_b = bcons * binsreg_summ(mu_m_hat, w=weights, std=False)[0]
    aux000num = imse_b*2*(ord-deriv)
    aux000dem = imse_v*(1+2*deriv)
    aux000 = aux000num/aux000dem
    aux00 = aux000**(1/(2*ord+1))
    aux0 = aux00 * eN**(1/(2*ord+1))
    aux1 = np.ceil(aux0)
    if not np.isnan(aux1): J_rot = int(aux1)
    else: J_rot = np.nan
    return J_rot, imse_b, imse_v

# locate h
def binsreg_locate(x, knot, type="all"):
    h = tl = np.nan
    loc_ind = FindInterval(x, knot)
    if type=="all" or type=="h":
        size    = np.diff(knot)
        h       = size[loc_ind]
    if type=="all" or type=="tl":
        pos     = knot[:-1]
        tl      = pos[loc_ind]
    return h, tl

def bernpoly(x, p):
    n = len(x)
    if p==0: bernx = np.ones(n)
    elif p==1: bernx = x - 0.5
    elif p==2: bernx = x**2-x+1/6
    elif p==3: bernx = x**3 - 1.5*x**2 + 0.5*x
    elif p==4: bernx = x**4 - 2*x**3 + x**2 - 1/30
    elif p==5: bernx = x**5 - 2.5*x**4 + 5/3*x**3 - 1/6*x
    elif p==6: bernx = x**6 - 3*x**5 + 2.5*x**4 - 0.5*x**2 + 1/42
    else: raise Exception("p>6 not allowed.")
    return bernx

# bias term
def bias(x, p, s, deriv, knot):
    h, tl = binsreg_locate(x, knot)
    bern = bernpoly((x-tl)/h, p+1-deriv) / np.math.factorial(p+1-deriv) * (h**(p+1-deriv))
    return bern

# IMSE V cons
def genV(y, x, w, p, s, deriv, knot, vce, cluster=None, weights=None):
    B  = binsreg_spdes(x=x, p=p, s=s, knot=knot, deriv=0)
    k  = ncol(B)
    if deriv>0: basis = binsreg_spdes(x=x, p=p, s=s, knot=knot, deriv=deriv)
    else: basis = B
    if w is not None: P = np.column_stack((B, w))
    else: P = B
    model  = binsreg_fit(y, P, weights = weights, cov_type=vce, cluster=cluster)
    pos = np.invert(np.isnan(model.params[:k]))
    k_new = np.sum(pos)
    vcv = model.cov_params()[:k_new,:k_new]
    m_s2 = binsreg_summ(np.sum(np.matmul(basis[:,pos], vcv) * basis[:,pos],1), w=weights, std=False)[0]
    return m_s2


def genB(y, x, w, p, s, deriv, knot, weights=None):
    B  = binsreg_spdes(x=x, p=p+1, s=s+1, knot=knot, deriv=0)  # use smoothest spline
    k  = ncol(B)
    if w is not None: P = np.column_stack((B, w))
    else: P = B
    beta = binsreg_fit(y, P, weights = weights).params[:k]
    pos  = np.invert(np.isnan(beta))
    basis = binsreg_spdes(x=x, p=p+1, s=s+1, knot=knot, deriv=p+1)
    mu_m_fit  = (np.matmul(basis[:,pos],beta[pos])).reshape(-1,1)

    bias_0 = mu_m_fit * bias(x, p, s, 0, knot)    # proj component, v=0!!!
    B = binsreg_spdes(x=x, p=p, s=s, knot=knot, deriv=0)
    beta = binsreg_fit(bias_0, B, weights = weights).params
    pos = np.invert(np.isnan(beta))
    if deriv > 0:
        basis = binsreg_spdes(x=x, p=p, s=s, knot=knot, deriv=deriv)
        bias_v = mu_m_fit * bias(x, p, s, deriv, knot)    # need to recalculate for v>0!!!
    else:
        basis = B
        bias_v = bias_0
    bias_l2 = bias_v - (np.matmul(basis[:,pos], beta[pos])).reshape(-1,1)
    bias_cons = binsreg_summ(bias_l2**2, w=weights, std=False)[0]

    return bias_cons


# DPI selector
def binsregselect_dpi(y, x, w, p, s, deriv, vce, nbinsrot, es=False, cluster=None, weights=None):
    J_rot = nbinsrot
    x = (x - np.min(x)) / (np.max(x) - np.min(x))
    ord = p+1

    if es: knot = genKnot_es(0, 1, J_rot)
    else: knot = genKnot_qs(x, J_rot)
    
    # bias constant
    imse_b = genB(y, x, w, p, s, deriv, knot, weights=weights) * J_rot**(2*(ord-deriv))

    # variance constant
    imse_v = genV(y, x, w, p, s, deriv, knot, vce, cluster, weights=weights) / (J_rot**(1+2*deriv))
    J_dpi = int(np.ceil((imse_b*2*(ord-deriv)/((1+2*deriv)*imse_v))**(1/(2*ord+1))))

    return J_dpi, imse_v, imse_b
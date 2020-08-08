/*******************************************************************************
BINSCATTER
Date: 13-MAR-2019 
Authors: Matias Cattaneo, Richard K. Crump, Max H. Farrell, Yingjie Feng
*******************************************************************************/
** hlp2winpdf, cdn(binsreg) replace
** hlp2winpdf, cdn(binsregtest) replace
** hlp2winpdf, cdn(binsregselect) replace
*******************************************************************************/
** net install binsreg, from(https://sites.google.com/site/nppackages/binsreg/stata) replace
********************************************************************************
clear all
set more off
set linesize 90

********************************************************************************
** Binsreg Sim Data: Generation / COMMENTED OUT
********************************************************************************
/*

local n=1000
set obs `n'
set seed 1234

* X-var
gen x=runiform()
* True fun
gen mu=4*x-2+0.8*exp(-256*(x-0.5)^2)
* A continuous covariate
gen w=runiform()*2-1
* A binary covariate
gen t=(runiform()>0.5)
* error term
gen eps=rnormal()
* A cluster id, only for illustration purpose
gen id=ceil(_n/2)
* Y-var
gen y=mu+w+t+eps

keep y x w t id
save binsreg_simdata, replace version(13)

*/

********************************************************************************
** Binscatter Sim Data: Setup and summary stats
********************************************************************************
sjlog using output/binsreg_out0, replace
use binsreg_simdata, clear
sum
sjlog close, replace

********************************************************************************
** BINSREG
********************************************************************************
* Default syntax
sjlog using output/binsreg_out1, replace
binsreg y x w
sjlog close, replace

graph export output/binsreg_fig1.pdf, replace

* Setting quantile-spaced bins to J=20, add a linear fit
sjlog using output/binsreg_out2, replace
binsreg y x w, nbins(20) polyreg(1)
sjlog close, replace

* Adding lines, ci, cb, polyreg
sjlog using output/binsreg_out3, replace
qui binsreg y x w, nbins(20) dots(0,0) line(3,3)
qui binsreg y x w, nbins(20) dots(0,0) line(3,3) ci(3,3)
qui binsreg y x w, nbins(20) dots(0,0) line(3,3) ci(3,3) cb(3,3)
qui binsreg y x w, nbins(20) dots(0,0) line(3,3) ci(3,3) cb(3,3) polyreg(4)
sjlog close, replace

binsreg y x w, nbins(20) dots(0,0) line(3,3)
graph export output/binsreg_fig2a.pdf, replace
binsreg y x w, nbins(20) dots(0,0) line(3,3) ci(3,3)
graph export output/binsreg_fig2b.pdf, replace
binsreg y x w, nbins(20) dots(0,0) line(3,3) ci(3,3) cb(3,3)
graph export output/binsreg_fig2c.pdf, replace
binsreg y x w, nbins(20) dots(0,0) line(3,3) ci(3,3) cb(3,3) polyreg(4)
graph export output/binsreg_fig2d.pdf, replace

* VCE option, factor vars, twoway options, graph data saving  
sjlog using output/binsreg_out4, replace
binsreg y x w i.t, dots(0,0) line(3,3) ci(3,3) cb(3,3) polyreg(4) ///
                   vce(cluster id) savedata(output/graphdat) replace ///
				   title("Binned Scatter Plot") 
sjlog close, replace

* Comparision by groups
sjlog using output/binsreg_out5, replace
binsreg y x w, by(t) dots(0,0) line(3,3) cb(3,3) ///
               bycolors(blue red) bysymbols(O T) 
sjlog close, replace
graph export output/binsreg_fig3.pdf, replace

********************************************************************************
** BINSREGTEST
********************************************************************************
* Basic syntax: linear?
sjlog using output/binsreg_out6, replace
binsregtest y x w, testmodelpoly(1)
sjlog close, replace

* Alternative: save parametric fit in another file
* If not available, first create empty file with grid points using binsregselect
sjlog using output/binsreg_out7, replace
qui binsregselect y x w, simsgrid(30) savegrid(output/parfitval) replace
qui reg y x w
use output/parfitval, clear
predict binsreg_fit_lm
save output/parfitval, replace
use binsreg_simdata, clear
binsregtest y x w, testmodelparfit(output/parfitval)
sjlog close, replace

* Shape restriction test: increasing?
sjlog using output/binsreg_out8, replace
binsregtest y x w, deriv(1) nbins(20) testshaper(0)
sjlog close, replace

* Test many things simultaneously
sjlog using output/binsreg_out9, replace
binsregtest y x w, nbins(20) testshaper(-2 0) testshapel(4) testmodelpoly(1) ///
                   nsims(1000) simsgrid(30)
sjlog close, replace


********************************************************************************
** BINSREGSELECT
********************************************************************************
* Basic syntax
sjlog using output/binsreg_out10, replace
binsregselect y x w
sjlog close, replace

* J ROT specified manually and require evenly-spaced binning
sjlog using output/binsreg_out11, replace
binsregselect y x w, nbinsrot(20) binspos(es)
sjlog close, replace

* Save grid for prediction purpose
sjlog using output/binsreg_out12, replace
binsregselect y x w, simsgrid(30) savegrid(output/parfitval) replace
sjlog close, replace

* Extrapolating the optimal number of bins to the full sample
sjlog using output/binsreg_out13, replace
binsregselect y x w if t==0, useeffn(1000)
sjlog close, replace

*******************************************************************************
** BINSREG, integrated with BINSREGTEST and BINSREGSELECT
*******************************************************************************
* Shut down checks to speed computation
sjlog using output/binsreg_out14, replace
qui binsreg y x w, masspoints(off)
sjlog close, replace

* Automatic bin selection, binscatter plotting, and testing
sjlog using output/binsreg_out15, replace
binsreg y x w, dots(0,0) line(3,3) ci(3,3) cb(3,3) polyreg(4) ///
               testmodelpoly(1) testshapel(4)
sjlog close, replace


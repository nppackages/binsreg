/*******************************************************************************
Binsreg: illustration file for plot
# Authors: M. D. Cattaneo, R. Crump, M. Farrell, Y. Feng and Ricardo Masini
# Last update: OCT 9, 2022 
*******************************************************************************/
*******************************************************************************/
** net install binsreg, from(https://raw.githubusercontent.com/nppackages/binsreg/master/stata) replace
********************************************************************************
clear all
set more off
set linesize 90

********************************************************************************
** Binscatter Sim Data: Setup and summary stats
********************************************************************************
use binsreg_simdata, clear
sum

********************************************************************************
** EXTERNAL PLOT USING BINSREG OUTPUT 
********************************************************************************

* Run binsreg, save data in the file "result.dta" 
binsreg y x w, line(3,3) ci(3,3) cb(3,3) polyreg(4) savedata(result)

* Read the data produced by binsreg 
use result, clear

* Add the dots
twoway (scatter dots_fit dots_x, sort mcolor(blue) msymbol(smcircle)), ytitle(Y) xtitle(X)

* Add the line
twoway (scatter dots_fit dots_x, sort mcolor(blue) msymbol(smcircle)) ///
       (line line_fit line_x, sort lcolor(blue)), ytitle(Y) xtitle(X) legend(off)

* Add the CI
twoway (scatter dots_fit dots_x, sort mcolor(blue) msymbol(smcircle)) ///
       (line line_fit line_x, sort lcolor(blue)) ///
	   (rcap CI_l CI_r CI_x, sort lcolor(blue)), ytitle(Y) xtitle(X) legend(off)

* Add the CB
twoway (scatter dots_fit dots_x, sort mcolor(blue) msymbol(smcircle)) ///
       (line line_fit line_x, sort lcolor(blue)) ///
	   (rcap CI_l CI_r CI_x, sort lcolor(blue)) ///
	   (rarea CB_l CB_r CB_x, sort lcolor(none%0) fcolor(blue%50) fintensity(50)), ytitle(Y) xtitle(X) legend(off)
	   
* Add the polyreg
twoway (scatter dots_fit dots_x, sort mcolor(blue) msymbol(smcircle)) ///
       (line line_fit line_x, sort lcolor(blue)) ///
	   (rcap CI_l CI_r CI_x, sort lcolor(blue)) ///
	   (rarea CB_l CB_r CB_x, sort lcolor(none%0) fcolor(blue%50) fintensity(50)) ///
	   (line poly_fit poly_x, sort lcolor(red)), ytitle(Y) xtitle(X) legend(off)

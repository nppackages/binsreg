*! version 0.2 13-MAR-2019

capture program drop binsregtest
program define binsregtest, eclass
    version 13
	 
	syntax varlist(min=2 numeric fv ts) [if] [in] [fw aw pw] [, deriv(integer 0) ///
		   testmodel(numlist integer max=2 >=0)  ///
		   testmodelparfit(string asis) testmodelpoly(string) ///
		   testshape(numlist integer max=2 >=0) ///
		   testshapel(numlist) testshaper(numlist) testshape2(numlist) ///
		   bins(numlist integer max=2 >=0) ///
		   nbins(integer 0) binspos(string) binsmethod(string) nbinsrot(string) ///
		   nsims(integer 500) simsgrid(integer 20) simsseed(integer 666) ///
		   dfcheck(numlist integer max=2 >=0) masspoints(string) ///
		   vce(passthru) ///
		   numdist(string) numclust(string)]
		   /* last line only for internal use */
	
	 * Regularization constant (for checking only)
	 local qrot=2
	 
	 **************************************
	 * Create weight local
     if ("`weight'"!="") {
	    local wt [`weight'`exp']
		local wtype=substr("`weight'",1,1)
	 }
	 
	 * Extract options
	 * default vce
	 if ("`vce'"=="") local vce "vce(robust)"
	 
	 * binning
	 tokenize `bins'
	 local binsp "`1'"
	 local binss "`2'"
	 if ("`binsp'"=="") local binsp=0
	 if ("`binss'"=="") local binss=`binsp'
	 if (`nbins'!=0) local binselectmethod "User-specified"
	 if ("`binspos'"=="es") local binspos "ES"
	 if ("`binspos'"=="qs") local binspos "QS"
	 if ("`binspos'"=="")   local binspos "QS"
	 if ("`binsmethod'"=="rot") local binsmethod "ROT"
	 if ("`binsmethod'"=="dpi") local binsmethod "DPI"
	 if ("`binsmethod'"=="")    local binsmethod "DPI"
	 
	 * Option for testing shape
	 tokenize `testshape'
	 local tsha_p "`1'"
	 local tsha_s "`2'"
	 if ("`tsha_p'"=="") local tsha_p 3
	 if ("`tsha_s'"=="") local tsha_s `tsha_p'
	 local val_L `testshapel'
	 local nL: word count `val_L'
	 local val_R `testshaper'
	 local nR: word count `val_R'
	 local val_T `testshape2'
	 local nT: word count `val_T'
	 local ntestshape=`nL'+`nR'+`nT'     /* number of tests (for shape) */
	 
	 * Option for testing model
	 if ("`testmodelpoly'"!="") {
	    confirm integer n `testmodelpoly'
		local testpolyp=`testmodelpoly'
	 }
	 tokenize `testmodel'
	 local tmod_p "`1'"
	 local tmod_s "`2'"
	 if ("`tmod_p'"=="") local tmod_p 3
	 if ("`tmod_s'"=="") local tmod_s `tmod_p'
	 
	 * mass check?
	 if ("`masspoints'"=="") {
	    local massadj "T"
		local localcheck "T"	    
	 }
	 else if ("`masspoints'"=="off") {
	    local massadj "F"
		local localcheck "F"
	 }
	 else if ("`masspoints'"=="noadjust") {
	    local massadj "F"
		local localcheck "T"
	 }
	 else if ("`masspoints'"=="nolocalcheck") {
	 	local massadj "T"
		local localcheck "F"
     }
	 else if ("`masspoints'"=="veryfew") {
	    di as error "veryfew() not allowed for testing."
		exit
	 }
	 
	 * extract dfcheck
	 if ("`dfcheck'"=="") local dfcheck 20 30
	 tokenize `dfcheck'
	 local dfcheck_n1 "`1'"
	 local dfcheck_n2 "`2'"
	 
	 * Error check
	 if (`"`testmodelparfit'"'==`""'&`ntestshape'==0&"`testmodelpoly'"=="") {
	     di as error "no tests specified."
		 exit
	 }
	 if (`tsha_p'<`tsha_s'|`tmod_p'<`tmod_s'|`binsp'<`binss') {
	     di as error "p cannot be smaller than s."
		 exit
	 }
	 if (`tsha_p'<=`binsp'|`tmod_p'<=`binsp') {
	     di as text in gr "warning: p for testing > p for bins() suggested."
	 }
	 if (`tsha_p'<`deriv'|`tmod_p'<`deriv') {
	    di as error "p for test cannot be smaller than deriv."
		exit
	 }
	 if ("`testmodelpoly'"!="") {
		if (`testpolyp'<`deriv') {
		   di as error "degree of polynomial model cannot be smaller than deriv."
		   exit
		}
	 }
	 
	 * Mark sample
	 preserve
	 marksample touse, nov         /* do not account for missing values !! */
	 qui keep if `touse'
	 
	 * Parse varlist into y_var, x_var and w_var
	 tokenize `varlist'
	 fvrevar `1', tsonly
	 local y_var "`r(varlist)'"
	 fvrevar `2', tsonly
	 local x_var "`r(varlist)'"
	 
	 macro shift 2
	 local w_var "`*'"
	 fvrevar `w_var', tsonly
	 local w_var "`r(varlist)'"
	 
	 marksample touse      /* now renew the mark to account for missing values */
	 qui keep if `touse'
     local eN=_N
	 local nsize=_N     /* # of rows in the original dataset */
	     
	 if ("`masspoints'"!="off"|"`binspos'"=="QS") {
	    if ("`:sortedby'"!="`x_var'") sort `x_var'
	 }
	 
	 if ("`wtype'"=="f") qui sum `x_var' `wt', meanonly
	 else                qui sum `x_var', meanonly
	 
	 local xmin=r(min)
	 local xmax=r(max)
	 local N=r(N)        /* sample size, with wt */
	 
	 tempname xvec binedges
	 mata: `xvec'=st_data(., "`x_var'")
	 * effective sample size
     local Ndist=.
	 if ("`massadj'"=="T") {
	    if ("`numdist'"!=""&"`numdist'"!=".") {
		   local Ndist=`numdist'
		}
		else {
		   mata: `binedges'=binsreg_uniq(`xvec', ., 1, "Ndist")
		   mata: mata drop `binedges'
		}   
		local eN=min(`eN', `Ndist')
	 }
	 local Nclust=.
	 if ("`vce'"!="") {
	    local vcetemp: subinstr local vce "vce(" "", all
		local vcetemp: subinstr local vcetemp ")" "", all
		tokenize `vcetemp'
		if ("`1'"=="cl"|"`1'"=="clu"|"`1'"=="clus"|"`1'"=="clust"| /// 
		    "`1'"=="cluste"|"`1'"=="cluster") {
		    if ("`numclust'"!=""&"`numclust'"!=".") {
			   local Nclust=`numclust'
			}
			else {
			   mata: st_local("Nclust", strofreal(rows(uniqrows(st_data(.,"`2'")))))
			}
		}
		local eN=min(`eN', `Nclust')
	 }
	 
	 **********************************
	 ********** Bins ******************
	 **********************************
	 * determine # of bins
	 if ("`binspos'"!="ES"&"`binspos'"!="QS") {
	   capture numlist "`binspos'", ascending
	   if (_rc==0) {
		  local knotlist `binspos'
		  local nbins: word count `knotlist'
		  local first: word 1 of `knotlist'
		  local last: word `nbins' of `knotlist'
		  if (`first'<`xmin'|`last'>`xmax') {
			 di as error "inner knots specified out of allowed range."
			 exit
		  }
		  else {
		     local nbins=`nbins'+1
			 local binspos "user"
		  }
	   }
	   else {
		  di as error "numeric list incorrectly specified in binspos()."
		  exit
	   }
	 }
	 
	 * if binsmethod is specified
	 if (`nbins'==0) {
		 * Check effective sample size
		 if ("`nbinsrot'"==""&(`eN'<=`dfcheck_n1'+`binsp'+1+`qrot')) {
		    * ROT inavailable, exit
			di as error "too few observations for bin selection."
	        exit
	     }
		 else {
			qui binsregselect `y_var' `x_var' `w_var' `wt', deriv(`deriv') bins(`binsp' `binss') ///
		                      binsmethod(`binsmethod') binspos(`binspos') nbinsrot(`nbinsrot') ///
							  `vce' masspoints(`masspoints') dfcheck(`dfcheck_n1' `dfcheck_n2') ///
							  numdist(`Ndist') numclust(`Nclust')
			if (e(nbinsrot_regul)==.) {
			   di as error "bin selection fails."
			   exit
			}
			if ("`binsmethod'"=="ROT") {
		       local nbins=e(nbinsrot_regul)
		    }
		    else if ("`binsmethod'"=="DPI") {
			   local nbins=e(nbinsdpi)
			   if (`nbins'==.) {
			      local nbins=e(nbinsrot_regul)
			      di as text in gr "warning: DPI selection fails. ROT choice used."
			   }
		    }
		 }
	 }
	 
	 *******************************************************
	 * Check if eff. sample size is large enough for testing
	 if ((`nbins'-1)*(`tsha_p'-`tsha_s'+1)+`tsha_p'+1+`dfcheck_n2'>=`eN') {
	    local tsha_fewobs "T"
		di as text in gr "warning: too small effective sample size for testing shape."
	 }
	 if ((`nbins'-1)*(`tmod_p'-`tmod_s'+1)+`tmod_p'+1+`dfcheck_n2'>=`eN') {
	    local tmod_fewobs "T"
		di as text "warning: too small effective sample size for testing models."
	 }
	 ********************************************************
	 	 	 
	 * Generate category variable for data and save knot in matrix
	 tempname kmat
	 tempvar xcat
	 qui gen `xcat'=. in 1
		
     if ("`pos'"=="ES") {
	    local stepsize=(`xmax'-`xmin')/`nbins'
	    forvalues i=1/`=`nbins'+1' {
		   mat `kmat'=(nullmat(`kmat') \ `=`xmin'+`stepsize'*(`i'-1)')
		}
	 }
	 else if ("`knotlist'"!="") {
	    foreach el of local knotlist {
		   mat `kmat'=(nullmat(`kmat') \ `el')
		}
        mat `kmat'=(`xmin' \ `kmat' \ `xmax')
	 }
	 else {
		if (`nbins'==1)  mat `kmat'=(`xmin' \ `xmax')
		else {		
	       binsreg_pctile `x_var' `wt', nq(`nbins')
		   mat `kmat'=(`xmin' \ r(Q) \ `xmax')
		}
	 }
	 
	 mata: st_matrix("`kmat'", (`xmin' \ uniqrows(st_matrix("`kmat'")[|2 \ `=`nbins'+1'|])))
	 binsreg_irecode `x_var', knotmat(`kmat') bin(`xcat')
	 if (`nbins'!=rowsof(`kmat')-1) {
	     di as text in gr "warnings: repeated knots. Some bins dropped."
		 local nbins=rowsof(`kmat')-1
	 }
	 
	 * Check for empty bins
	 if ("`localcheck'"=="T") {
	   mata: `binedges'=binsreg_uniq(`xvec', st_data(.,"`xcat'"), `nbins', "uniqmin")
	   mata: mata drop `binedges'
	   if (`ntestshape'!=0) {
		  if (`uniqmin'<`tsha_p'+1) {
		     local tsha_fewobs "T"
		     di as text in gr "warning: some bins have too few distinct x-values for testing."
		  }
	   }
	   if (`"`testmodelparfit'"'!=`""'|"`testmodelpoly'"!="") {
		  if (`uniqmin'<`tmod_p'+1) {
		     local tmod_fewobs "T"
		     di as text in gr "warning: some bins have too few distinct x-values for testing."
		  }
	   }
	 }
	 
	 ********************************************************
	 * Set seed
	 set seed `simsseed'
	 local uni_last=`simsgrid'*`nbins'+`nbins'-1
	 
	 tempname Xm uni_grid uni_basis tstat            /* objects in MATA */
	 mata: `uni_grid'=binsreg_grids("`kmat'", `simsgrid')
	 
	 *******************************
	 ******* Testing Shape *********
	 *******************************
	 tempname stat_shape pval_shape       /* test stat and p value */
	 if (`ntestshape'!=0&"`tsha_fewobs'"!="T") {
	    * Regression
	    local nseries=(`tsha_p'-`tsha_s'+1)*(`nbins'-1)+`tsha_p'+1
	    local tsha_series ""
	    forvalues i=1/`nseries' {
	       tempvar sp`i'
	       local tsha_series `tsha_series' `sp`i''
		   qui gen `sp`i''=. in 1
	    }
		
		tempname tsha_b tsha_V
		mata: binsreg_st_spdes(`xvec', "`tsha_series'", "`kmat'", st_data(.,"`xcat'"), `tsha_p', 0, `tsha_s')
	    capture reg `y_var' `tsha_series' `w_var' `wt', nocon `vce' 
	    * store results
	    if (_rc==0) {
		    matrix `tsha_b'=e(b)
		    matrix `tsha_V'=e(V)
			mata: binsreg_checkdrop("`tsha_b'", "`tsha_V'", `nseries')
	    }
	    else {
	        error  _rc
	   	    exit _rc
        }
		   
		* Predict		
	    * fitted values
  	    mata: `uni_basis'=binsreg_spdes(`uni_grid'[,1], "`kmat'", `uni_grid'[,3], `tsha_p', `deriv', `tsha_s'); ///
	          `Xm'=binsreg_pred(`uni_basis', st_matrix("`tsha_b'")[|1 \ `nseries'|]', ///
		                        st_matrix("`tsha_V'")[|1,1 \ `nseries',`nseries'|], "all")
		
		mata: `tstat'=J(`ntestshape',2,.)
		  
		forval i=1/`ntestshape' {
		   if (`i'<=`nL') {
			  local val: word `i' of `val_L'
			  mata: `tstat'[`i',.]=(max((`Xm'[,1]:-`val'):/`Xm'[,2]), 1)			
		   }
		   else if (`i'<=`nL'+`nR') {
			  local val: word `=`i'-`nL'' of `val_R'
			  mata: `tstat'[`i',.]=(min((`Xm'[,1]:-`val'):/`Xm'[,2]), 2)
		   }
		   else {
			  local val: word `=`i'-`nL'-`nR'' of `val_T'
			  mata: `tstat'[`i',.]=(max(abs((`Xm'[,1]:-`val'):/`Xm'[,2])), 3)
		   }
		}
		mata: st_matrix("`stat_shape'", `tstat')
		
		* p value
		mata: binsreg_pval(`uni_basis', `Xm'[,2], "`tsha_V'", "`stat_shape'", `nsims', `nseries', ///
		                   ".", 0, "`pval_shape'", ".")
		drop `tsha_series'
		mata: mata drop `Xm' `uni_basis' `tstat'
		
		if ("`testshapel'"!="") {
	 	   tempname stat_shapeL pval_shapeL
	       mat `stat_shapeL'=`stat_shape'[1..`nL',1]
		   mat `pval_shapeL'=`pval_shape'[1..`nL',1]
		}
		if ("`testshaper'"!="") {
		   tempname stat_shapeR pval_shapeR
	       mat `stat_shapeR'=`stat_shape'[`=`nL'+1'..`=`nL'+`nR'',1]
		   mat `pval_shapeR'=`pval_shape'[`=`nL'+1'..`=`nL'+`nR'',1]
		}
		if ("`testshape2'"!="") {
		   tempname stat_shape2 pval_shape2
	       mat `stat_shape2'=`stat_shape'[`=`nL'+`nR'+1'..`ntestshape',1]
		   mat `pval_shape2'=`pval_shape'[`=`nL'+`nR'+1'..`ntestshape',1]
	    }
	 }
	 else {
	    local tsha_p=.
		local tsha_s=.
	 }
	 
	 *************************************
	 ****** Testing Models ***************
	 *************************************
	 tempname stat_poly pval_poly               /* for testing poly reg */
	 tempname stat_model pval_model             /* for testing models   */
	 if ((`"`testmodelparfit'"'!=`""'|"`testmodelpoly'"!="")&"`tmod_fewobs'"!="T") {
	    ***********************************************
	 	* Regression: for BOTH
	    local nseries=(`tmod_p'-`tmod_s'+1)*(`nbins'-1)+`tmod_p'+1
	    local tmod_series ""
	    forvalues i=1/`nseries' {
	       tempvar sp`i'
	       local tmod_series `tmod_series' `sp`i''
		   qui gen `sp`i''=. in 1
	    }
		
		tempname tmod_b tmod_V
		capture confirm matrix `tsha_b' `tsha_V'
		if (_rc==0&`tmod_p'==`tsha_p'& `tmod_s'==`tsha_s') {
		    matrix `tmod_b'=`tsha_b'
		    matrix `tmod_V'=`tsha_V'
		} 
		else {
		   mata: binsreg_st_spdes(`xvec', "`tmod_series'", "`kmat'", st_data(.,"`xcat'"), `tmod_p', 0, `tmod_s')
	       capture reg `y_var' `tmod_series' `w_var' `wt', nocon `vce'
	 
	       * store results
	       if (_rc==0) {
		       matrix `tmod_b'=e(b)
		       matrix `tmod_V'=e(V)
			   mata: binsreg_checkdrop("`tmod_b'", "`tmod_V'", `nseries')
	       }
	       else {
	           error  _rc
	   	       exit _rc
           }
		}
		qui drop `tmod_series'
		
		********************************************************
		* If a test for poly reg is requested
		if ("`testmodelpoly'"!="") {
		   	* fitted values
  	       mata: `uni_basis'=binsreg_spdes(`uni_grid'[,1], "`kmat'", `uni_grid'[,3], `tmod_p', `deriv', `tmod_s'); ///
	             `Xm'=binsreg_pred(`uni_basis', st_matrix("`tmod_b'")[|1 \ `nseries'|]', ///
			                    st_matrix("`tmod_V'")[|1,1 \ `nseries',`nseries'|], "all")
		   
		   tempvar poly_fit
	       local poly_series ""
	       forval i=0/`testpolyp' {
		      tempvar x_var_`i'
			  qui gen `x_var_`i''=`x_var'^`i'
	          local poly_series `poly_series' `x_var_`i''
		   }
		 
		   capture reg `y_var' `poly_series' `w_var' `wt', nocon
		   * store results
		   tempname poly_b
	       if (_rc==0) {
	 	       matrix `poly_b'=e(b)
			   matrix `poly_b'=`poly_b'[1, `=`deriv'+1'..`=`testpolyp'+1']
	       }
	       else {
	           error  _rc
	   	       exit _rc
           }
		 
		   * Data for derivative
		   tempname polym
		   mata: `polym'=J(`uni_last',0,.)
		   forval i=`deriv'/`testpolyp' {
			   mata: `polym'=(`polym', `uni_grid'[,1]:^(`i'-`deriv')*factorial(`i')/factorial(`i'-`deriv'))
	       }
		   
		   mata: `polym'=`polym'*st_matrix("`poly_b'")'; ///
	   	         st_matrix("`stat_poly'", (max(abs((`Xm'[,1]-`polym'):/`Xm'[,2])),3)); ///
			     binsreg_pval(`uni_basis', `Xm'[,2], "`tmod_V'", "`stat_poly'", ///
				              `nsims', `nseries', ".", 0, "`pval_poly'", ".")
							  
		   mata: mata drop `Xm' `polym' `uni_basis'
		}
		
		******************************************************************
		* if the model is stored in another file
		if (`"`testmodelparfit'"'!=`""') {
	       use `"`testmodelparfit'"', clear
		   qui ds binsreg_fit*
		   local varls=r(varlist)
		   local nfitval: word count `varls'
		   tempvar  uni_xcat uni_fit uni_se

	       qui gen `uni_fit'=. in 1
	       qui gen `uni_se'=. in 1
		   qui gen `uni_xcat'=. in 1
		   binsreg_irecode `x_var', knotmat(`kmat') bin(`uni_xcat')
		
		   mata: `uni_basis'=binsreg_spdes(st_data(.,"`x_var'"), "`kmat'", st_data(.,"`uni_xcat'"), ///
		                                  `tmod_p', `deriv', `tmod_s'); ///
	             `Xm'=binsreg_pred(`uni_basis', st_matrix("`tmod_b'")[|1 \ `nseries'|]', ///
				               st_matrix("`tmod_V'")[|1,1 \ `nseries',`nseries'|], "all"); ///
				 `tstat'=J(`nfitval',2,.)

		   matrix `stat_model'=J(`nfitval',2,.)
		   
		   local counter=1
	       foreach var of local varls {
		      mata: `tstat'[`counter',]=(max(abs((`Xm'[,1]-st_data(.,"`var'")):/`Xm'[,2])), 3)
		      local ++counter
		   }
		   mata: st_matrix("`stat_model'", `tstat')
		   
		   * p values
		   mata: binsreg_pval(`uni_basis', `Xm'[,2], "`tmod_V'", "`stat_model'", `nsims', ///
		                      `nseries', ".", 0, "`pval_model'", ".")
		   mata: mata drop `Xm' `tstat' `uni_basis'
		}
	 }
	 else {
	    local tmod_p=.
		local tmod_s=.
	 }
	 mata: mata drop `uni_grid' `xvec'
	 ****** End of testing *****************************************
	 
	 ******************************
	 ******* Display **************
	 ******************************
	 if ("`knotlist'"!="") {
	     local binselectmethod "User-specified"
		 local placement "User-specified"
	 }
	 else {
	    if ("`binsmethod'"=="DPI") local binselectmethod "IMSE-optimal plug-in choice"
		if ("`binsmethod'"=="ROT") local binselectmethod "IMSE-optimal rule-of-thumb choice"
		if ("`pos'"=="ES") local placement "Evenly-spaced"
		if ("`pos'"=="QS") local placement "Quantile-spaced"
	 }
	 
	 di ""
	 di in smcl in gr "Hypothesis tests based on binscatter estimates"
	 di in smcl in gr "Bin selection method: `binselectmethod'"
	 di in smcl in gr "Placement: `placement'"
	 di in smcl in gr "Derivative: `deriv'"
	 di ""
	 di in smcl in gr "{hline 30}{c TT}{hline 15}"
	 di in smcl in gr "{lalign 1:# of observations}"   _col(30) " {c |} " _col(32) as result %7.0f `N'
	 di in smcl in gr "{lalign 1:# of distinct values}"   _col(30) " {c |} " _col(32) as result %7.0f `Ndist'
	 di in smcl in gr "{lalign 1:# of clusters}"   _col(30) " {c |} " _col(32) as result %7.0f `Nclust'	 
	 di in smcl in gr "{hline 30}{c +}{hline 15}"
	 di in smcl in gr "{lalign 1:Bin selection:}"             _col(30) " {c |} "
	 if ("`binselectmethod'"=="User-specified") {
	    di in smcl in gr "{ralign 29:Degree of polynomial}"         _col(30) " {c |} " _col(39) as result %7.0f "."
	    di in smcl in gr "{ralign 29:# of smoothness constraints}"  _col(30) " {c |} " _col(39) as result %7.0f "."
	 }
	 else {
	    di in smcl in gr "{ralign 29:Degree of polynomial}"         _col(30) " {c |} " _col(32) as result %7.0f `binsp'
	    di in smcl in gr "{ralign 29:# of smoothness constraints}"  _col(30) " {c |} " _col(32) as result %7.0f `binss'
	 }
	 di in smcl in gr "{ralign 29:# of bins}"                       _col(30) " {c |} " _col(32) as result %7.0f `nbins'
	 di in smcl in gr "{hline 30}{c BT}{hline 15}"

	 if ("`tsha_fewobs'"!="T") {
	    if (`ntestshape'!=0) {
		   di ""
	       di in smcl in gr "Shape Restriction Tests:"
	       di in smcl in gr "Degree: `tsha_p'" _col(15) "# of smoothness constraints: `tsha_s'"
	    
	    }
	    if ("`testshapel'"!="") {
		   di ""
	       di in smcl in gr "{hline 19}{c TT}{hline 30}"
	       di in smcl in gr "H0: sup mu <="  _col(20) in gr ///
		                    "{c |}" _col(22) "sup T"  _col(40) "p value"
		   di in smcl in gr "{hline 19}{c +}{hline 30}"
		   forval i=1/`nL' {
	          local val: word `i' of `testshapel'   
		      local stat=`stat_shapeL'[`i',1]
		      local pval=`pval_shapeL'[`i',1]
	          di in smcl in yellow "{rcenter 19:`val'}" _col(20) in gr "{c |}" ///
		                    _col(22) as result %7.3f `stat' ///
		                    _col(40) as result %7.3f `pval'
	       }
		   di in smcl in gr "{hline 19}{c BT}{hline 30}"
	    }
	    
	    if ("`testshaper'"!="") {
		   di ""
	       di in smcl in gr "{hline 19}{c TT}{hline 30}"
		   di in smcl in gr "H0: inf mu >="  _col(20) in gr ///
		                    "{c |}" _col(22) "inf T"  _col(40) "p value"
		   di in smcl in gr "{hline 19}{c +}{hline 30}"	    
		   forval i=1/`nR' {
	          local val: word `i' of `testshaper'
		      local stat=`stat_shapeR'[`i',1]
		      local pval=`pval_shapeR'[`i',1]
	          di in smcl in yellow "{rcenter 19:`val'}" _col(20) in gr "{c |}" ///
		                    _col(22) as result %7.3f `stat' ///
		                    _col(40) as result %7.3f `pval'
	       }
		   di in smcl in gr "{hline 19}{c BT}{hline 30}"
	    }
	    if ("`testshape2'"!="") {
		   di ""
	       di in smcl in gr "{hline 19}{c TT}{hline 30}"
	       di in smcl in gr "H0: mu ="  _col(20) in gr ///
		                    "{c |}" _col(22) "sup |T|"  _col(40) "p value"
		   di in smcl in gr "{hline 19}{c +}{hline 30}"
		   forval i=1/`nT' {
	          local val: word `i' of `testshape2'   
		      local stat=`stat_shape2'[`i',1]
		      local pval=`pval_shape2'[`i',1]
	          di in smcl in yellow "{rcenter 19:`val'}" _col(20) in gr "{c |}" ///
		                    _col(22) as result %7.3f `stat' ///
		                    _col(40) as result %7.3f `pval'
	       } 
		   di in smcl in gr "{hline 19}{c BT}{hline 30}"
	    }	 
	 }
	 
	 if ("`tmod_fewobs'"!="T") {
	    if ("`testmodelpoly'"!=""|`"`testmodelparfit'"'!=`""') {
	       di ""
	       di in smcl in gr "Model specification Tests:"
	       di in smcl in gr "Degree: `tmod_p'" _col(15) "# of smoothness constraints: `tmod_s'"
	    }
        if ("`testmodelpoly'"!="") {
		   di ""
	       local stat=`stat_poly'[1,1]
		   local pval=`pval_poly'[1,1]
		   di in smcl in gr "{hline 19}{c TT}{hline 30}"
		   di in smcl in gr "H0: mu =" _col(20) in gr ///
		                 "{c |}" _col(22) "sup |T|"  _col(40) "p value"
		   di in smcl in gr "{hline 19}{c +}{hline 30}"
	       di in smcl in gr "poly. degree  " as result `testpolyp' _col(20) in gr "{c |}" ///
		                    _col(22) as result %7.3f `stat' ///
		                    _col(40) as result %7.3f `pval'
		   di in smcl in gr "{hline 19}{c BT}{hline 30}"
	    }
	    if (`"`testmodelparfit'"'!=`""') {
	       di ""
		   di in smcl in gr `"Input file: `testmodelparfit'.dta"'
	       di in smcl in gr "{hline 19}{c TT}{hline 30}"
	       di in smcl in gr "H0: mu ="  _col(20) in gr ///
		                    "{c |}" _col(22) "sup |T|"  _col(40) "p value"
		   di in smcl in gr "{hline 19}{c +}{hline 30}"
		   forval i=1/`nfitval' {
	          local val: word `i' of `varls'
		      local stat=`stat_model'[`i',1]
		      local pval=`pval_model'[`i',1]
	          di in smcl in yellow "{rcenter 19:`val'}" _col(20) in gr "{c |}" ///
		                    _col(22) as result %7.3f `stat' ///
		                    _col(40) as result %7.3f `pval'
	       }
		   di in smcl in gr "{hline 19}{c BT}{hline 30}"
	    }
	 }
	 
	 ****************************
	 ******* Return *************
	 ****************************
	 ereturn clear
	 ereturn scalar N=`N'
	 ereturn scalar Ndist=`Ndist'
	 ereturn scalar Nclust=`Nclust'
	 ereturn scalar nbins=`nbins'
	 ereturn scalar p=`binsp'
	 ereturn scalar s=`binss'
	 ereturn scalar testshape_p=`tsha_p'
	 ereturn scalar testshape_s=`tsha_s'
	 ereturn scalar testmodel_p=`tmod_p'
	 ereturn scalar testmodel_s=`tmod_s'
	 
	 if ("`tsha_fewobs'"!="T") {
	    if ("`testshapel'"!="") {
		   ereturn local testvalueL `testshapel'
	       ereturn matrix stat_shapeL=`stat_shapeL'
	       ereturn matrix pval_shapeL=`pval_shapeL'
	    }
	    if ("`testshaper'"!="") {
		   ereturn local testvalueR `testshaper'
	       ereturn matrix stat_shapeR=`stat_shapeR'
	       ereturn matrix pval_shapeR=`pval_shapeR'  
	    }
	    if ("`testshape2'"!="") {
		   ereturn local testvalue2 `testshape2'
	 	   ereturn matrix stat_shape2=`stat_shape2'
	       ereturn matrix pval_shape2=`pval_shape2'
	    }
	 }
	 
	 if ("`tmod_fewobs'"!="T") {
	    if ("`testmodelpoly'"!="") {
	       ereturn scalar testpolyp=`testpolyp'
	       ereturn scalar stat_poly=`stat_poly'[1,1]
		   ereturn scalar pval_poly=`pval_poly'[1,1]
	    }
	    if (`"`testmodelparfit'"'!=`""') {
	       ereturn local testvarlist `varls'
	       ereturn matrix stat_model=`stat_model'
		   ereturn matrix pval_model=`pval_model'
	    }
	 }

end


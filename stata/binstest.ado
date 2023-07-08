*! version 1.3 03-Jul-2023

capture program drop binstest
program define binstest, eclass
    version 13
	 
	syntax varlist(min=2 numeric fv ts) [if] [in] [fw aw pw] [, deriv(integer 0) at(string asis) nolink ///
	       estmethod(string) estmethodopt(string asis) absorb(string asis) reghdfeopt(string asis) ///
		   testmodel(string)  ///
		   testmodelparfit(string asis) testmodelpoly(string) ///
		   testshape(string) ///
		   testshapel(numlist) testshaper(numlist) testshape2(numlist) ///
		   lp(string) ///
		   bins(numlist integer max=2 >=0) nbins(string) ///
		   pselect(numlist integer >=0) sselect(numlist integer >=0) ///
		   binspos(string) binsmethod(string) nbinsrot(string) randcut(numlist max=1 >=0 <=1) ///
		   nsims(integer 500) simsgrid(integer 20) simsseed(numlist integer max=1 >=0) ///
		   dfcheck(numlist integer max=2 >=0) masspoints(string) usegtools(string) ///
		   vce(passthru) asyvar(string) ///
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
	 * which model?	 
	 if ("`absorb'"!="") {
	    if ("`estmethod'"!="") {
		   if ("`estmethod'"!="reghdfe") {
		      di as error "absorb() can only be combined with estmethod(reghdfe)."
			  exit
		   }
		}
	    else local estmethod "reghdfe"
	 }
	 if ("`estmethod'"=="") local estmethod "reg"
	 tokenize `estmethod'
	 local estmethod `1'
	 if ("`estmethod'"=="reg") {
	    local estcmd "reg"
	 } 
	 else if ("`estmethod'"=="qreg") {
	    local estcmd "qreg"
		local quantile `2'
		if ("`quantile'"=="") local quantile=0.5
	 }
	 else if ("`estmethod'"=="logit") {
	    local estcmd "logit"
	 }
	 else if ("`estmethod'"=="probit") {
	    local estcmd "probit"
	 }
	 else if ("`estmethod'"=="reghdfe") {
	    local estcmd "reghdfe"
	 }
	 	 	 
	 * report the results for the cond. mean model?
	 if ("`link'"!="") local transform "F"
	 else              local transform "T"

	 * default vce
	 if ("`vce'"=="") local vce "vce(robust)"
	 * vce for bin selection
	 if ("`estmethod'"=="qreg") {
	    if ("`vce'"=="vce(iid)") local vce_select "vce(ols)"
		else                     local vce_select "vce(robust)"
	 }
	 else if ("`estmethod'"=="logit"|"`estmethod'"=="probit") {
	    if ("`vce'"=="oim"|"`vce'"=="opg") local vce_select "vce(ols)"
		else                               local vce_select "`vce'"
	 }
	 else if ("`estmethod'"=="reg"|"`estmethod'"=="reghdfe") {
	    local vce_select "`vce'"
	 }
	 
	 * use bootstrap cmd? cluster specified?
	 local vcetemp: subinstr local vce "vce(" "", all
     local vcetemp: subinstr local vcetemp ")" "", all
	 tokenize "`vcetemp'", parse(", ")
	 if ("`1'"=="boot" | "`1'"=="bootstrap") {
		local boot "on"
		local repstemp `3'
		if ("`repstemp'"=="") local repstemp reps(20)
		local repstemp: subinstr local repstemp "reps(" "", all
        local reps: subinstr local repstemp ")" "", all
		if ("`estmethod'"=="qreg") {
		   local estcmd "bsqreg"
		   if ("`weight'"!="") {
		      di as error "Weights not allowed for bootstrapping."
		      exit
		   }
		}
	 }
	 else if ("`1'"=="cl"|"`1'"=="clu"|"`1'"=="clus"|"`1'"=="clust"| /// 
		 "`1'"=="cluste"|"`1'"=="cluster") {
		if ("`3'"==""|"`3'"==",") local clusterON "T"           /* cluster is specified */
		local clustervar `2'
		local boot "off"
	 } 
	 else {
		local boot "off"
	 }  

	 if ("`asyvar'"=="") local asyvar "off"
	 
	 if ("`binspos'"=="es") local binspos "ES"
	 if ("`binspos'"=="qs") local binspos "QS"
	 if ("`binspos'"=="")   local binspos "QS"
	 if ("`binsmethod'"=="rot") local binsmethod "ROT"
	 if ("`binsmethod'"=="dpi") local binsmethod "DPI"
	 if ("`binsmethod'"=="")    local binsmethod "DPI"

	 
	 * analyze options related to J, p and s
	 if ("`testshape'"!="T"&"`testshape'"!="F"&"`testshape'"!="") {
	    numlist "`testshape'", integer max(2) range(>=0)
		local testshape=r(numlist)
	 }
	 if ("`testmodel'"!="T"&"`testmodel'"!="F"&"`testmodel'"!="") {
	    numlist "`testmodel'", integer max(2) range(>=0)
		local testmodel=r(numlist)
	 }

	 if ("`testshape'"=="F")   local testshape ""
	 if ("`testmodel'"=="F")   local testmodel ""
	 
	 local selection ""
	 
	 * analyze nbins
	 if ("`nbins'"=="T") local nbins=0
	 local len_nbins=0
	 if ("`nbins'"!=""&"`nbins'"!="F") {
	    numlist "`nbins'", integer sort
	    local nbins=r(numlist)
	    local len_nbins: word count `nbins'
	 }
	 * shut down selection if knot is specified by users
	 if ("`binspos'"!="ES"&"`binspos'"!="QS") {
		if ("`nbins'"!=""|"`pselect'"!=""|"`sselect'"!="") {
		   di as error "nbins(), pselect() or sselect() incorrectly specified."
		   exit
		}
	 }

	 
	 * analyze numlist in pselect and sselect
	 local len_p=0
	 local len_s=0
	 
	 if ("`pselect'"!="") {
	    numlist "`pselect'", integer range(>=`deriv') sort
		local plist=r(numlist)
	 }
	 
	 if ("`sselect'"!="") {
	    numlist "`sselect'", integer range(>=0) sort
		local slist=r(numlist)
	 }
	 	 	 
	 local len_p: word count `plist'
	 local len_s: word count `slist'
	 
	 if (`len_p'==1&`len_s'==0) {
	    local slist `plist'
		local len_s=1
	 }
	 if (`len_p'==0&`len_s'==1) {
	    local plist `slist'
		local len_p=1
     }
	 
	 tokenize `bins'
	 local binsp "`1'"
	 local binss "`2'"
	 if ("`binsp'"=="") local binsp=.
	 if ("`binss'"=="") local binss `binsp'
	 if ("`bins'"!="") {
	 	if ("`nbins'"!=""&"`nbins'"!="T"&"`nbins'"!="0"&`len_nbins'<=1) {
			di as error "nbins() or bins() is incorrectly specified."
			exit
		}
	 }

	 * 1st case: select J
	 if (("`bins'"!=""|"`nbins'"=="0"|`len_nbins'>1|"`nbins'"=="T"|"`nbins'"=="")&("`binspos'"=="ES"|"`binspos'"=="QS")) {
	    local selection "J"
	 }
	 
	 if ("`selection'"=="J") {
	 	if (`len_p'>1|`len_s'>1) {
		   if ("`nbins'"=="") {
		      di as error "nbins() must be specified for degree/smoothness selection."
			  exit
		   }
		   else {
			  di as error "Only one p and one s are allowed to select # of bins."
			  exit
		   }
		}
	 	if ("`plist'"=="") local plist=`deriv'
		if ("`slist'"=="") local slist=`plist'
		if ("`bins'"=="") {
		   local binsp `plist'
		   local binss `slist'
		}
		local len_p=1
	    local len_s=1
	    if ("`testshape'"=="T"|"`testshape'"=="") local testshape `=`binsp'+1' `=`binss'+1'
		if ("`testmodel'"=="T"|"`testmodel'"=="") local testmodel `=`binsp'+1' `=`binss'+1'
     }                                                                          
	 
	 * 2nd case: select P (the special case with nbins() pselect() will be modified in the next step) 
	 if ("`selection'"!="J" & ("`testshape'"==""|"`testshape'"=="T"|"`testmodel'"==""|"`testmodel'"=="T")) {
	    local pselectOK "T"     
	 }

	 if ("`pselectOK'"=="T" & `len_nbins'==1 & (`len_p'>1|`len_s'>1)) {
	    local selection "P"
	    *if ("`plist'"=="") {
		 *  numlist "`=max(`deriv', 0)'/4"
		  * local plist=r(numlist)
		*}
	 }                                                                          
	 
	 * 3rd case: user-specified J and p
	 *if ("`testshape'"!="T"&"`testmodel'"!="T") local userOK "T"
	 if ((`len_p'<=1&`len_s'<=1) & "`selection'"!="J") {
		local selection "NA" 
	    if ("`testshape'"=="") {
		   if ("`bins'"!="")   local testshape `=`binsp'+1' `=`binss'+1'                 
		   else {
		      if (`len_p'==1&`len_s'==1) local testshape `=`plist'+1' `=`slist'+1'
		      else                       local testshape `=`deriv'+1' `=`deriv'+1'
		   }
		}
		if ("`testmodel'"=="") {
           if ("`bins'"!="")  local testmodel `=`binsp'+1' `=`binss'+1' 
		   else {
		      if (`len_p'==1&`len_s'==1) local testmodel `=`plist'+1' `=`slist'+1'   
		      else                       local testmodel `=`deriv'+1' `=`deriv'+1'
		   }
		}
	 }
	 	 
	 * exclude all other cases
	 if ("`selection'"=="") {
	    di as error "Degree, smoothness, or # of bins are not correctly specified."
        exit
	 }
	 
	 * Option for testing shape
	 tokenize `testshape'
	 local tsha_p "`1'"
	 local tsha_s "`2'"
	 if ("`tsha_p'"==""|"`tsha_p'"=="T") local tsha_p=.
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
	 if ("`tmod_p'"==""|"`tmod_p'"=="T") local tmod_p=.
	 if ("`tmod_s'"=="") local tmod_s `tmod_p'
	 
	 
	 
	 * Add warnings about degrees for estimation and inference
	 if ("`selection'"=="J") {
	    if ("`tsha_p'"!=".") {
		   if (`tsha_p'<=`binsp') {
		      local tsha_p=`binsp'+1
			  local tsha_s=`tsha_p'
			  di as text "Warning: Degree for testshape() has been changed. It must be greater than the degree for bin selection."
		   }
		}
	    if ("`tmod_p'"!=".") {
		   if (`tmod_p'<=`binsp') {
		      local tmod_p=`binsp'+1
			  local tmod_s=`tmod_p'
			  di as text "Warning: Degree for testmodel() has been changed. It must be greater than the degree for bin selection."
		   }
		}
	 }
	 if ("`selection'"=="NA") {
		di as text "Warning: Testing procedures are valid when nbins() is much larger than the IMSE-optimal choice. Compare your choice with the IMSE-optimal one obtained by binsregselect."
	 }

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
	 
	 * evaluate at w from another dataset?
	 if (`"`at'"'!=`""'&`"`at'"'!=`"mean"'&`"`at'"'!=`"median"'&`"`at'"'!=`"0"') local atwout "user"
	 
	 * default for lp metric
	 if ("`lp'"=="") local lp "inf"
	 
	 * use gtools commands instead?
	 if ("`usegtools'"=="off") local usegtools ""
	 if ("`usegtools'"=="on")  local usegtools usegtools
	 if ("`usegtools'"!="") {
	    capture which gtools
		if (_rc) {
		   di as error "Gtools package not installed."
		   exit
		}
		local localcheck "F"
		local sel_gtools "on"
	 }
	 else local sel_gtools "off"
	 
	 * use reghdfe?
	 if ("`absorb'"!="") {
	    capture which reghdfe
		if (_rc) {
		   di as error "reghdfe not installed."
		   exit
		}
	 }
	 
	 * Error check
	 if (`"`testmodelparfit'"'==`""'&`ntestshape'==0&"`testmodelpoly'"=="") {
	     di as error "No tests specified."
		 exit
	 }
	 if (`tsha_p'<`tsha_s'|`tmod_p'<`tmod_s'|`binsp'<`binss') {
	     di as error "p cannot be smaller than s."
		 exit
	 }
	 if ("`tsha_p'"!="."&"`binsp'"!=".") {
	    if (`tsha_p'<=`binsp') {
	       di as text in gr "Warning: p for testing <= p for bins() not suggested."
	    }
	 }
	 if ("`tmod_p'"!="."&"`binsp'"!=".") {
	    if (`tmod_p'<=`binsp') {
	        di as text in gr "Warning: p for testing <= p for bins() not suggested."
	    }
	 }
	 if (`tsha_p'<`deriv'|`tmod_p'<`deriv') {
	    di as error "p for test cannot be smaller than deriv."
		exit
	 }
	 if ("`testmodelpoly'"!="") {
		if (`testpolyp'<`deriv') {
		   di as error "Degree of polynomial model cannot be smaller than deriv."
		   exit
		}
	 }
	 if (`nsims'<2000|`simsgrid'<50) {
	    di as text "Note: Setting at least nsims(2000) and simsgrid(50) is recommended to obtain the final results."
	 }
	 
	 * Mark sample
	 preserve
	 
	 * Parse varlist into y_var, x_var and w_var
	 tokenize `varlist'
	 fvrevar `1', tsonly
	 local y_var "`r(varlist)'"
	 fvrevar `2', tsonly
	 local x_var "`r(varlist)'"
	 
	 macro shift 2
	 local w_var "`*'"
	 
	 * read eval point for w from another file
	 if ("`atwout'"=="user") {
	    append using `at'
	 }
	 
	 fvrevar `w_var', tsonly
	 local w_var "`r(varlist)'"
	 local nwvar: word count `w_var'
	 
	 * Save the last obs in a vector and then drop it
	 tempname wuser                  /* a vector used to keep eval for w */
	 if ("`atwout'"=="user") {
	    mata: st_matrix("`wuser'", st_data(`=_N', "`w_var'")) 
	    qui drop in `=_N'
	 }
	 
	 * Get positions of factor vars
     local indexlist ""
     local i = 1
     foreach v in `w_var' {
        if strpos("`v'", ".") == 0 {
           local indexlist  `indexlist' `i'
        }
        local ++i
     }
	 
	 * add a default for at
	 if (`"`at'"'==""&`nwvar'>0) {
	    local at "mean"
	 }
	 
	 marksample touse      /* now renew the mark to account for missing values */
	 qui keep if `touse'
     local eN=_N
	 local nsize=_N     /* # of rows in the original dataset */
	     
	 if ("`usegtools'"==""&("`masspoints'"!="off"|"`binspos'"=="QS")) {
	    if ("`:sortedby'"!="`x_var'") sort `x_var', stable
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
		   if ("`usegtools'"=="") {
	          mata: `binedges'=binsreg_uniq(`xvec', ., 1, "Ndist")
		      mata: mata drop `binedges'
		   }
		   else {
		      qui gunique `x_var'
			  local Ndist=r(unique)
		   }
		}   
		local eN=min(`eN', `Ndist')
	 }
	 
	 local Nclust=.
	 if ("`clusterON'"=="T") {
		if ("`numclust'"!=""&"`numclust'"!=".") {
		   local Nclust=`numclust'
	    }
		else {
		   if ("`usegtools'"=="") {
		      mata: st_local("Nclust", strofreal(rows(uniqrows(st_data(.,"`clustervar'")))))
		   }
		   else {
		      qui gunique `clustervar'
			  local Nclust=r(unique)
		   }
	    }
	    if ("`estmethod'"=="qreg") {
			local vce "vce(robust)"
			di as text in gr "Warning: vce(cluster) not allowed. vce(robust) used instead."
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
			 di as error "Inner knots specified out of allowed range."
			 exit
		  }
		  else {
		     local nbins=`nbins'+1
			 local binspos "user"
		  }
	   }
	   else {
		  di as error "Numeric list incorrectly specified in binspos()."
		  exit
	   }
	 }
	 
	 * if binsmethod is specified
	 local imse_bsq_rot=.
	 local imse_var_rot=.
	 local imse_bsq_dpi=.
	 local imse_var_dpi=.
	 if ("`selection'"!="NA") {
		 * Check effective sample size
		 if ("`binsp'"==".") local binspcheck=6
		 else                local binspcheck=`binsp'
		 if ("`nbinsrot'"==""&(`eN'<=`dfcheck_n1'+`binspcheck'+1+`qrot')) {
		    * ROT unavailable, exit
			di as error "Too few observations for bin selection."
	        exit
	     }
		 else {
		   local randcut1k `randcut'
		   if ("`randcut'"=="" & `N'>5000) {
	           local randcut1k=max(5000/`N', 0.01)
		       di as text in gr "Warning: To speed up computation, bin/degree selection uses a subsample of roughly max(5000, 0.01n) observations if n>5000. To use the full sample, set randcut(1)."
	       }

		   if ("`selection'"=="J") {
		      qui binsregselect `y_var' `x_var' `w_var' `wt', deriv(`deriv') bins(`binsp' `binss') nbins(`nbins') ///
				   			    absorb(`absorb') reghdfeopt(`reghdfeopt') ///
								binsmethod(`binsmethod') binspos(`binspos') nbinsrot(`nbinsrot') ///
								`vce_select' masspoints(`masspoints') dfcheck(`dfcheck_n1' `dfcheck_n2') ///
								numdist(`Ndist') numclust(`Nclust') randcut(`randcut1k') usegtools(`sel_gtools')
			  if (e(nbinsrot_regul)==.) {
			      di as error "bin selection fails."
				  exit
			  }
			  if ("`binsmethod'"=="ROT") {
				  local nbins=e(nbinsrot_regul)
				  local imse_bsq_rot=e(imse_bsq_rot)
				  local imse_var_rot=e(imse_var_rot)
			  }
			  else if ("`binsmethod'"=="DPI") {
				 local nbins=e(nbinsdpi)
				 local imse_bsq_dpi=e(imse_bsq_dpi)
				 local imse_var_dpi=e(imse_var_dpi)
				 if (`nbins'==.) {
				    local nbins=e(nbinsrot_regul)
					local imse_bsq_rot=e(imse_bsq_rot)
				    local imse_var_rot=e(imse_var_rot)
					di as text in gr "Warning: DPI selection fails. ROT choice used."
				 }
			  }
		   }
		   else if ("`selection'"=="P") {
		   	  qui binsregselect `y_var' `x_var' `w_var' `wt', deriv(`deriv') nbins(`nbins') ///
				   			    absorb(`absorb') reghdfeopt(`reghdfeopt') ///
								pselect(`plist') sselect(`slist') ///
								binsmethod(`binsmethod') binspos(`binspos') nbinsrot(`nbinsrot') ///
								`vce_select' masspoints(`masspoints') dfcheck(`dfcheck_n1' `dfcheck_n2') ///
								numdist(`Ndist') numclust(`Nclust') randcut(`randcut1k') usegtools(`sel_gtools')
			  if (e(prot_regul)==.) {
			      di as error "Bin selection fails."
				  exit
			  }
			  if ("`binsmethod'"=="ROT") {
				  local binsp=e(prot_regul)
				  local binss=e(srot_regul)
				  local imse_bsq_rot=e(imse_bsq_rot)
				  local imse_var_rot=e(imse_var_rot)
			  }
			  else if ("`binsmethod'"=="DPI") {
				 local binsp=e(pdpi)
				 local binss=e(sdpi)
				 local imse_bsq_dpi=e(imse_bsq_dpi)
				 local imse_var_dpi=e(imse_var_dpi)
				 if (`binsp'==.) {
				    local binsp=e(prot_regul)
					local binss=e(srot_regul)
					local imse_bsq_rot=e(imse_bsq_rot)
				    local imse_var_rot=e(imse_var_rot)
					di as text in gr "Warning: DPI selection fails. ROT choice used."
				 }
			  }
			  if ("`testshape'"=="T"|"`testshape'"=="") {
			     local tsha_p=`binsp'+1
				 local tsha_s=`binss'+1
			  }
			  else {
			     if (`tsha_p'<=`binsp') {
				   local tsha_p=`binsp'+1
				   local tsha_s=`tsha_p'
				   di as text "Warning: Degree for testshape() has been changed. It must be greater than the IMSE-optimal degree."
				 }
			  }
			  if ("`testmodel'"=="T"|"`testmodel'"=="") {
			     local tmod_p=`binsp'+1
				 local tmod_s=`binss'+1
			  }
			  else {
			     if (`tmod_p'<=`binsp') {
			        local tmod_p=`binsp'+1
				    local tmod_s=`tmod_p'
				    di as text "Warning: Degree for testmodel() has been changed. It must be greater than the IMSE-optimal degree."
			     }
			  }
		   }
		}
	 }
	 
	 *******************************************************
	 * Check if eff. sample size is large enough for testing
	 if ((`nbins'-1)*(`tsha_p'-`tsha_s'+1)+`tsha_p'+1+`dfcheck_n2'>=`eN') {
	    local tsha_fewobs "T"
		di as text in gr "Warning: Too small effective sample size for testing shape."
	 }
	 if ((`nbins'-1)*(`tmod_p'-`tmod_s'+1)+`tmod_p'+1+`dfcheck_n2'>=`eN') {
	    local tmod_fewobs "T"
		di as text "Warning: Too small effective sample size for testing models."
	 }
	 ********************************************************
	 	 	 
	 * Generate category variable for data and save knot in matrix
	 tempname kmat
	 tempvar xcat
	 qui gen `xcat'=. in 1
		
     if ("`binspos'"=="ES") {
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
	       binsreg_pctile `x_var' `wt', nq(`nbins') `usegtools'
		   mat `kmat'=(`xmin' \ r(Q) \ `xmax')
		}
	 }
	 
	 mata: st_matrix("`kmat'", (`xmin' \ uniqrows(st_matrix("`kmat'")[|2 \ `=`nbins'+1'|])))
	 
	 binsreg_irecode `x_var', knotmat(`kmat') bin(`xcat') ///
		                      `usegtools' nbins(`nbins') pos(`binspos') knotliston(T)
	 
	 if (`nbins'!=rowsof(`kmat')-1) {
	     di as text in gr "Warning: Repeated knots. Some bins dropped."
		 local nbins=rowsof(`kmat')-1
	 }
	 
	 * Check for empty bins
	 if ("`localcheck'"=="T") {
	   mata: st_local("Ncat", strofreal(rows(uniqrows(st_data(.,"`xcat'")))))
	   if (`nbins'==`Ncat') {
		  mata: `binedges'=binsreg_uniq(`xvec', st_data(.,"`xcat'"), `nbins', "uniqmin")
		  mata: mata drop `binedges'
	   }
	   else {
		  local uniqmin=0
		  di as text in gr "Warning: There are empty bins. Specify a smaller number in nbins()."
	   }
	   
	   if (`ntestshape'!=0) {
		  if (`uniqmin'<`tsha_p'+1) {
		     local tsha_fewobs "T"
		     di as text in gr "Warning: Some bins have too few distinct x-values for testing."
		  }
	   }
	   if (`"`testmodelparfit'"'!=`""'|"`testmodelpoly'"!="") {
		  if (`uniqmin'<`tmod_p'+1) {
		     local tmod_fewobs "T"
		     di as text in gr "Warning: Some bins have too few distinct x-values for testing."
		  }
	   }
	 }
	 
	 ********************************************************
	 * Set seed
	 if ("`simsseed'"!="") set seed `simsseed'
	 local uni_last=`simsgrid'*`nbins'+`nbins'-1
	 
	 tempname Xm Xm0 fit se fit0 uni_grid uni_basis tstat vcov           /* objects in MATA */
	 mata: `uni_grid'=binsreg_grids("`kmat'", `simsgrid')
	 mata: `Xm0'=.; `fit'=.; `fit0'=0; `se'=.; `vcov'=.
	 
	 * adjust w vars
	 tempname wval
	 if (`nwvar'>0) {
	   if (`"`at'"'==`"mean"'|`"`at'"'==`"median"') {
	      matrix `wval'=J(1, `nwvar', 0)
	 	  tempname wvaltemp mataobj
		  mata: `mataobj'=.
		  foreach wpos in `indexlist' {
			 local wname: word `wpos' of `w_var'
			 if ("`usegtools'"=="") {
		        if ("`wtype'"!="") qui tabstat `wname' `conds' [aw`exp'], stat(`at') save
			    else               qui tabstat `wname' `conds', stat(`at') save
			    mat `wvaltemp'=r(StatTotal)
			 }
			 else {
			    qui gstats tabstat `wname' `conds' `wt', stat(`at') matasave("`mataobj'")
				mata: st_matrix("`wvaltemp'", `mataobj'.getOutputCol(1))
			 }
			 mat `wval'[1,`wpos']=`wvaltemp'[1,1]
		  }
		  mata: mata drop `mataobj'
	   }
	   else if (`"`at'"'==`"0"') {
   		  matrix `wval'=J(1,`nwvar',0)
	   }
	   else if ("`atwout'"=="user") {
		  matrix `wval'=`wuser'
	   }
	}
	
	* define a w vector (possibly a constant) in MATA
	tempname wvec wvec0
	mata: `wvec'=J(1,0,.); `wvec0'=J(1,0,.)
	if (`nwvar'>0) {
	   mata: `wvec0'=st_matrix("`wval'")
	   if (`deriv'==0&"`asyvar'"=="off") mata: `wvec'=(`wvec', `wvec0')
	   else                              mata: `wvec'=(`wvec', J(1,`nwvar',0))
	}
	if ("`estmethod'"=="qreg"|"`estmethod'"=="reghdfe") {
       mata: `wvec0'=(`wvec0', 1)
	   if (`deriv'==0) mata: `wvec'=(`wvec', 1)
	   else            mata: `wvec'=(`wvec', 0)
	}
	
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
	    if ("`estmethod'"!="qreg"&"`estmethod'"!="reghdfe") {
		    capture `estcmd' `y_var' `tsha_series' `w_var' `wt', nocon `vce' `estmethodopt'
		}
		else if ("`estmethod'"=="qreg") {
		   if ("`boot'"=="on") capture bsqreg `y_var' `tsha_series' `w_var', quantile(`quantile') reps(`reps')
		   else                capture qreg `y_var' `tsha_series' `w_var' `wt', quantile(`quantile') `vce' `estmethodopt'
		}
		else {
		   capture `estcmd' `y_var' `tsha_series' `w_var' `wt', absorb(`absorb') `reghdfeopt' `vce' 
		}
		
	    * store results
	    if (_rc==0) {
		    matrix `tsha_b'=e(b)
		    matrix `tsha_V'=e(V)
			if ("`estmethod'"!="qreg"&"`estmethod'"!="reghdfe") mata: binsreg_checkdrop("`tsha_b'", "`tsha_V'", `nseries')
			else                                                mata: binsreg_checkdrop("`tsha_b'", "`tsha_V'", `nseries', "T")
			matrix `tsha_b'=`tsha_b''
	    }
	    else {
	        error  _rc
	   	    exit _rc
        }
		   
		* Predict		
	    * fitted values & standard errors
		mata: `uni_basis'=binsreg_spdes(`uni_grid'[,1], "`kmat'", `uni_grid'[,3], `tsha_p', `deriv', `tsha_s')
		if (("`estmethod'"=="logit"|"`estmethod'"=="probit")&"`transform'"=="T") {
		   if (`deriv'==0) {
		      mata: `fit0'=(`uni_basis', J(rows(`uni_basis'),1,1)#`wvec0')*st_matrix("`tsha_b'")
			  if ("`estmethod'"=="logit") {
			     mata: `fit'=logistic(`fit0'); ///
				       `se'=logisticden(`fit0'):* ///
					        binsreg_pred((`uni_basis', J(rows(`uni_basis'),1,1)#`wvec'),.,st_matrix("`tsha_V'"),"se")[,2]
			  }
			  else {
			     mata: `fit'=normal(`fit0'); ///
				       `se'=normalden(`fit0'):* ///
					        binsreg_pred((`uni_basis', J(rows(`uni_basis'),1,1)#`wvec'),.,st_matrix("`tsha_V'"),"se")[,2]
			  }
		   }
		   if (`deriv'==1) {
		      mata: `Xm0'=binsreg_spdes(`uni_grid'[,1], "`kmat'", `uni_grid'[,3], `tsha_p', 0, `tsha_s'); ///
			        `Xm0'=(`Xm0', J(rows(`Xm0'),1,1)#`wvec0'); ///
					`fit0'=`Xm0'*st_matrix("`tsha_b'"); ///
					`Xm'=(`uni_basis', J(rows(`uni_basis'),1,1)#`wvec')
			  if ("`estmethod'"=="logit") {
			     mata: `fit'=binsreg_pred(`Xm',st_matrix("`tsha_b'"),.,"xb")[,1]
				 if ("`asyvar'"=="off") {
					mata: `Xm'=logisticden(`fit0'):*(1:-2*logistic(`fit0')):*`fit':*`Xm0' + ///
					           logisticden(`fit0'):*`Xm'; /// 
				          `se'=sqrt(rowsum((`Xm'*st_matrix("`tsha_V'")):*`Xm'))
				 }
				 else {
				    mata: `se'=logisticden(`fit0'):*(binsreg_pred(`Xm',.,st_matrix("`tsha_V'"),"se")[,2])
				 }
				 mata: `fit'=logisticden(`fit0'):*`fit'
			  }
			  else {
			     mata: `fit'=binsreg_pred(`Xm',st_matrix("`tsha_b'"),.,"xb")[,1]
				 if ("`asyvar'"=="off") {
					 mata:`Xm'=(-`fit0'):*normalden(`fit0'):*`fit':*`Xm0' + ///
                                normalden(`fit0'):*`Xm'; ///
						  `se'=sqrt(rowsum((`Xm'*st_matrix("`tsha_V'")):*`Xm'))
				 }
				 else {
				    mata: `se'=normalden(`fit0'):*(binsreg_pred(`Xm',.,st_matrix("`tsha_V'"),"se")[,2])
				 }
				 mata: `fit'=normalden(`fit0'):*`fit'
			  }
		   }
		   mata: `Xm'=(`fit', `se')
		}
		else {
		   mata: `Xm'=(`uni_basis', J(rows(`uni_basis'),1,1)#`wvec'); ///
		         `Xm'=binsreg_pred(`Xm', st_matrix("`tsha_b'"), st_matrix("`tsha_V'"), "all")
		}
		
		* Test statistics
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
			  if ("`lp'"=="inf") {
			     mata: `tstat'[`i',.]=(max(abs((`Xm'[,1]:-`val'):/`Xm'[,2])), 3)
			  }
			  else {
			     mata: `tstat'[`i',.]=(mean(abs((`Xm'[,1]:-`val'):/`Xm'[,2]):^`lp')^(1/`lp'), 3)
			  }
		   }
		}
		mata: st_matrix("`stat_shape'", `tstat')
		
		* p value
		if ("`estmethod'"=="qreg"|"`estmethod'"=="reghdfe") {
		   if (`deriv'==0) mata: `uni_basis'=(`uni_basis', J(rows(`uni_basis'),1,1))
		   else            mata: `uni_basis'=(`uni_basis', J(rows(`uni_basis'),1,0))
		   mata: `vcov'=st_matrix("`tsha_V'"); ///
		         `vcov'= (`vcov'[|1,1 \ `nseries', `nseries'|], `vcov'[|1,cols(`vcov') \ `nseries', cols(`vcov')|] \ ///
				          `vcov'[|cols(`vcov'), 1 \ cols(`vcov'), `nseries'|], `vcov'[cols(`vcov'), cols(`vcov')]); ///
			     st_matrix("`vcov'", `vcov')
		}
		
		if ("`estmethod'"!="qreg"&"`estmethod'"!="reghdfe") {
		   mata: `Xm'=binsreg_pred(`uni_basis', ., st_matrix("`tsha_V'")[|1,1 \ `nseries',`nseries'|], "se"); ///
		         binsreg_pval(`uni_basis', `Xm'[,2], "`tsha_V'", "`stat_shape'", `nsims', `nseries', ///
		                      ".", 0, "`pval_shape'", ".", "`lp'")
		}
		else {
		   mata: `Xm'=binsreg_pred(`uni_basis', ., `vcov', "se"); ///
		         binsreg_pval(`uni_basis', `Xm'[,2], "`vcov'", "`stat_shape'", `nsims', `=`nseries'+1', ///
		                      ".", 0, "`pval_shape'", ".", "`lp'")
		}
		
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
		
		tempname tmod_b tmod_V
		capture confirm matrix `tsha_b' `tsha_V'
		if (_rc==0&`tmod_p'==`tsha_p'& `tmod_s'==`tsha_s') {
		    matrix `tmod_b'=`tsha_b'
		    matrix `tmod_V'=`tsha_V'
		} 
		else {
		   local tmod_series ""
	       forvalues i=1/`nseries' {
	          tempvar sp`i'
	          local tmod_series `tmod_series' `sp`i''
		      qui gen `sp`i''=. in 1
	       }
		
		   mata: binsreg_st_spdes(`xvec', "`tmod_series'", "`kmat'", st_data(.,"`xcat'"), `tmod_p', 0, `tmod_s')
		   if ("`estmethod'"!="qreg"&"`estmethod'"!="reghdfe") {
	          capture `estcmd' `y_var' `tmod_series' `w_var' `wt', nocon `vce' `estmethodopt'
		   }
		   else if ("`estmethod'"=="qreg") {
		      if ("`boot'"=="on") capture bsqreg `y_var' `tmod_series' `w_var', quantile(`quantile') reps(`reps')
			  else                capture qreg `y_var' `tmod_series' `w_var' `wt', quantile(`quantile') `vce' `estmethodopt'
		   }
		   else {
		      capture `estcmd' `y_var' `tmod_series' `w_var' `wt', absorb(`absorb') `reghdfeopt' `vce' 
		   }
	 
	       * store results
	       if (_rc==0) {
		       matrix `tmod_b'=e(b)
		       matrix `tmod_V'=e(V)
			   if ("`estmethod'"!="qreg"&"`estmethod'"!="reghdfe") mata: binsreg_checkdrop("`tmod_b'", "`tmod_V'", `nseries')
			   else                                                mata: binsreg_checkdrop("`tmod_b'", "`tmod_V'", `nseries', "T")
	           matrix `tmod_b'=`tmod_b''
		   }
	       else {
	           error  _rc
	   	       exit _rc
           }
		   
		   drop `tmod_series'
		}
		
		
		********************************************************
		* If a test for poly reg is requested
		if ("`testmodelpoly'"!="") {
		   	* fitted values
  	       mata: `uni_basis'=binsreg_spdes(`uni_grid'[,1], "`kmat'", `uni_grid'[,3], `tmod_p', `deriv', `tmod_s')
		   
		   if (("`estmethod'"=="logit"|"`estmethod'"=="probit")&"`transform'"=="T") {
		      if (`deriv'==0) {
		         mata: `fit0'=(`uni_basis', J(rows(`uni_basis'),1,1)#`wvec0')*st_matrix("`tmod_b'")
			     if ("`estmethod'"=="logit") {
			        mata: `fit'=logistic(`fit0'); ///
				          `se'=logisticden(`fit0'):* ///
					           binsreg_pred((`uni_basis', J(rows(`uni_basis'),1,1)#`wvec'),.,st_matrix("`tmod_V'"),"se")[,2]
			     }
			     else {
			        mata: `fit'=normal(`fit0'); ///
				          `se'=normalden(`fit0'):* ///
					           binsreg_pred((`uni_basis', J(rows(`uni_basis'),1,1)#`wvec'),.,st_matrix("`tmod_V'"),"se")[,2]
			     }
		      }
		      if (`deriv'==1) {
		         mata: `Xm0'=binsreg_spdes(`uni_grid'[,1], "`kmat'", `uni_grid'[,3], `tmod_p', 0, `tmod_s'); ///
			           `Xm0'=(`Xm0', J(rows(`Xm0'),1,1)#`wvec0'); ///
					   `fit0'=`Xm0'*st_matrix("`tmod_b'"); ///
					   `Xm'=(`uni_basis', J(rows(`uni_basis'),1,1)#`wvec')
			     if ("`estmethod'"=="logit") {
			        mata: `fit'=binsreg_pred(`Xm',st_matrix("`tmod_b'"),.,"xb")[,1]
				    if ("`asyvar'"=="off") {
				       mata: `Xm'=logisticden(`fit0'):*(1:-2*logistic(`fit0')):*`fit':*`Xm0' + ///
				                  logisticden(`fit0'):*`Xm'; ///
				             `se'=sqrt(rowsum((`Xm'*st_matrix("`tmod_V'")):*`Xm'))
				    }
				    else {
				       mata: `se'=logisticden(`fit0'):*(binsreg_pred(`Xm',.,st_matrix("`tmod_V'"),"se")[,2])
				    }
				    mata: `fit'=logisticden(`fit0'):*`fit'
			     }
			     else {
			        mata: `fit'=binsreg_pred(`Xm',st_matrix("`tmod_b'"),.,"xb")[,1]
				    if ("`asyvar'"=="off") {
					    mata:`Xm'=(-`fit0'):*normalden(`fit0'):*`fit':*`Xm0' + ///
                                  normalden(`fit0'):*`Xm'; ///
						     `se'=sqrt(rowsum((`Xm'*st_matrix("`tmod_V'")):*`Xm'))
				    }
				    else {
				       mata: `se'=normalden(`fit0'):*(binsreg_pred(`Xm',.,st_matrix("`tmod_V'"),"se")[,2])
				    }
				    mata: `fit'=normalden(`fit0'):*`fit'
			     }
		      }
		      mata: `Xm'=(`fit', `se')
		   }
		   else {
		      mata: `Xm'=(`uni_basis', J(rows(`uni_basis'),1,1)#`wvec'); ///
		            `Xm'=binsreg_pred(`Xm', st_matrix("`tmod_b'"), st_matrix("`tmod_V'"), "all")
		   }
		   		   
		   * Polynomial fit
		   tempvar poly_fit
	       local poly_series ""
		   *if ("`estmethod'"!="qreg"&"`estmethod'"!="reghdfe") local ini=0
		   *else                                                local ini=1
	       forval i=1/`testpolyp' {
		      tempvar x_var_`i'
			  qui gen `x_var_`i''=`x_var'^`i'
	          local poly_series `poly_series' `x_var_`i''
		   }
		 
		   if ("`estmethod'"!="qreg"&"`estmethod'"!="reghdfe") {
		      capture `estcmd' `y_var' `poly_series' `w_var' `wt', `estmethodopt'
		   }
		   else if ("`estmethod'"=="qreg") {
		      capture qreg `y_var' `poly_series' `w_var' `wt', quantile(`quantile') `estmethodopt'
		   }
		   else {
		      capture `estcmd' `y_var' `poly_series' `w_var' `wt', absorb(`absorb') `reghdfeopt'
		   }
		   
		   * store results
		   tempname poly_b poly_adjw
	       if (_rc==0) {
	 	       matrix `poly_b'=e(b)
			   if (`nwvar'>0&`deriv'==0) matrix `poly_adjw'=`wval'*`poly_b'[1, `=`testpolyp'+1'..`=`testpolyp'+`nwvar'']'
               else                      matrix `poly_adjw'=0                        
               
			   if (`deriv'==0) {
			      if (`testpolyp'>0) matrix `poly_b'=(`poly_b'[1, `=`testpolyp'+`nwvar'+1'], `poly_b'[1,1..`testpolyp'])
				  else               matrix `poly_b'=`poly_b'[1, `=`testpolyp'+`nwvar'+1']
               }
			   else {
			      matrix `poly_b'=`poly_b'[1, `deriv'..`testpolyp']
			   }
			   *if ("`estmethod'"=="qreg") matrix `poly_b'=(`poly_b'[1,colsof(`poly_b')], `poly_b'[1, 1..`testpolyp'])
			   *matrix `poly_b'=`poly_b'[1, `=`deriv'+1'..`=`testpolyp'+1']
	       }
	       else {
	           error  _rc
	   	       exit _rc
           }
		 
		   * Data for derivative
		   tempname polym polym0
		   mata: `polym'=J(`uni_last',0,.)
		   forval i=`deriv'/`testpolyp' {
		      mata: `polym'=(`polym', `uni_grid'[,1]:^(`i'-`deriv')*factorial(`i')/factorial(`i'-`deriv'))
	       }
		   
		   mata: `polym'=`polym'*st_matrix("`poly_b'")':+st_matrix("`poly_adjw'")
		   
		   if (("`estmethod'"=="logit"|"`estmethod'"=="probit")&"`transform'"=="T") {
		      mata: `polym0'=J(rows(`uni_grid'),0,.)
			  if (`deriv'==1) {
			     forval i=1/`testpolyp' {
                    mata: `polym0'=(`polym0', `uni_grid'[,1]:^`i')
                 }
				 if (`nwvar'>0) mata: `polym0'=(`polym0', J(rows(`polym0'),1,1)#st_matrix("`wval'"))
				 mata: `polym0'=(`polym0', J(rows(`polym0'),1,1))
			  }
			  
			  if ("`estmethod'"=="logit") {
			     if (`deriv'==0) mata: `polym'=logistic(`polym')
				 if (`deriv'==1) mata: `polym'=logisticden(`polym0'*st_matrix("e(b)")'):*`polym'
			  }
			  else {
			     if (`deriv'==0) mata: `polym'=normal(`polym')
			     if (`deriv'==1) mata: `polym'=normalden(`polym0'*st_matrix("e(b)")'):*`polym'
			  }
			  mata: mata drop `polym0'
		   }
		   
	   	   if ("`lp'"=="inf") {
		      mata: st_matrix("`stat_poly'", (max(abs((`Xm'[,1]-`polym'):/`Xm'[,2])),3))
		   }
		   else {
		      mata: st_matrix("`stat_poly'", (mean(abs((`Xm'[,1]-`polym'):/`Xm'[,2]):^`lp')^(1/`lp'),3))
		   }
		   
		   * p value
		   if ("`estmethod'"=="qreg"|"`estmethod'"=="reghdfe") {
		      if (`deriv'==0) mata: `uni_basis'=(`uni_basis', J(rows(`uni_basis'),1,1))
		      else            mata: `uni_basis'=(`uni_basis', J(rows(`uni_basis'),1,0))
		      mata: `vcov'=st_matrix("`tmod_V'"); ///
		            `vcov'= (`vcov'[|1,1 \ `nseries', `nseries'|], `vcov'[|1,cols(`vcov') \ `nseries', cols(`vcov')|] \ ///
				             `vcov'[|cols(`vcov'), 1 \ cols(`vcov'), `nseries'|], `vcov'[cols(`vcov'), cols(`vcov')]); ///
			        st_matrix("`vcov'", `vcov')
		   }
		   
		   if ("`estmethod'"!="qreg"&"`estmethod'"!="reghdfe") {
		      mata: `Xm'=binsreg_pred(`uni_basis', ., st_matrix("`tmod_V'")[|1,1 \ `nseries',`nseries'|], "se"); ///
			        binsreg_pval(`uni_basis', `Xm'[,2], "`tmod_V'", "`stat_poly'", ///
			  	                 `nsims', `nseries', ".", 0, "`pval_poly'", ".", "`lp'")
		   }
		   else {
		      mata: `Xm'=binsreg_pred(`uni_basis', ., `vcov', "se"); ///
		            binsreg_pval(`uni_basis', `Xm'[,2], "`vcov'", "`stat_poly'", ///
			  	                 `nsims', `=`nseries'+1', ".", 0, "`pval_poly'", ".", "`lp'")
		   }
							  
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
		   binsreg_irecode `x_var', knotmat(`kmat') bin(`uni_xcat') ///
		                            `usegtools' nbins(`nbins') pos(`binspos') knotliston(T)
		   
		   mata: `uni_basis'=binsreg_spdes(st_data(.,"`x_var'"), "`kmat'", st_data(.,"`uni_xcat'"), ///
		                                  `tmod_p', `deriv', `tmod_s')
										  
		   if (("`estmethod'"=="logit"|"`estmethod'"=="probit")&"`transform'"=="T") {
		      if (`deriv'==0) {
		         mata: `fit0'=(`uni_basis', J(rows(`uni_basis'),1,1)#`wvec0')*st_matrix("`tmod_b'")
			     if ("`estmethod'"=="logit") {
			        mata: `fit'=logistic(`fit0'); ///
				          `se'=logisticden(`fit0'):* ///
					           binsreg_pred((`uni_basis', J(rows(`uni_basis'),1,1)#`wvec'),.,st_matrix("`tmod_V'"),"se")[,2]
			     }
			     else {
			        mata: `fit'=normal(`fit0'); ///
				          `se'=normalden(`fit0'):* ///
					           binsreg_pred((`uni_basis', J(rows(`uni_basis'),1,1)#`wvec'),.,st_matrix("`tmod_V'"),"se")[,2]
			     }
		      }
		      if (`deriv'==1) {
		         mata: `Xm0'=binsreg_spdes(st_data(.,"`x_var'"), "`kmat'", st_data(.,"`uni_xcat'"), `tmod_p', 0, `tmod_s'); ///
			           `Xm0'=(`Xm0', J(rows(`Xm0'),1,1)#`wvec0'); ///
					   `fit0'=`Xm0'*st_matrix("`tmod_b'"); ///
					   `Xm'=(`uni_basis', J(rows(`uni_basis'),1,1)#`wvec')
			     if ("`estmethod'"=="logit") {
			        mata: `fit'=binsreg_pred(`Xm',st_matrix("`tmod_b'"),.,"xb")[,1]
				    if ("`asyvar'"=="off") {
				       mata: `Xm'=logisticden(`fit0'):*(1:-2*logistic(`fit0')):*`fit':*`Xm0' + ///
				                  logisticden(`fit0'):*`Xm'; ///
				             `se'=sqrt(rowsum((`Xm'*st_matrix("`tmod_V'")):*`Xm'))
				    }
				    else {
				       mata: `se'=logisticden(`fit0'):*(binsreg_pred(`Xm',.,st_matrix("`tmod_V'"),"se")[,2])
				    }
				    mata: `fit'=logisticden(`fit0'):*`fit'
			     }
			     else {
			        mata: `fit'=binsreg_pred(`Xm',st_matrix("`tmod_b'"),.,"xb")[,1]
				    if ("`asyvar'"=="off") {
					    mata:`Xm'=(-`fit0'):*normalden(`fit0'):*`fit':*`Xm0' + ///
                                  normalden(`fit0'):*`Xm'; ///
						     `se'=sqrt(rowsum((`Xm'*st_matrix("`tmod_V'")):*`Xm'))
				    }
				    else {
				       mata: `se'=normalden(`fit0'):*(binsreg_pred(`Xm',.,st_matrix("`tmod_V'"),"se")[,2])
				    }
				    mata: `fit'=normalden(`fit0'):*`fit'
			     }
		      }
		      mata: `Xm'=(`fit', `se')
		   }
		   else {
		      mata: `Xm'=(`uni_basis', J(rows(`uni_basis'),1,1)#`wvec'); ///
		            `Xm'=binsreg_pred(`Xm', st_matrix("`tmod_b'"), st_matrix("`tmod_V'"), "all")
		   }
		   		   
		   mata: `tstat'=J(`nfitval',2,.)
		   local counter=1
		   if ("`lp'"=="inf") {
	          foreach var of local varls {
		         mata: `tstat'[`counter',]=(max(abs((`Xm'[,1]-st_data(.,"`var'")):/`Xm'[,2])), 3)
		         local ++counter
		      }
		   }
		   else {
		      foreach var of local varls {
		         mata: `tstat'[`counter',]=(mean(abs((`Xm'[,1]-st_data(.,"`var'")):/`Xm'[,2]):^`lp')^(1/`lp'), 3)
		         local ++counter
		      }
		   }
		   mata: st_matrix("`stat_model'", `tstat')
		   
		   * p values
		   if ("`estmethod'"=="qreg"|"`estmethod'"=="reghdfe") {
		      if (`deriv'==0) mata: `uni_basis'=(`uni_basis', J(rows(`uni_basis'),1,1))
		      else            mata: `uni_basis'=(`uni_basis', J(rows(`uni_basis'),1,0))
		      mata: `vcov'=st_matrix("`tmod_V'"); ///
		            `vcov'= (`vcov'[|1,1 \ `nseries', `nseries'|], `vcov'[|1,cols(`vcov') \ `nseries', cols(`vcov')|] \ ///
				             `vcov'[|cols(`vcov'), 1 \ cols(`vcov'), `nseries'|], `vcov'[cols(`vcov'), cols(`vcov')]); ///
			        st_matrix("`vcov'", `vcov')
		   }
		   if ("`estmethod'"!="qreg"&"`estmethod'"!="reghdfe") {
		      mata: `Xm'=binsreg_pred(`uni_basis', ., st_matrix("`tmod_V'")[|1,1 \ `nseries',`nseries'|], "se"); ///
		            binsreg_pval(`uni_basis', `Xm'[,2], "`tmod_V'", "`stat_model'", `nsims', ///
		                         `nseries', ".", 0, "`pval_model'", ".", "`lp'")
		   }
		   else {
		      mata: `Xm'=binsreg_pred(`uni_basis', ., `vcov', "se"); ///
		            binsreg_pval(`uni_basis', `Xm'[,2], "`vcov'", "`stat_model'", `nsims', ///
		                         `=`nseries'+1', ".", 0, "`pval_model'", ".", "`lp'")
		   }
		   
		   mata: mata drop `Xm' `tstat' `uni_basis'
		}
	 }
	 else {
	    local tmod_p=.
		local tmod_s=.
	 }
	 mata: mata drop `uni_grid' `xvec' `Xm0' `fit' `se' `fit0' `wvec' `wvec0' `vcov'
	 
	 ****** End of testing *****************************************
	 
	 ******************************
	 ******* Display **************
	 ******************************
	 if ("`knotlist'"!=""|"`selection'"=="NA") {
	    local binselectmethod "User-specified"
		local placement "User-specified"
	 }
	 else {
	 	if ("`binsmethod'"=="DPI") local binselectmethod "IMSE-optimal plug-in choice"
	    if ("`binsmethod'"=="ROT") local binselectmethod "IMSE-optimal rule-of-thumb choice"
	    if ("`selection'"=="J") local binselectmethod "`binselectmethod' (select # of bins)"
	    if ("`selection'"=="P") local binselectmethod "`binselectmethod' (select degree and smoothness)"
	    if ("`binspos'"=="ES") local placement "Evenly-spaced"
	    if ("`binspos'"=="QS") local placement "Quantile-spaced"
	 }
	 
	 di ""
	 di in smcl in gr "Hypothesis tests based on binscatter estimates"
	 di in smcl in gr "Estimation method: `estmethod'"
	 di in smcl in gr "Bin selection method: `binselectmethod'"
	 di in smcl in gr "Placement: `placement'"
	 di in smcl in gr "Derivative: `deriv'"
	 di ""
	 di in smcl in gr "{hline 30}{c TT}{hline 15}"
	 di in smcl in gr "{lalign 1:# of observations}"   _col(30) " {c |} " _col(32) as result %7.0f `N'
	 di in smcl in gr "{lalign 1:# of distinct values}"   _col(30) " {c |} " _col(32) as result %7.0f `Ndist'
	 di in smcl in gr "{lalign 1:# of clusters}"   _col(30) " {c |} " _col(32) as result %7.0f `Nclust'	 
	 di in smcl in gr "{hline 30}{c +}{hline 15}"
	 di in smcl in gr "{lalign 1:Bin/Degree selection:}"             _col(30) " {c |} "
*	 if ("`binselectmethod'"=="User-specified") {
*	    di in smcl in gr "{ralign 29:Degree of polynomial}"         _col(30) " {c |} " _col(39) as result %7.0f "."
*	    di in smcl in gr "{ralign 29:# of smoothness constraints}"  _col(30) " {c |} " _col(39) as result %7.0f "."
*	 }
*	 else {
	 di in smcl in gr "{ralign 29:Degree of polynomial}"         _col(30) " {c |} " _col(32) as result %7.0f `binsp'
	 di in smcl in gr "{ralign 29:# of smoothness constraints}"  _col(30) " {c |} " _col(32) as result %7.0f `binss'
*	 }
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
		   if ("`lp'"=="inf") { 
	          di in smcl in gr "H0: mu ="  _col(20) in gr ///
		                       "{c |}" _col(22) "sup |T|"  _col(40) "p value"
		   }
		   else {
		      di in smcl in gr "H0: mu ="  _col(20) in gr ///
		                       "{c |}" _col(22) "L`lp' of T"  _col(40) "p value"
		   }
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
		   if ("`lp'"=="inf") {
		      di in smcl in gr "H0: mu =" _col(20) in gr ///
		                       "{c |}" _col(22) "sup |T|"  _col(40) "p value"
		   }
		   else {
		      di in smcl in gr "H0: mu =" _col(20) in gr ///
		                       "{c |}" _col(22) "L`lp' of T"  _col(40) "p value"
		   }
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
		   if ("`lp'"=="inf") {
	          di in smcl in gr "H0: mu ="  _col(20) in gr ///
		                       "{c |}" _col(22) "sup |T|"  _col(40) "p value"
		   }
		   else {
		      di in smcl in gr "H0: mu ="  _col(20) in gr ///
		                       "{c |}" _col(22) "L`lp' of T"  _col(40) "p value"
		   }
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
	 
	 ereturn scalar imse_var_rot=`imse_var_rot'
	 ereturn scalar imse_bsq_rot=`imse_bsq_rot'
	 ereturn scalar	imse_var_dpi=`imse_var_dpi'
	 ereturn scalar imse_bsq_dpi=`imse_bsq_dpi'
	 
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


*! version 1.3 03-Jul-2023

capture program drop binspwc
program define binspwc, eclass
    version 13
	 
	syntax varlist(min=2 numeric fv ts) [if] [in] [fw aw pw] , by(varname) [deriv(integer 0) at(string asis) nolink ///
	       estmethod(string) estmethodopt(string asis) absorb(string asis) reghdfeopt(string asis) ///
		   pwc(string) testtype(string) lp(string) ///
		   bins(numlist integer max=2 >=0) bynbins(string) binspos(string) ///
		   pselect(numlist integer >=0) sselect(numlist integer >=0) ///
		   binsmethod(string) nbinsrot(string) samebinsby randcut(numlist max=1 >=0 <=1) ///
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
	 
	 if ("`testtype'"=="") {
	    local testtype "2"
	 }
	 
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
	 
	 * Extract options
	 * default vce
	 if ("`vce'"=="") local vce "vce(robust)"
	 local vcetemp: subinstr local vce "vce(" "", all
     local vcetemp: subinstr local vcetemp ")" "", all
	 tokenize "`vcetemp'", parse(", ")
	 if ("`1'"=="cl"|"`1'"=="clu"|"`1'"=="clus"|"`1'"=="clust"| /// 
		 "`1'"=="cluste"|"`1'"=="cluster") {
		if ("`3'"==""|"`3'"==",") local clusterON "T"           /* Mark cluster is specified */
		local clustervar `2'
		if ("`estmethod'"=="qreg") {
		   local vce "vce(robust)"
		   di as text in gr "Warning: vce(cluster) not allowed. vce(robust) used instead."
		}
	 }
	 
	 * use bootstrap cmd?
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
	 else {
		local boot "off"
	 }  
	 
	 if ("`asyvar'"=="") local asyvar "off"

	 
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
	 
	 if ("`binspos'"=="es") local binspos "ES"
	 if ("`binspos'"=="qs") local binspos "QS"
	 if ("`binspos'"=="")   local binspos "QS"
	 if ("`binsmethod'"=="rot") local binsmethod "ROT"
	 if ("`binsmethod'"=="dpi") local binsmethod "DPI"
	 if ("`binsmethod'"=="")    local binsmethod "DPI"

	 * degree, smoothness and binning
	 if ("`pwc'"!="T"&"`pwc'"!="F"&"`pwc'"!="") {
	    numlist "`pwc'", integer max(2) range(>=0)
		local pwc=r(numlist)
	 }
	 
	 if ("`pwc'"=="F") local pwc ""
	 if ("`pwc'"=="T"&"`binspos'"!="ES"&"`binspos'"!="QS") local pwc ""
	 
	 local selection ""
	 
	 * analyze nbins
	 local lenbynbins=0
	 if ("`bynbins'"!="T"&"`bynbins'"!="") {
	    numlist "`bynbins'", integer range(>=0)
	    local bynbins=r(numlist)
		local lenbynbins: word count `bynbins'
	    if (`lenbynbins'==1) {
	       local nbins_all=`bynbins'
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
	 
	 if ("`binspos'"!="ES"&"`binspos'"!="QS") {
	    if ("`bynbins'"!=""|"`pselect'"!=""|"`sselect'"!="") {
           di as error "bynbins(), pselect() or sselect() incorrectly specified."
		   exit
		}
	 }

	 tokenize `bins'
	 local binsp "`1'"
	 local binss "`2'"
	 if ("`binsp'"=="") local binsp=.
	 if ("`binss'"=="") local binss=`binsp'
	 if ("`bins'"!="") {
	 	if ("`bynbins'"!=""&"`bynbins'"!="T"&"`bynbins'"!="0") {
			di as error "bynbins() or bins() is incorrectly specified."
			exit
		}
	 }
	 
	 * 1st case: select J
	 if (("`bins'"!=""|"`bynbins'"=="0"|"`bynbins'"=="T"|"`bynbins'"=="")&("`binspos'"=="ES"|"`binspos'"=="QS")) local selection "J"
	 
	 if ("`selection'"=="J") {
	 	if (`len_p'>1|`len_s'>1) {
		   di as error "Only one p and one s are allowed to select # of bins."
		   exit
		}
		if ("`plist'"=="") local plist=`deriv'
		if ("`slist'"=="") local slist=`plist'
		if ("`bins'"=="") {
		   local binsp `plist'
		   local binss `slist'
		}
		local len_p=1
	    local len_s=1
		if ("`pwc'"=="T"|"`pwc'"=="") local pwc `=`binsp'+1' `=`binss'+1'
     }                                                                          
	 
	 * 2nd case: select P (the special case with nbins() pselect() will be modified in the next step) 
	 if ("`selection'"!="J" & ("`pwc'"=="T"|"`pwc'"=="")) local pselectOK "T"  

	 if (("`pselectOK'"=="T") & ("`bynbins'"!=""&"`bynbins'"!="T") & (`len_p'>1|`len_s'>1)) {
	    local selection "P"
	 }                                                                          
	 
	 * 3rd case: user-specified J and p
	 if ((`len_p'<=1&`len_s'<=1) & "`selection'"!="J") {
		local selection "NA" 
		if ("`pwc'"=="") {
		   if ("`bins'"!="")  local pwc `=`binsp'+1' `=`binss'+1'
		   else {
		   	 if (`len_p'==1&`len_s'==1) local pwc `=`plist'+1' `=`slist'+1'
			 else                       local pwc `=`deriv'+1' `=`deriv'+1'
		   }
		}
	 }
	 
	 * exclude all other cases
	 if ("`selection'"=="") {
	    di as error "Degree, smoothness, or # of bins are not correctly specified."
        exit
	 }
	 
	 if ("`selection'"=="NA") local binselectmethod "User-specified"
	 else {
	    if ("`binsmethod'"=="DPI")  local binselectmethod "IMSE-optimal plug-in choice"
        if ("`binsmethod'"=="ROT")  local binselectmethod "IMSE-optimal rule-of-thumb choice"
		if ("`selection'"=="J") local binselectmethod "`binselectmethod' (select # of bins)"
		if ("`selection'"=="P") local binselectmethod "`binselectmethod' (select degree and smoothness)"
	 }
	 
	 * option for comparison
	 tokenize `pwc'
	 local tsha_p "`1'"
	 local tsha_s "`2'"
	 if ("`tsha_p'"==""|"`tsha_p'"=="T") local tsha_p=.
	 if ("`tsha_s'"=="") local tsha_s `tsha_p'
	 
	 
	 * Add warnings about degrees for estimation and inference
	 if ("`selection'"=="J") {
	    if ("`tsha_p'"!=".") {
		   if (`tsha_p'<=`binsp') {
		      local tsha_p=`binsp'+1
			  local tsha_s=`tsha_p'
			  di as text "Warning: Degree for pwc() has been changed. It must be greater than the degree for bin selection."
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
	 
	 * error check
	 if (`tsha_p'<`tsha_s'|`binsp'<`binss') {
	     di as error "p cannot be smaller than s."
		 exit
	 }
	 if ("`tsha_p'"!="."&"`binsp'"!=".") {
	    if (`tsha_p'<=`binsp') {
	       di as text in gr "Warning: p for testing <= p for bins() not suggested."
	    }
	 }
	 if (`tsha_p'<`deriv') {
	    di as error "p for test cannot be smaller than deriv."
		exit
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
	 markout `touse' `by', strok
	 qui keep if `touse'
     local eN=_N
	 *local nsize=_N     /* # of rows in the original dataset */
	     
	 if ("`usegtools'"==""&("`masspoints'"!="off"|"`binspos'"=="QS")) {
	    if ("`:sortedby'"!="`x_var'") sort `x_var', stable
	 }
	 
	 *************************************************************
	 * Check number of unique byvals & create local storing byvals
	 local byvarname `by'        
     capture confirm numeric variable `by'
     if _rc {
		local bystring "T"
        * generate a numeric version
        tempvar by
        tempname bylabel
        qui egen `by'=group(`byvarname'), lname(`bylabel')
     }
                
     local bylabel `:value label `by'' /* catch value labels for numeric by-vars too */ 
                
     tempname byvalmatrix
     qui tab `by', nofreq matrow(`byvalmatrix')
     
	 * save by-value in a local and calculate group mins and maxs
     local bynum=r(r)
	 if (`bynum'==1) {
	    di as error "More than one group is required."
		exit
	 }
	 tempname xminmat xmaxmat Nmat
	 matrix `xminmat'=J(`bynum', 1, .)
	 matrix `xmaxmat'=`xminmat'
	 matrix `Nmat'=`xminmat'
	 local Ntotal=0
     forvalues i=1/`bynum' {
	    local byv `=`byvalmatrix'[`i',1]' 
        local byvals `byvals' `byv'
	    if ("`wtype'"=="f") sum `x_var' if `by'==`byv' `wt', meanonly
	    else                sum `x_var' if `by'==`byv', meanonly
	    mat `xminmat'[`i',1]=r(min)
	    mat `xmaxmat'[`i',1]=r(max)
	    mat `Nmat'[`i',1]=r(N)        /* sample size, with wt */
		local Ntotal=`Ntotal'+r(N)
     }
	 mata: st_local("Ntotal", strofreal(sum(st_matrix("`Nmat'"))))
	 
	 * define a common support for eval points
	 mata: st_local("max_xmin", strofreal(max(st_matrix("`xminmat'")))); ///
	       st_local("min_xmax", strofreal(min(st_matrix("`xmaxmat'")))); ///
	 	   st_local("xmin", strofreal(min(st_matrix("`xminmat'")))); ///
		   st_local("xmax", strofreal(max(st_matrix("`xmaxmat'"))))
	 
	 
	 * Temp name in MATA
	 tempname xvec yvec byvec cluvec binedges
	 mata: `xvec'=st_data(., "`x_var'"); `yvec'=st_data(.,"`y_var'"); `byvec'=.; `cluvec'=.
	 
	 *******************************************************
	 *** Mass point counting *******************************
	 tempname Ndistlist Nclustlist mat_imse_var_rot mat_imse_bsq_rot mat_imse_var_dpi mat_imse_bsq_dpi
	 mat `Ndistlist'=J(`bynum',1,.)
	 mat `Nclustlist'=J(`bynum',1,.)
	 * Matrices saving imse
	 mat `mat_imse_var_rot'=J(`bynum',1,.)
	 mat `mat_imse_bsq_rot'=J(`bynum',1,.)
	 mat `mat_imse_var_dpi'=J(`bynum',1,.)
	 mat `mat_imse_bsq_dpi'=J(`bynum',1,.)

	 if (`bynum'>1) mata: `byvec'=st_data(.,"`by'")
	 if ("`clusterON'"=="T") mata: `cluvec'=st_data(.,"`clustervar'")
	 
	 ********************************************************
	 ********** Bins, based on FULL sample ******************
	 ********************************************************
	 * knotlist: inner knot seq; knotlistON: local, knot available before loop
	 
	 tempname fullkmat   /* matrix name for saving knots based on the full sample */
	 
	 * Extract user-specified knot list
	 if ("`binspos'"!="ES"&"`binspos'"!="QS") {
		capture numlist "`binspos'", ascending
		if (_rc==0) {
		   local knotlistON "T"
		   local knotlist `binspos'
		   local nbins_all: word count `knotlist'
		   local first: word 1 of `knotlist'
		   local last: word `nbins_all' of `knotlist'
		   if (`first'<=`max_xmin'|`last'>=`min_xmax') {
			  di as error "Inner knots specified out of allowed range."
			  exit
		   }
		   else {
		      local nbins_all=`nbins_all'+1
              local pos "user"
				   
			  foreach el of local knotlist {
		         mat `fullkmat'=(nullmat(`fullkmat') \ `el')
		      }
              mat `fullkmat'=(`xmin' \ `fullkmat' \ `xmax')
		   }
	    }
	    else {
		   di as error "Numeric list incorrectly specified in binspos()."
		   exit
		}
	 }
	 
	 * Bin selection using the whole sample if
	 if ("`selection'"!="NA" & "`samebinsby'"!="") {
	    local selectfullON "T"
	 }
	 
	 if ("`selectfullON'"=="T") {
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
	    * # of clusters
	    local Nclust=.
	    if ("`clusterON'"=="T") {
		   if ("`numclust'"!=""&"`numclust'"!=".") {
			  local Nclust=`numclust'
		   }
		   else {
		      if ("`usegtools'"=="") {
		         mata: st_local("Nclust", strofreal(rows(uniqrows(`cluvec'))))
		      }
		      else {
		         qui gunique `clustervar'
			     local Nclust=r(unique)
		      }
		   }
		   local eN=min(`eN', `Nclust')   /* effective sample size */
	    }
	 
		* Check effective sample size
		if ("`binsp'"==".") local binspcheck=6
		else                local binspcheck=`binsp'
		if ("`nbinsrot'"==""&(`eN'<=`dfcheck_n1'+`binspcheck'+1+`qrot')) {
		    * ROT inavailable, exit
			di as error "Too few observations for bin selection."
			exit
	    }
		else {
			local randcut1k `randcut'
		    if ("`randcut'"=="" & `Ntotal'>5000) {
	           local randcut1k=max(5000/`Ntotal', 0.01)
		       di as text in gr "Warning: To speed up computation, bin/degree selection uses a subsample of roughly max(5000, 0.01n) observations if n>5000. To use the full sample, set randcut(1)."
	        }
			if ("`selection'"=="J") {
		      qui binsregselect `y_var' `x_var' `w_var' `wt', deriv(`deriv') bins(`binsp' `binss') nbins() ///
				   			    absorb(`absorb') reghdfeopt(`reghdfeopt') ///
								binsmethod(`binsmethod') binspos(`binspos') nbinsrot(`nbinsrot') ///
								`vce_select' masspoints(`masspoints') dfcheck(`dfcheck_n1' `dfcheck_n2') ///
								numdist(`Ndist') numclust(`Nclust') randcut(`randcut1k') usegtools(`sel_gtools')
			  if (e(nbinsrot_regul)==.) {
			      di as error "Bin selection fails."
				  exit
			  }
			  if ("`binsmethod'"=="ROT") {
				  local nbins_all=e(nbinsrot_regul)
				  mat `mat_imse_var_rot'=J(`bynum',1,e(imse_var_rot))
	              mat `mat_imse_bsq_rot'=J(`bynum',1,e(imse_bsq_rot))
			  }
			  else if ("`binsmethod'"=="DPI") {
				 local nbins_all=e(nbinsdpi)
                 mat `mat_imse_var_dpi'=J(`bynum',1,e(imse_var_dpi))
	             mat `mat_imse_bsq_dpi'=J(`bynum',1,e(imse_bsq_dpi))
				 if (`nbins_all'==.) {
				    local nbins_all=e(nbinsrot_regul)
					mat `mat_imse_var_rot'=J(`bynum',1,e(imse_var_rot))
	                mat `mat_imse_bsq_rot'=J(`bynum',1,e(imse_bsq_rot))
					di as text in gr "Warning: DPI selection fails. ROT choice used."
				 }
			  }
		   }
		   else if ("`selection'"=="P") {
		   	  qui binsregselect `y_var' `x_var' `w_var' `wt', deriv(`deriv') nbins(`nbins_all') ///
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
				  mat `mat_imse_var_rot'=J(`bynum',1,e(imse_var_rot))
	              mat `mat_imse_bsq_rot'=J(`bynum',1,e(imse_bsq_rot))
			  }
			  else if ("`binsmethod'"=="DPI") {
				 local binsp=e(pdpi)
				 local binss=e(sdpi)
				 mat `mat_imse_var_dpi'=J(`bynum',1,e(imse_var_dpi))
	             mat `mat_imse_bsq_dpi'=J(`bynum',1,e(imse_bsq_dpi))
				 if (`binsp'==.) {
				    local binsp=e(prot_regul)
					local binss=e(srot_regul)
					mat `mat_imse_var_rot'=J(`bynum',1,e(imse_var_rot))
	                mat `mat_imse_bsq_rot'=J(`bynum',1,e(imse_bsq_rot))
					di as text in gr "Warning: DPI selection fails. ROT choice used."
				 }
			  }
			  local tsha_p=`binsp'+1
			  local tsha_s=`binss'+1
		   }
		}
	 }
	 
	 if ("`selectfullON'"=="T"|("`selection'"=="NA"&"`samebinsby'"!="")) {
		 * Save in a knot list
         local knotlistON "T"
		 if ("`binspos'"=="ES") {
	        local stepsize=(`xmax'-`xmin')/`nbins_all'
	        forvalues i=1/`=`nbins_all'+1' {
		       mat `fullkmat'=(nullmat(`fullkmat') \ `=`xmin'+`stepsize'*(`i'-1)')
		    }
	     }
	     else if ("`binspos'"=="QS") {
			 if (`nbins_all'==1)  mat `fullkmat'=(`xmin' \ `xmax')
		     else {		
	           binsreg_pctile `x_var' `wt', nq(`nbins_all') `usegtools'
		       mat `fullkmat'=(`xmin' \ r(Q) \ `xmax')
		     }
	     }
	 }
	 
	 *** Placement name, for display ************
	 if ("`pos'"=="user") {
	     local binselectmethod "User-specified"
		 local placement "User-specified"
	 }
	 else if ("`binspos'"=="ES") {
	     local placement "Evenly-spaced"
	 }
	 else if ("`binspos'"=="QS") {
	     local placement "Quantile-spaced"
	 }
	 
	 
	 ********************************************************
	 * Set seed
	 if ("`simsseed'"!="") set seed `simsseed'
	 
	 * generate eval points
	 tempname Xm uni_grid uni_grid_bin uni_basis num denom nummat tstat pmat xsub ysub byindex xcatsub   /* objects in MATA */
	 mata: `tstat'=J(`=`bynum'*(`bynum'-1)/2',3,.); `pmat'=J(`=`bynum'*(`bynum'-1)/2',1,.)
	 
	 tempname Xm0 fit fit0 se vcov
	 mata: `Xm0'=.; `fit'=.; `fit0'=0; `se'=.; `vcov'=.
	 
	 tempvar xcat bycond
	 qui gen `xcat'=. in 1 
	 qui gen `bycond'=. in 1
	 
	 * matrix names, for returns
	 tempname nbinslist teststat pvalue pwc_plist pwc_slist
	 
	 * prepare grid
	 mata: `uni_grid'=rangen(`max_xmin', `min_xmax', `simsgrid'+2); ///
	       `uni_grid'=`uni_grid'[|2 \ `=`simsgrid'+1'|]           /* only keep inner points, simsgrid>=1 */
	 
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
		        if ("`wtype'"!="") qui tabstat `wname' [aw`exp'], stat(`at') save
			    else               qui tabstat `wname', stat(`at') save
			    mat `wvaltemp'=r(StatTotal)
			 }
			 else {
			    qui gstats tabstat `wname' `wt', stat(`at') matasave("`mataobj'")
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
	 

	local byvalnamelist ""              /* save group name (value) */
	local counter=1
	local counter2=1
	***************************************************************************
	******************* Now, enter the loop ***********************************
	***************************************************************************
	foreach byval in `byvals' {
		local conds "if `by'==`byval'"     /* with "if" */
		qui replace `bycond'=(`by'==`byval')
		
	    if ("`bylabel'"=="") local byvalname=`byval'
		else {
		   local byvalname `: label `bylabel' `byval''
		}
		local byvalnamelist `" `byvalnamelist' `"`byvalname'"' "'
		
		mata: `byindex'=`byvec':==`byval'
		mata: `xsub'=select(`xvec',`byindex'); `ysub'=select(`yvec', `byindex')
		
		************************************
		* Calculate various sample sizes
		* Subsample size
		if ("`wtype'"=="f") sum `x_var' `conds' `wt', meanonly
		else                sum `x_var' `conds', meanonly
		
		local xmin=r(min)
	    local xmax=r(max)
		local N=r(N)
		
	    * Effective sample size
		if ("`wtype'"!="f") local eN=r(N)
		else {
		   qui count `conds'
	       local eN=r(N)
		}
		
	    local Ndist=.
	    if ("`massadj'"=="T") {
		   if ("`usegtools'"=="") {
		      mata: `binedges'=binsreg_uniq(`xsub', ., 1, "Ndist")
		      mata: mata drop `binedges'
		   }
		   else {
		      qui gunique `x_var' `conds'
			  local Ndist=r(unique)
		   }
		   local eN=min(`eN', `Ndist')
		   mat `Ndistlist'[`counter',1]=`Ndist'
	    }
		
	    * # of clusters
	    local Nclust=.
	    if ("`clusterON'"=="T") {
		   if ("`usegtools'"=="") {
		      mata: st_local("Nclust", strofreal(rows(uniqrows(select(`cluvec', `byindex')))))
		   }
		   else {
			  qui gunique `clustervar' `conds'
			  local Nclust=r(unique)
		   }
		   local eN=min(`eN', `Nclust')   /* effective SUBsample size */
		   mat `Nclustlist'[`counter',1]=`Nclust'
	    }
	    
		*********************************************************
		************** Prepare bins, within loop ****************
		*********************************************************
		if ("`pos'"!="user") local pos `binspos'              /* initialize pos */
		
		* Selection?
		local nbins ""
		if (`lenbynbins'>1)    local nbins: word `counter' of `bynbins'
		if ("`nbins_all'"!="") local nbins=`nbins_all'    /* add the universal nbins */
		
	    if ("`selection'"!="NA"&"`knotlistON'"!="T") {
		   * Check effective sample size
		   if ("`binsp'"==".") local binspcheck=6
		   else                local binspcheck=`binsp'
		   if ("`nbinsrot'"==""&(`eN'<=`dfcheck_n1'+`binspcheck'+1+`qrot')) {
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
			     qui binsregselect `y_var' `x_var' `w_var' `conds' `wt', deriv(`deriv') ///
				                   bins(`binsp' `binss') nbins() ///
		                           absorb(`absorb') reghdfeopt(`reghdfeopt') /// 
					   			   binsmethod(`binsmethod') binspos(`pos') nbinsrot(`nbinsrot') ///
							       `vce_select' masspoints(`masspoints') dfcheck(`dfcheck_n1' `dfcheck_n2') ///
							       numdist(`Ndist') numclust(`Nclust') randcut(`randcut1k') usegtools(`sel_gtools')
			     if (e(nbinsrot_regul)==.) {
			        di as error "Bin selection fails."
			        exit
			     }
			     if ("`binsmethod'"=="ROT") {
		            local nbins=e(nbinsrot_regul)
				    mat `mat_imse_bsq_rot'[`counter',1]=e(imse_bsq_rot)
				    mat `mat_imse_var_rot'[`counter',1]=e(imse_var_rot)
		         }
		         else if ("`binsmethod'"=="DPI") {
			        local nbins=e(nbinsdpi)
				    mat `mat_imse_bsq_dpi'[`counter',1]=e(imse_bsq_dpi)
				    mat `mat_imse_var_dpi'[`counter',1]=e(imse_var_dpi)
			        if (`nbins'==.) {
			           local nbins=e(nbinsrot_regul)
				       mat `mat_imse_bsq_rot'[`counter',1]=e(imse_bsq_rot)
				       mat `mat_imse_var_rot'[`counter',1]=e(imse_var_rot)
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
					mat `mat_imse_bsq_rot'[`counter',1]=e(imse_bsq_rot)
				    mat `mat_imse_var_rot'[`counter',1]=e(imse_var_rot)
			     }
			     else if ("`binsmethod'"=="DPI") {
				    local binsp=e(pdpi)
				    local binss=e(sdpi)
                    mat `mat_imse_bsq_dpi'[`counter',1]=e(imse_bsq_dpi)
				    mat `mat_imse_var_dpi'[`counter',1]=e(imse_var_dpi)
				    if (`binsp'==.) {
				       local binsp=e(prot_regul)
					   local binss=e(srot_regul)
                       mat `mat_imse_bsq_rot'[`counter',1]=e(imse_bsq_rot)
				       mat `mat_imse_var_rot'[`counter',1]=e(imse_var_rot)
					   di as text in gr "Warning: DPI selection fails. ROT choice used."
				    }
			     }
				 local tsha_p=`binsp'+1
				 local tsha_s=`binss'+1
		      }				  
		   }
	    }
		
		*******************************************************
	    * Check if eff. sample size is large enough for testing
		if ((`nbins'-1)*(`tsha_p'-`tsha_s'+1)+`tsha_p'+1+`dfcheck_n2'>=`eN') {
		   di as text in gr "Warning: Too small effective sample size for testing shape."
	    }
		
	    * Generate category variable for data and save knot in matrix
		tempname kmat
		if ("`knotlistON'"=="T") {
		   mat `kmat'=`fullkmat'
		}
		else {
           if ("`pos'"=="ES") {
	          local stepsize=(`xmax'-`xmin')/`nbins'
	          forvalues i=1/`=`nbins'+1' {
		         mat `kmat'=(nullmat(`kmat') \ `=`xmin'+`stepsize'*(`i'-1)')
		      }
		   }
	       else {		
		      if (`nbins'==1)  mat `kmat'=(`xmin' \ `xmax')
		      else {		
	             binsreg_pctile `x_var' `conds' `wt', nq(`nbins') `usegtools'
		         mat `kmat'=(`xmin' \ r(Q) \ `xmax')
		      }
		   }
		}
	    
		* Renew knot list
		mata: st_matrix("`kmat'", (`xmin' \ uniqrows(st_matrix("`kmat'")[|2 \ `=`nbins'+1'|])))
	    if (`nbins'!=rowsof(`kmat')-1) {
	       di as text in gr "Warning: Repeated knots. Some bins dropped."
		   local nbins=rowsof(`kmat')-1
	    }
		binsreg_irecode `x_var' `conds', knotmat(`kmat') bin(`xcat') ///
		                                 `usegtools' nbins(`nbins') pos(`pos') knotliston(`knotlistON')
		
		mata: `xcatsub'=st_data(., "`xcat'")
		mata: `xcatsub'=select(`xcatsub', `byindex')
        
		* Now, save nbins, p and s in a matrix !!!
		mat `nbinslist'=(nullmat(`nbinslist') \ `nbins')
		mat `pwc_plist'=(nullmat(`pwc_plist') \ `tsha_p')
		mat `pwc_slist'=(nullmat(`pwc_slist') \ `tsha_s')
		
		* Check for empty bins
	    if ("`localcheck'"=="T") {
		   mata: st_local("Ncat", strofreal(rows(uniqrows(`xcatsub'))))
		   if (`nbins'==`Ncat') {
		      mata: `binedges'=binsreg_uniq(`xsub', `xcatsub', `nbins', "uniqmin")
			  mata: mata drop `binedges'
		   }
		   else {
		      local uniqmin=0
			  di as text in gr "Warning: There are empty bins. Specify a smaller number in nbins()."
		   }
	       
		   if (`uniqmin'<`tsha_p'+1) {
		      di as text in gr "Warning: Some bins have too few distinct x-values for testing."
		   }
	    }		
		
		************************************************************
        ************************************************************
	    * Regression
	    local nseries=(`tsha_p'-`tsha_s'+1)*(`nbins'-1)+`tsha_p'+1
	    local tsha_series ""
	    forvalues i=1/`nseries' {
	       tempvar sp`i'
	       local tsha_series `tsha_series' `sp`i''
		   qui gen `sp`i''=. in 1
	    }
		
		tempname tsha_b tsha_V
		mata: binsreg_st_spdes(`xsub', "`tsha_series'", "`kmat'", `xcatsub', `tsha_p', 0, `tsha_s', "`bycond'")
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
		mata: `uni_grid_bin'`counter'=binspwc_locate(`uni_grid', st_matrix("`kmat'"))
		
	    * fitted values
		mata: `uni_basis'=binsreg_spdes(`uni_grid', "`kmat'", `uni_grid_bin'`counter', `tsha_p', `deriv', `tsha_s')
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
		      mata: `Xm0'=binsreg_spdes(`uni_grid', "`kmat'", `uni_grid_bin'`counter', `tsha_p', 0, `tsha_s'); ///
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
		

		* num: fitted value; denom: standard error
		mata: `num'`counter'=`Xm'[,1]; ///
		      `denom'`counter'=`Xm'[,2]
			    
		* For p value
		if ("`estmethod'"=="qreg"|"`estmethod'"=="reghdfe") {
		   if (`deriv'==0) mata: `uni_basis'=(`uni_basis', J(rows(`uni_basis'),1,1))
		   else            mata: `uni_basis'=(`uni_basis', J(rows(`uni_basis'),1,0))
		   mata: `vcov'=st_matrix("`tsha_V'"); ///
		         `vcov'= (`vcov'[|1,1 \ `nseries', `nseries'|], `vcov'[|1,cols(`vcov') \ `nseries', cols(`vcov')|] \ ///
				          `vcov'[|cols(`vcov'), 1 \ cols(`vcov'), `nseries'|], `vcov'[cols(`vcov'), cols(`vcov')]); ///
			     st_matrix("`vcov'", `vcov')
		}
		
		if ("`estmethod'"!="qreg"&"`estmethod'"!="reghdfe") {
		   mata: `se'`counter'=binsreg_pred(`uni_basis', ., st_matrix("`tsha_V'")[|1,1 \ `nseries',`nseries'|], "se")[,2]; ///
		         `nummat'`counter'=binspwc_nummat(`uni_basis', "`tsha_V'", `nseries')
		}
		else {
		   mata: `se'`counter'=binsreg_pred(`uni_basis', ., `vcov', "se")[,2]; ///
		         `nummat'`counter'=binspwc_nummat(`uni_basis', "`vcov'", `=`nseries'+1')
		}
				
		
		* pairwise comparison
		if (`counter'>1) {
		   forval gr=1/`=`counter'-1' {
		      * calculate test stat
			  if ("`testtype'"=="l") {
			     mata: `tstat'[`counter2',.]=(max((`num'`counter'-`num'`gr'):/ ///
				                                   sqrt((`denom'`counter':^2)+(`denom'`gr':^2))), `counter', `gr')
			  }
			  else if ("`testtype'"=="r") {
			  	 mata: `tstat'[`counter2',.]=(min((`num'`counter'-`num'`gr'):/ ///
				                                   sqrt((`denom'`counter':^2)+(`denom'`gr':^2))), `counter', `gr')
			  }
			  else {
			     if ("`lp'"=="inf") {
			     mata: `tstat'[`counter2',.]=(max(abs((`num'`counter'-`num'`gr'):/ ///
				                                       sqrt((`denom'`counter':^2)+(`denom'`gr':^2)))), `counter', `gr')
				 }
				 else {
				 mata: `tstat'[`counter2',.]=(mean(((`num'`counter'-`num'`gr'):/ ///
				                                   sqrt((`denom'`counter':^2)+(`denom'`gr':^2))):^`lp')^(1/`lp'), `counter', `gr')
				 }
			  }
			  
			  * calculate p val
			  mata: `pmat'[`counter2',1]=binspwc_pval(`nummat'`counter', `nummat'`gr', `se'`counter', `se'`gr', ///
			                                             `tstat'[`counter2',1], `nsims', "`testtype'", "`lp'")
									
			  local ++counter2
		   }
		}
		
		drop `tsha_series'
	 
	    local ++counter
	 
	 }
	 
	 mata: st_matrix("`teststat'", `tstat'); st_matrix("`pvalue'", `pmat')
	 
	 * drop objects in MATA
	 mata: mata drop `Xm' `uni_grid' `uni_basis' `tstat' `pmat' `xsub' `ysub' `byindex' `xcatsub' ///
	                 `xvec' `yvec' `byvec' `cluvec' `Xm0' `fit' `fit0' `se' `vcov' `wvec' `wvec0'
	
	 mata: mata drop `uni_grid_bin'* `num'* `denom'* `nummat'* `se'*
	 
	 
	 ******************************
	 ******* Display **************
	 ******************************
	 di ""
	 di in smcl in gr "Pairwise group comparison based on binscatter estimates"
	 di in smcl in gr "Estimation method: `estmethod'"
	 di in smcl in gr "Derivative: `deriv'"
	 di in smcl in gr "Group variable: `byvarname'"
	 di in smcl in gr "Bin/Degree selection method: `binselectmethod'"
	 di in smcl in gr "Placement: `placement'"

	 *di in smcl in gr "{hline 30}{c TT}{hline 15}"
	 *di in smcl in gr "{lalign 1:Bin/Degree selection:}"            _col(30) " {c |} "
	 *di in smcl in gr "{ralign 29:Degree of polynomial}"         _col(30) " {c |} " _col(32) as result %7.0f `binsp'
	 *di in smcl in gr "{ralign 29:# of smoothness constraints}"  _col(30) " {c |} " _col(32) as result %7.0f `binss'
	 *di in smcl in gr "{hline 30}{c +}{hline 15}"
	 *di in smcl in gr "{lalign 1:Hypothesis test:}"             _col(30) " {c |} "
	 *di in smcl in gr "Pairwise Group Comparison, by `byvarname':"
	 *di in smcl in gr "{ralign 29:Degree of polynomial}"         _col(30) " {c |} " _col(32) as result %7.0f `tsha_p'
	 *di in smcl in gr "{ralign 29:# of smoothness constraints}"  _col(30) " {c |} " _col(32) as result %7.0f `tsha_s'
     *di in smcl in gr "{hline 30}{c BT}{hline 15}"

	 
     forval i=1/`=`counter2'-1' {
	    local g1=`=`teststat'[`i',2]'
		local g2=`=`teststat'[`i',3]'
	    local group1: word `g1' of `byvalnamelist'
		local group2: word `g2' of `byvalnamelist'
	    
		di ""
		di in smcl in gr "Group `group1' vs. Group `group2'"
		di in smcl in gr "{hline 30}{c TT}{hline 10}{c TT}{hline 10}" 
		di in smcl in gr "{lalign 1:Group `byvarname'=}"   _col(30) " {c |}" _col(38) "`group1'"  _col(42) "{c |}" _col(48) "`group2'"
	    di in smcl in gr "{hline 30}{c +}{hline 10}{c +}{hline 10}"
	    di in smcl in gr "{lalign 1:# of observations}"   _col(30) " {c |} " _col(32) as result %7.0f `=`Nmat'[`g1',1]' _col(42) in gr "{c |}" _col(44) as result %7.0f `=`Nmat'[`g2',1]'
	    di in smcl in gr "{lalign 1:# of distinct values}"   _col(30) " {c |} " _col(32) as result %7.0f `=`Ndistlist'[`g1',1]' _col(42) in gr "{c |}" _col(44) as result %7.0f `=`Ndistlist'[`g2',1]'
	    di in smcl in gr "{lalign 1:# of clusters}"   _col(30) " {c |} " _col(32) as result %7.0f `=`Nclustlist'[`g1',1]' _col(42) in gr "{c |}" _col(44) as result %7.0f `=`Nclustlist'[`g2',1]'	 
	   	di in smcl in gr "{lalign 1:Degree of polynomial}"         _col(30) " {c |} " _col(32) as result %7.0f `=`pwc_plist'[`g1',1]' _col(42) in gr "{c |}" _col(44) as result %7.0f `=`pwc_plist'[`g2',1]'
	    di in smcl in gr "{lalign 1:# of smoothness constraints}"  _col(30) " {c |} " _col(32) as result %7.0f `=`pwc_slist'[`g1',1]' _col(42) in gr "{c |}" _col(44) as result %7.0f `=`pwc_slist'[`g2',1]'
	    di in smcl in gr "{lalign 1:# of bins}"       _col(30) " {c |} " _col(32) as result %7.0f `=`nbinslist'[`g1',1]' _col(42) in gr "{c |}" _col(44) as result %7.0f `=`nbinslist'[`g2',1]'
	    di in smcl in gr "{hline 30}{c BT}{hline 21}"

		di ""
		di in smcl in gr "diff = group `group1' - group `group2'"
		if ("`testtype'"=="l") {   
	       di in smcl in gr "{hline 19}{c TT}{hline 30}"
	       di in smcl in gr "H0:"  _col(20) in gr ///
		                    "{c |}" _col(22) "sup T"  _col(40) "p value"
		   di in smcl in gr "{hline 19}{c +}{hline 30}"
		   local stat=`teststat'[`i',1]
		   local pval=`pvalue'[`i',1]
	       di in smcl in gr "diff<=0" _col(20) in gr "{c |}" ///
		                    _col(22) as result %7.3f `stat' ///
		                    _col(40) as result %7.3f `pval'
	       
		   di in smcl in gr "{hline 19}{c BT}{hline 30}"
	    }
	    else if ("`testtype'"=="r") {
	       di in smcl in gr "{hline 19}{c TT}{hline 30}"
		   di in smcl in gr "H0:"  _col(20) in gr ///
		                    "{c |}" _col(22) "inf T"  _col(40) "p value"
		   di in smcl in gr "{hline 19}{c +}{hline 30}"	    
		   local stat=`teststat'[`i',1]
		   local pval=`pvalue'[`i',1]
	       di in smcl in gr "diff>=0" _col(20) in gr "{c |}" ///
		                    _col(22) as result %7.3f `stat' ///
		                    _col(40) as result %7.3f `pval'
		   di in smcl in gr "{hline 19}{c BT}{hline 30}"
	    }
		else {
	       di in smcl in gr "{hline 19}{c TT}{hline 30}"
		   if ("`lp'"=="inf") { 
	          di in smcl in gr "H0:"  _col(20) in gr ///
		                       "{c |}" _col(22) "sup |T|"  _col(40) "p value"
		   }
		   else {
		      di in smcl in gr "H0:"  _col(20) in gr ///
		                       "{c |}" _col(22) "L`lp' of T"  _col(40) "p value"
		   }
		   di in smcl in gr "{hline 19}{c +}{hline 30}"
		   local stat=`teststat'[`i',1]
		   local pval=`pvalue'[`i',1]
	       di in smcl in gr "diff=0" _col(20) in gr "{c |}" ///
		                    _col(22) as result %7.3f `stat' ///
		                    _col(40) as result %7.3f `pval' 
		   di in smcl in gr "{hline 19}{c BT}{hline 30}"
	    }
		di ""
	 }
	 
	 *******************************************************
	 	 
	 ********** Return ***************
	 *********************************
	 ereturn clear
	 * # of observations
	 ereturn scalar N=`Ntotal'
	 ereturn scalar p=`binsp'
	 ereturn scalar s=`binss'
	 
	 * by pair
	 ereturn matrix pval=`pvalue'
	 ereturn matrix stat=`teststat'

	 * by group:
	 ereturn matrix nbins_by=`nbinslist'
	 ereturn matrix pwc_plist=`pwc_plist'
	 ereturn matrix pwc_slist=`pwc_slist'
	 ereturn matrix Nclust_by=`Nclustlist'    
	 ereturn matrix Ndist_by=`Ndistlist'
	 ereturn matrix N_by=`Nmat'
	 
	 ereturn matrix imse_var_rot=`mat_imse_var_rot'
	 ereturn matrix imse_bsq_rot=`mat_imse_bsq_rot'
	 ereturn matrix	imse_var_dpi=`mat_imse_var_dpi'
	 ereturn matrix imse_bsq_dpi=`mat_imse_bsq_dpi'

	 
	 * local: corresponding by-values
	 ereturn local byvalue `byvalnamelist'
end
	 
mata:
   // calculate numerator matrix for simulation
   real matrix binspwc_nummat(real matrix X, string scalar covname, real scalar k)
   {
     real matrix cov, num, U, V, sv
	 
	 cov=st_matrix(covname)[|1,1\k,k|]
	 if (rank(cov)==k) {
	     num=X*cholesky(cov)
	 }
	 else {
	     svd(cov, U=., sv=., V=.)
		 pragma unused V
		 num=X*U*diag(sv:^0.5)*U'
	 }
	 
	 return(num)
   }
 
 
   // calculate p val for pairwise comparison
   real scalar binspwc_pval(real matrix nummatA, real matrix nummatB, /// 
                               real vector denomA, real vector denomB,  
					           real scalar stat, real scalar rep, ///
					           string scalar type, string scalar metric)
   {
     real scalar kA, kB, i, lp, pval
	 real vector t
	 
	 kA=cols(nummatA)
	 kB=cols(nummatB)
	 if (metric!="inf") {
	    lp=strtoreal(metric)
	 }
	 pval=0
   
   	 for (i=1; i<=rep; i++) {
	     t=(nummatA*rnormal(kA,1,0,1)-nummatB*rnormal(kB,1,0,1)):/sqrt(denomA:^2+denomB:^2)
		 if (type=="l") {
		    pval=pval+(max(t)>=stat)
	     }
		 else if (type=="r") {
   		    pval=pval+(min(t)<=stat)
		 }
		 else {
		    if (metric=="inf") {
			   pval=pval+(max(abs(t))>=stat)
		    }
			else {
			   pval=pval+(mean(abs(t):^lp)^(1/lp)>=stat)
			}
		 }
	   }
	   
	   return(pval/rep)
   }
 

   // locate a point 
   real vector binspwc_locate(real vector x, real vector kmat)
   {
     real vector bin, index
	 real scalar n, nbin, i, lb, ub
	 
	 n=rows(x)
	 nbin=rows(kmat)-1
	 bin=J(n,1,.)
	 for (i=1; i<=nbin; i++) {
	    lb=kmat[i,1]
		ub=kmat[i+1,1]
		if (i<nbin) {
		   index=selectindex((x:>=lb):&(x:<ub))
	       bin[index]=J(rows(index),1,i)
		}
		else {
		   index=selectindex((x:>=lb):&(x:<=ub))
		   bin[index]=J(rows(index),1,i)
		}
	 }
     
	 return(bin)
   }
   	
end

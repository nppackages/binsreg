*! version 0.2 13-MAR-2019 

capture program drop binsreg
program define binsreg, eclass
     version 13
	 
	 syntax varlist(min=2 numeric fv ts) [if] [in] [fw aw pw] [, deriv(integer 0) ///
	        dots(numlist integer max=2 >=0) dotsgrid(string) dotsplotopt(string asis) ///
			line(numlist integer max=2 >=0) linegrid(integer 20) lineplotopt(string asis) ///
			ci(numlist integer max=2 >=0) cigrid(string) ciplotopt(string asis) /// 
			cb(numlist integer max=2 >=0) cbgrid(integer 20) cbplotopt(string asis) ///
			polyreg(string) polyreggrid(integer 20) polyregcigrid(integer 0) polyregplotopt(string asis) ///
			by(varname) bycolors(string asis) bysymbols(string asis) bylpatterns(string asis) ///
			testmodel(numlist integer max=2 >=0) testmodelparfit(string asis) testmodelpoly(string) /// 
			testshape(numlist integer max=2 >=0) testshapel(numlist) ///
			testshaper(numlist) testshape2(numlist) ///
			nbins(integer 0) binspos(string) binsmethod(string) nbinsrot(string) samebinsby ///
			nsims(integer 500) simsgrid(integer 20) simsseed(integer 666) ///
			dfcheck(numlist integer max=2 >=0) masspoints(string) ///
			vce(passthru) level(real 95)   ///
			noplot savedata(string asis) replace *]
	 
	 *********************************************
	 * Regularization constant (for checking only)
	 local qrot=2
	 
	 **************************************
	 * Create weight local
     if ("`weight'"!="") {
	    local wt [`weight'`exp']
		local wtype=substr("`weight'",1,1)
	 }
	 
	 **********************
	 ** Extract options ***
	 **********************	 	 
	 * default vce, clustered?
	 if ("`vce'"=="") local vce "vce(robust)"
	 local vcetemp: subinstr local vce "vce(" "", all
     local vcetemp: subinstr local vcetemp ")" "", all
	 tokenize `vcetemp'
	 if ("`1'"=="cl"|"`1'"=="clu"|"`1'"=="clus"|"`1'"=="clust"| /// 
		 "`1'"=="cluste"|"`1'"=="cluster") {
		local clusterON "T"           /* Mark cluster is specified */
		local clustervar `2'
	 }
	 
	 * Extract p and s, grid, set default
	 * dots
	 tokenize `dots'
	 local dots_p "`1'"
	 local dots_s "`2'"
	 if ("`dots_p'"=="") local dots_p 0
	 if ("`dots_s'"=="") local dots_s `dots_p'
	 if ("`dotsgrid'"=="") local dotsgrid "mean"
	 local dotsngrid_mean=0
	 if (strpos("`dotsgrid'","mean")!=0) {
	    local dotsngrid_mean=1
		local dotsgrid: subinstr local dotsgrid "mean" "", all
	 }
	 if (wordcount("`dotsgrid'")==0) local dotsngrid=0
	 else {
	    confirm integer n `dotsgrid'
	    local dotsngrid `dotsgrid'
     }
	 local dotsntot=`dotsngrid_mean'+`dotsngrid'
	 
	 * line
	 tokenize `line'
	 local line_p "`1'"
	 local line_s "`2'"
	 local linengrid `linegrid'
	 if ("`line_p'"=="") {
	    local linengrid=0
		local line_p=.
	 }
	 if ("`line_p'"!=""&"`line_s'"=="") local line_s `line_p'
	 
	 * ci
	 if ("`cigrid'"=="") local cigrid "mean"
	 local cingrid_mean=0
	 if (strpos("`cigrid'","mean")!=0) {
	    local cingrid_mean=1
		local cigrid: subinstr local cigrid "mean" "", all
	 }
	 if (wordcount("`cigrid'")==0) local cingrid=0
	 else {
	    confirm integer n `cigrid'
	    local cingrid `cigrid'
	 }
	 local cintot=`cingrid_mean'+`cingrid'
	 
	 tokenize `ci'
	 local ci_p "`1'"
	 local ci_s "`2'"
	 if ("`ci_p'"=="") {
	    local cintot=0
	    local ci_p=.
	 }
	 if ("`ci_p'"!=""&"`ci_s'"=="") local ci_s `ci_p'

	 
	 * cb
	 tokenize `cb'
	 local cb_p "`1'"
	 local cb_s "`2'"
	 local cbngrid `cbgrid'
	 if ("`cb_p'"=="") {
	    local cbngrid=0
	    local cb_p=.
	 }
	 if ("`cb_p'"!=""&"`cb_s'"=="") local cb_s `cb_p'
	 
	 * poly fit
	 local polyregngrid `polyreggrid'
	 local polyregcingrid `polyregcigrid'
	 if ("`polyreg'"!="") {
	    confirm integer n `polyreg'
	 }
	 else {
	    local polyregngrid=0
	 }
	 
	 * Simuls
	 local simsngrid=`simsgrid'
	 
	 * Options for testing shape
	 tokenize `testshape'
	 local tsha_p "`1'"
	 local tsha_s "`2'"
	 if ("`tsha_p'"=="") local tsha_p 3
	 if ("`tsha_s'"=="") local tsha_s `tsha_p'
	 local nL: word count `testshapel'
	 local nR: word count `testshaper'
	 local nT: word count `testshape2'
	 local ntestshape=`nL'+`nR'+`nT'     /* number of tests (for shape) */
	 
	 * Options for testing model
	 tokenize `testmodel'
	 local tmod_p "`1'"
	 local tmod_s "`2'"
	 if ("`tmod_p'"=="") local tmod_p 3
	 if ("`tmod_s'"=="") local tmod_s `tmod_p'
	 
	 * Record if nbins specified by users, set default
	 local nbins_all=`nbins'              /* local save common nbins */
	 if (`nbins'!=0) local binselectmethod "User-specified"
	 if ("`binsmethod'"=="rot") local binsmethod "ROT"
	 if ("`binsmethod'"=="dpi") local binsmethod "DPI"
	 if ("`binsmethod'"=="")    local binsmethod "DPI"
	 if ("`binspos'"=="es") local binspos "ES"
	 if ("`binspos'"=="qs") local binspos "QS"
	 if ("`binspos'"=="")   local binspos "QS"
	 
	 * Mass point check?
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
	    local fewmasspoints "T"    /* count mass point, but turn off checks */
	 }
	 
	 * extract dfcheck
	 if ("`dfcheck'"=="") local dfcheck 20 30
	 tokenize `dfcheck'
	 local dfcheck_n1 "`1'"
	 local dfcheck_n2 "`2'"
     
	 
	 *************************
	 **** error checks *******
	 *************************
	 if (`deriv'<0) {
	    di as error "derivative incorrectly specified."
		exit
	 }
	 if (`dotsngrid'<0|`linengrid'<0|`cingrid'<0|`cbngrid'<0|`simsngrid'<0) {
	    di as error "number of evaluation points incorrectly specified."
		exit
	 }
	 if (`nbins'<0) {
	    di as error "number of bins incorrectly specified."
		exit
	 }		
	 if (`level'>100|`level'<0) {
	     di as error "confidence level incorrectly specified."
		 exit
	 }
	 if (`dots_p'<`dots_s') {
	     di as error "p cannot be smaller than s."
		 exit
	 }
	 if (`dots_p'<`deriv') {
	     di as error "p for dots cannot be less than deriv."
		 exit
	 }
	 if ("`line_p'"!=".") {
	    if (`line_p'<`line_s') {
		   di as error "p cannot be smaller than s."
		   exit
		}
		if (`line_p'<`deriv') {
		   di as error "p for line cannot be less than deriv."
		   exit
		}
	 }
	 if ("`ci_p'"!=".") {
	    if (`ci_p'<`ci_s') {
		   di as error "p cannot be smaller than s."
		   exit
		}
		if (`ci_p'<`deriv') {
		   di as error "p for CI cannot be less than deriv."
		   exit
		}
	 }
	 if ("`cb_p'"!=".") {
	    if (`cb_p'<`cb_s') {
		   di as error "p cannot be smaller than s."
		   exit
		}
		if (`cb_p'<`deriv') {
		   di as error "p for CB cannot be less than deriv."
		   exit
		}
	 }
	 if ("`polyreg'"!="") {
	    if (`polyreg'<`deriv') {
		   di as error "polyreg() cannot be less than deriv()."
		   exit
		}
	 }
     
	 if (`ntestshape'!=0) {
	     if (`tsha_p'<=`dots_p') {
		    di as text in gr "warning: p for testing > p for dots suggested."
	     }
	     if (`tsha_p'<`tsha_s') {
		    di as error "p cannot be smaller than s."
		    exit
	     }
	 }
	 if (`"`testmodelparfit'"'!=`""'|"`testmodelpoly'"!="") {
	     if (`"`testmodelparfit'"'!=`""') confirm file `"`testmodelparfit'.dta"'
	     if (`tmod_p'<=`dots_p') {
		    di as text in gr "warning: p for testing > p for dots suggested."
		 }
		 if (`tmod_p'<`tmod_s') {
		    di as error "p cannot be smaller than s."
			exit
		 }
	 }
	 if (`"`savedata'"'!=`""') {
	    if ("`replace'"=="") {
		   confirm new file `"`savedata'.dta"'
		}
		if ("`plot'"!="") {
		    di as error "plot cannot be turned off if graph data are requested."
			exit
		}
	 }
	 if (`polyregcingrid'!=0&"`polyreg'"=="") {
	    di as error "polyreg() is missing."
		exit
	 }
	 if ("`binsmethod'"!="DPI"&"`binsmethod'"!="ROT") {
	    di as error "binsmethod incorrectly specified."
		exit
	 }
	 ******** END error checking ***************************

	 * Mark sample
	 preserve
	 marksample touse, nov   /* do not account for missing values !! */
	 qui keep if `touse'
	 
	 * Parse varlist into y_var, x_var and w_var
	 tokenize `varlist'
	 fvrevar `1', tsonly
	 local y_var "`r(varlist)'"
	 local y_varname "`1'"
	 fvrevar `2', tsonly
	 local x_var "`r(varlist)'"
	 local x_varname "`2'"
	 
	 macro shift 2
	 local w_var "`*'"
	 fvrevar `w_var', tsonly
	 local w_var "`r(varlist)'"
	 
	 marksample touse
	 markout `touse' `by', strok	 
	 qui keep if `touse'
	 local nsize=_N             /* # of rows in the original dataset */
	 
	 if ("`masspoints'"!="off"|"`binspos'"=="QS") {
	    if ("`:sortedby'"!="`x_var'") {
		  di as text in gr "Sorting dataset on `x_varname'..."
		  di as text in gr "Note: This step is omitted if dataset already sorted by `x_varname'."
		  sort `x_var'
		}
		local sorted "sorted"
	 }
	 
	 if ("`wtype'"=="f") qui sum `x_var' `wt', meanonly
	 else 	             qui sum `x_var', meanonly

	 local xmin=r(min)
	 local xmax=r(max)
	 local Ntotal=r(N)  /* total sample size, with wt */
	 
	 * Effective sample size
	 local eN=`nsize'
	 * DO NOT check mass points and clusters outside loop unless needed
	 
	 * Check number of unique byvals & create local storing byvals
	 local byvarname `by'
     if "`by'"!="" {        
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
                
        local bynum=r(r)
        forvalues i=1/`bynum' {
           local byvals `byvals' `=`byvalmatrix'[`i',1]'
        }
     }
     else local bynum=1
	 
	 * Default colors, symbols, linepatterns
     if (`"`bycolors'"'==`""') local bycolors ///
                navy maroon forest_green dkorange teal cranberry lavender ///
                khaki sienna emidblue emerald brown erose gold bluishgray
	 if (`"`bysymbols'"'==`""') local bysymbols ///
	            O D T S + X A a | V o d s t x
	 if (`"`bylpatterns'"'==`""') {
	    forval i=1/`bynum' {
		   local bylpatterns `bylpatterns' solid 
		}
	 }        

	 * Temp name in MATA
	 tempname xvec yvec byvec cluvec binedges
	 mata: `xvec'=st_data(., "`x_var'"); `yvec'=st_data(.,"`y_var'"); `byvec'=.; `cluvec'=.
	 
	 *******************************************************
	 *** Mass point counting *******************************
	 tempname Ndistlist Nclustlist
	 mat `Ndistlist'=J(`bynum',1,.)
	 mat `Nclustlist'=J(`bynum',1,.)
	 if (`bynum'>1) mata: `byvec'=st_data(.,"`by'")
	 if ("`clusterON'"=="T") mata: `cluvec'=st_data(.,"`clustervar'")
	 
	 ********************************************************
	 ********** Bins, based on FULL sample ******************
	 ********************************************************
	 * knotlist: inner knot seq; knotlistON: local, knot available before loop
	 
	 tempname fullkmat   /* matrix name for saving knots based on the full sample */
	 
	 if ("`binsmethod'"=="DPI")  local binselectmethod "IMSE-optimal plug-in choice"
     if ("`binsmethod'"=="ROT")  local binselectmethod "IMSE-optimal rule-of-thumb choice"
	 
	 * Extract user-specified knot list
	 if ("`binspos'"!="ES"&"`binspos'"!="QS") {
		capture numlist "`binspos'", ascending
		if (_rc==0) {
		   local knotlistON "T"
		   local knotlist `binspos'
		   local nbins: word count `knotlist'
		   local first: word 1 of `knotlist'
		   local last: word `nbins' of `knotlist'
		   if (`first'<=`xmin'|`last'>=`xmax') {
			  di as error "inner knots specified out of allowed range."
			  exit
		   }
		   else {
		      local nbins=`nbins'+1
              local pos "user"
				   
			  foreach el of local knotlist {
		         mat `fullkmat'=(nullmat(`fullkmat') \ `el')
		      }
              mat `fullkmat'=(`xmin' \ `fullkmat' \ `xmax')
		   }
	    }
	    else {
		   di as error "numeric list incorrectly specified in binspos()."
		   exit
		}
	 }
	 
	 * Discrete x?
	 if ("`fewmasspoints'"!="") local fullfewobs "T"
	 
	 * Bin selection using the whole sample if
	 if ("`fullfewobs'"==""&(`nbins'==0)&(("`by'"=="")|(("`by'"!="")&("`samebinsby'"!="")))) {
	    local selectfullON "T"
	 }
	 
	 if ("`selectfullON'"=="T") {
	    local Ndist=.
	    if ("`massadj'"=="T") {
	       mata: `binedges'=binsreg_uniq(`xvec', ., 1, "Ndist")
		   mata: mata drop `binedges'
	       local eN=min(`eN', `Ndist')
	    }
	    * # of clusters
	    local Nclust=.
	    if ("`clusterON'"=="T") {
		   mata: st_local("Nclust", strofreal(rows(uniqrows(`cluvec'))))
		   local eN=min(`eN', `Nclust')   /* effective sample size */
	    }
	 
		* Check effective sample size
		if ("`nbinsrot'"==""&(`eN'<=`dfcheck_n1'+`dots_p'+1+`qrot')) {
		    di as text in gr "warning: too small effective sample size for bin selection." ///
			                 _newline _skip(9) "# of mass points or clusters used and by() option ignored."
	        local by ""                 
			local byvals ""
			local fullfewobs "T"
			local binspos "QS"                /* forced to be QS */
	    }
		else {
			qui binsregselect `y_var' `x_var' `w_var' `wt', deriv(`deriv') bins(`dots_p' `dots_s') ///
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
	 
	 if (("`selectfullON'"=="T"|(`nbins'!=0&"`samebinsby'"!=""))&"`fullfewobs'"=="") {
		 * Save in a knot list
         local knotlistON "T"
		 if ("`binspos'"=="ES") {
	        local stepsize=(`xmax'-`xmin')/`nbins'
	        forvalues i=1/`=`nbins'+1' {
		       mat `fullkmat'=(nullmat(`fullkmat') \ `=`xmin'+`stepsize'*(`i'-1)')
		    }
	     }
	     else {
			 if (`nbins'==1)  mat `kmat'=(`xmin' \ `xmax')
		     else {		
	           binsreg_pctile `x_var' `wt', nq(`nbins')
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
	 
	 * NOTE: ALL checkings are put within the loop
	 
	 * Set seed
	 set seed `simsseed'
	 
	 * alpha quantile (for two-sided CI)
	 local alpha=(100-(100-`level')/2)/100
	 
	 * Shut down test if there are many groups
	 if (`ntestshape'!=0|`"`testmodelparfit'"'!=`""'|"`testmodelpoly'"!="") { 
	    if (`bynum'!=1) {
	       local test_fewobs "T"
		   di as text in gr "warning: tests cannot be implemented with by()." 
	    }
	 }
	 
	 ***************************************************************************
	 *************** Preparation before loop************************************
	 ***************************************************************************
	 
	 ********** Prepare vars for plotting ********************
	 * names for mata objection storing graph data
	 * plotmat: final output (defined outside); 
	 * plotmatby: output for each group
	 tempname plotmat plotmatby xsub ysub byindex xcatsub
	 tempname Xm mata_fit mata_se          /* temp name for mata obj */
	 
	 * count the number of requested columns, record the positions
	 local ncolplot=1                /* 1st col reserved for group */
	 if ("`plot'"=="") {
	    if (`dotsntot'!=0) {
		   local dots_start=`ncolplot'+1
		   local dots_end=`ncolplot'+4
		   local ncolplot=`ncolplot'+4
		}
	    if (`linengrid'!=0&"`fullfewobs'"=="") {
		   local line_start=`ncolplot'+1
		   local line_end=`ncolplot'+4
		   local ncolplot=`ncolplot'+4
		}
        if (`polyregngrid'!=0) {
		   local poly_start=`ncolplot'+1
		   local poly_end=`ncolplot'+4
		   local ncolplot=`ncolplot'+4
		   if (`polyregcingrid'!=0) {
		      local polyci_start=`ncolplot'+1
		      local polyci_end=`ncolplot'+5
		      local ncolplot=`ncolplot'+5
		   }
		}
		if (`cintot'!=0) {
		   local ci_start=`ncolplot'+1
		   local ci_end=`ncolplot'+5
		   local ncolplot=`ncolplot'+5
		}
		if (`cbngrid'!=0&"`fullfewobs'"=="") {
		   local cb_start=`ncolplot'+1
		   local cb_end=`ncolplot'+5
		   local ncolplot=`ncolplot'+5
		}
	 }
	 mata: `plotmat'=J(0,`ncolplot',.)
	 
	 * mark the (varying) last row (for plotting)
	 local bylast=0               
     *******************************************************************	 
	 * temp var: bin id
	 tempvar xcat
	 qui gen `xcat'=. in 1
	 
	 * matrix names, for returns
	 tempname Nlist nbinslist cvallist
	 
	 * local vars, for plotting
	 local counter_by=1
	 local plotnum=0                     /* count the number of series, for legend */
     if ("`by'"=="") local noby="noby"
	 local byvalnamelist ""              /* save group name (value) */
	 local plotcmd ""                    /* plotting cmd */
	 
	 ***************************************************************************
	 ******************* Now, enter the loop ***********************************
	 ***************************************************************************
	 foreach byval in `byvals' `noby' {
		if ("`by'"!="") {
		    local conds "if `by'==`byval'"     /* with "if" */
	        if ("`bylabel'"=="") local byvalname=`byval'
		    else {
			   local byvalname `: label `bylabel' `byval''
			}
			local byvalnamelist `byvalnamelist' `byvalname'
		}
		if (`bynum'>1) {
		   mata: `byindex'=`byvec':==`byval'
		   mata: `xsub'=select(`xvec',`byindex'); `ysub'=select(`yvec', `byindex')
		}
		else {
		   mata: `xsub'=`xvec'; `ysub'=`yvec'
		}
		
		* Subsample size
		if ("`wtype'"=="f") sum `x_var' `conds' `wt', meanonly
		else                sum `x_var' `conds', meanonly
		
		local xmin=r(min)
	    local xmax=r(max)
	    local N=r(N)
		mat `Nlist'=(nullmat(`Nlist') \ `N')
		
	    * Effective sample size
		if (`bynum'==1) local eN=`nsize'
		else {
		   if ("`wtype'"!="f") local eN=r(N)
		   else {
		      qui count `conds'
	          local eN=r(N)
		   }
		}
	    
	    local Ndist=.
	    if ("`massadj'"=="T") {
		   mata: `binedges'=binsreg_uniq(`xsub', ., 1, "Ndist")
		   mata: mata drop `binedges'
		   local eN=min(`eN', `Ndist')
	    }
		
	    * # of clusters
	    local Nclust=.
	    if ("`clusterON'"=="T") {
		   if (`bynum'==1) {
		   	  mata: st_local("Nclust", strofreal(rows(uniqrows(`cluvec'))))
		   }
		   else {
		      mata: st_local("Nclust", strofreal(rows(uniqrows(select(`cluvec', `byvec'==`byval')))))
		   }
		   local eN=min(`eN', `Nclust')   /* effective SUBsample size */
	    }
	    
		*********************************************************
		************** Prepare bins, within loop ****************
		*********************************************************
		if ("`pos'"!="user") local pos `binspos'              /* initialize pos */
		* Selection?
	    if (`nbins_all'==0&"`knotlistON'"!="T"&"`fullfewobs'"=="") {
		   * Check effective sample size
		   if ("`nbinsrot'"==""&(`eN'<=`dfcheck_n1'+`dots_p'+1+`qrot')) {
		      di as text in gr "warning: too small effective sample size for bin selection." ///
			                   _newline _skip(9) "# of mass points or clusters used."
	          local fewobs "T"
			  local nbins=`eN'
			  local pos "QS"                /* forced to be QS */
	       }
		   else {
			  qui binsregselect `y_var' `x_var' `w_var' `conds' `wt', deriv(`deriv') bins(`dots_p' `dots_s') ///
		                        binsmethod(`binsmethod') binspos(`pos') nbinsrot(`nbinsrot') ///
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
		
		if (`nbins_all'!=0) local nbins=`nbins_all'    /* add the universal nbins */
		if ("`fullfewobs'"!="") {
		   local fewobs "T"
		   local nbins=`eN'
		}
		
	    ******************************************************
	    * Check effective sample size for each case **********
		******************************************************
	    if ("`fewobs'"!="T") {
		   if ((`nbins'-1)*(`dots_p'-`dots_s'+1)+`dots_p'+1+`dfcheck_n2'>=`eN') {
		      local fewobs "T"           /* even though ROT available, treat it as few obs case */
		      local nbins=`eN'
		      local pos "QS"
		      di as text in gr "warning: too small effective sample size for dots. # of mass points or clusters used."
	       }
		   if ("`line_p'"!=".") {
	          if ((`nbins'-1)*(`line_p'-`line_s'+1)+`line_p'+1+`dfcheck_n2'>=`eN') {
	             local line_fewobs "T"
		         di as text in gr "warning: too small effective sample size for line."
	          }
		   }
		   if ("`ci_p'"!=".") {
	          if ((`nbins'-1)*(`ci_p'-`ci_s'+1)+`ci_p'+1+`dfcheck_n2'>=`eN') {
	             local ci_fewobs "T"
		         di as text in gr "warning: too small effective sample size for CI."
	          }
		   }
		   if ("`cb_p'"!=".") {
	          if ((`nbins'-1)*(`cb_p'-`cb_s'+1)+`cb_p'+1+`dfcheck_n2'>=`eN') {
	             local cb_fewobs "T"
		         di as text in gr "warning: too small effective sample size for CB."
	          }
		   }
		   if (`ntestshape'!=0&`bynum'==1) {
	          if ((`nbins'-1)*(`tsha_p'-`tsha_s'+1)+`tsha_p'+1+`dfcheck_n2'>=`eN') {
	             local tsha_fewobs "T"
		         di as text in gr "warning: too small effective sample size for testing shape."
		         local test_fewobs "T"
	          }
		   }
		   if ((`"`testmodelparfit'"'!=`""'|"`testmodelpoly'"!="")&`bynum'==1) {
	          if ((`nbins'-1)*(`tmod_p'-`tmod_s'+1)+`tmod_p'+1+`dfcheck_n2'>=`eN') {
	             local tmod_fewobs "T"
		         di as text in gr "warning: too small effective sample size for testing models."
	             local test_fewobs "T"
	          }
		   }
	    }
		
	    if ("`polyreg'"!="") {
	       if (`polyreg'+1>=`eN') {
	          local polyreg_fewobs "T"
		      di as text in gr "warning: too small effective sample size for polynomial fit."
	       }
	    }
 
	    * Generate category variable for data and save knot in matrix
		tempname kmat
	    
		if ("`knotlistON'"=="T") {
		   mat `kmat'=`fullkmat'
		   if ("`fewobs'"=="T"&"`eN'"!="`Ndist'") {
		      if (`nbins'==1)  mat `kmat'=(`xmin' \ `xmax')
		      else {		
	            binsreg_pctile `x_var' `conds' `wt', nq(`nbins')
		        mat `kmat'=(`xmin' \ r(Q) \ `xmax')
		      }
		   }
		}
		else {
		   if ("`fewmasspoints'"==""&("`fewobs'"!="T"|"`eN'"!="`Ndist'")) {
             if ("`pos'"=="ES") {
	            local stepsize=(`xmax'-`xmin')/`nbins'
	            forvalues i=1/`=`nbins'+1' {
		           mat `kmat'=(nullmat(`kmat') \ `=`xmin'+`stepsize'*(`i'-1)')
		        }
		     }
	         else {		
		        if (`nbins'==1)  mat `kmat'=(`xmin' \ `xmax')
		        else {		
	              binsreg_pctile `x_var' `conds' `wt', nq(`nbins')
		          mat `kmat'=(`xmin' \ r(Q) \ `xmax')
		        }
		     }
		   }
		}
	
		* Renew knot list if few mass points
		if (("`fewobs'"=="T"&"`eN'"=="`Ndist'")|"`fewmasspoints'"!="") {
		   qui tab `x_var' `conds', matrow(`kmat')
		   if ("`fewmasspoints'"!="") {
		      local nbins=rowsof(`kmat')
			  local Ndist=`nbins'
			  local eN=`Ndist'
		   }
		}
		else {
		   mata: st_matrix("`kmat'", (`xmin' \ uniqrows(st_matrix("`kmat'")[|2 \ `=`nbins'+1'|])))
	       if (`nbins'!=rowsof(`kmat')-1) {
	          di as text in gr "warnings: repeated knots. Some bins dropped."
		      local nbins=rowsof(`kmat')-1
	       }
		   binsreg_irecode `x_var' `conds', knotmat(`kmat') bin(`xcat')
		   mata: `xcatsub'=st_data(., "`xcat'")
		   if (`bynum'>1) {
		      mata: `xcatsub'=select(`xcatsub', `byindex')
		   }
        }

	    *************************************************
	    **** Check for empty bins ***********************
		*************************************************
		mata: `binedges'=.               /* initialize */
		if ("`fewobs'"!="T"&"`localcheck'"=="T") {
		   mata: `binedges'=binsreg_uniq(`xsub', `xcatsub', `nbins', "uniqmin")
		   
		   if ("`dots_p'"!=".") {
		      if (`uniqmin'<`dots_p'+1) {
		         local dots_fewobs "T"
		         di as text in gr "warning: some bins have too few distinct x-values for dots."
		      }
	       }
	       if ("`line_p'"!=".") {
		      if (`uniqmin'<`line_p'+1) {
		         local line_fewobs "T"
		         di as text in gr "warning: some bins have too few distinct x-values for line."
		      }
	       }
	       if ("`ci_p'"!=".") {
		      if (`uniqmin'<`ci_p'+1) {
		         local ci_fewobs "T"
		         di as text in gr "warning: some bins have too few distinct x-values for CI."
		      }
	       }
	       if ("`cb_p'"!=".") {
		      if (`uniqmin'<`cb_p'+1) {
		         local cb_fewobs "T"
		         di as text in gr "warning: some bins have too few distinct x-values for CB."
		      }
	       }
	       if (`ntestshape'!=0&"`by'"=="") {
	          if ("`tsha_p'"!=".") {
		         if (`uniqmin'<`tsha_p'+1) {
		            local test_fewobs "T"
		            di as text in gr "warning: some bins have too few distinct x-values for testing."
	             }
		      }
	       }
	       if ((`"`testmodelparfit'"'!=`""'|"`testmodelpoly'"!="")&"`by'"=="") {
	          if ("`tmod_p'"!=".") {
		         if (`uniqmin'<`tmod_p'+1) {
		            local test_fewobs "T"
		            di as text in gr "warning: some bins have too few distinct x-values for testing."
	             }
		      }
	       }
		}
		
		* Now, save nbins in a list !!!
		mat `nbinslist'=(nullmat(`nbinslist') \ `nbins')
		
		**********************************************************
	    **** Count the number of rows needed (within loop!) ******
		**********************************************************
		local byfirst=`bylast'+1
		local byrange=0
		if ("`fewobs'"!="T") {
		   local dots_nr=`dotsngrid_mean'*`nbins'
		   if (`dotsngrid'!=0) local dots_nr=`dots_nr'+`dotsngrid'*`nbins'+`nbins'-1
		   local ci_nr=`cingrid_mean'*`nbins'
		   if (`cingrid'!=0)   local ci_nr=`ci_nr'+`cingrid'*`nbins'+`nbins'-1
		   if (`linengrid'!=0) local line_nr=`linengrid'*`nbins'+`nbins'-1
		   if (`cbngrid'!=0)   local cb_nr=`cbngrid'*`nbins'+`nbins'-1
		   if (`polyregngrid'!=0) {
		      local poly_nr=`polyregngrid'*`nbins'+`nbins'-1
		      if (`polyregcingrid'!=0) local polyci_nr=`polyregcingrid'*`nbins'+`nbins'-1
		   }	   
		   local byrange=max(`dots_nr'+0,`line_nr'+0,`ci_nr'+0,`cb_nr'+0, `poly_nr'+0, `polyci_nr'+0)
		}
		else {
		   if ("`eN'"=="`Ndist'") {
		      if (`polyregngrid'!=0) {
			     local poly_nr=`polyregngrid'*(`nbins'-1)+`nbins'-1-1
				 if (`polyregcingrid'!=0) local polyci_nr=`polyregcingrid'*(`nbins'-1)+`nbins'-1-1
			  }
		   }
		   else {
		      if (`polyregngrid'!=0) {
		         local poly_nr=`polyregngrid'*`nbins'+`nbins'-1
			     if (`polyregcingrid'!=0) local polyci_nr=`polyregcingrid'*`nbins'+`nbins'-1
			  }
		   }
		   local byrange=max(`nbins', `poly_nr'+0, `polyci_nr'+0)
		}
		local bylast=`bylast'+`byrange'
		mata: `plotmatby'=J(`byrange',`ncolplot',.)
		if ("`byval'"!="noby") {
		   mata: `plotmatby'[.,1]=J(`byrange',1,`byval')
		}
 
		************************************************
	    **** START: prepare data for plotting***********
	    ************************************************
		
		*************************************************
	    ********** dots and ci for few obs. case ********
	    *************************************************
	    if (`dotsntot'!=0&"`plot'"==""&"`fewobs'"=="T") {
		   di as text in gr "warning: dots(0 0) is used."
		   
		   local dots_first=`byfirst'
		   local dots_last=`byfirst'-1+`nbins'
	       
		   mata: `plotmatby'[|1,`dots_start'+2 \ `nbins',`dots_start'+2|]=range(1,`nbins',1)
		
		   if ("`eN'"=="`Ndist'") {
		      mata: `plotmatby'[|1,`dots_start' \ `nbins',`dots_start'|]=st_matrix("`kmat'"); ///
			        `plotmatby'[|1,`dots_start'+1 \ `nbins',`dots_start'+1|]=J(`nbins',1,1)
			  
			  * Renew knot commalist, each value forms a group
		      local xknot ""
		      forvalues i=1/`nbins' {
		         local xknot `xknot' `kmat'[`i',1]
		      }
		      local xknotcommalist : subinstr local xknot " " ",", all
		      qui replace `xcat'=1+irecode(`x_var',`xknotcommalist') `conds'
		   }
		   else {
		      tempname grid
		      mat `grid'=(`kmat'[1..`nbins',1]+`kmat'[2..`nbins'+1,1])/2
		      mata: `plotmatby'[|1,`dots_start' \ `nbins',`dots_start'|]=st_matrix("`grid'");///
			        `plotmatby'[|1,`dots_start'+1 \ `nbins',`dots_start'+1|]=J(`nbins',1,0)
		   }
		   
	       local nseries=`nbins'
	       capture reg `y_var' ibn.`xcat' `w_var' `conds' `wt', nocon `vce'
		   tempname fewobs_b fewobs_V
		   if (_rc==0) {
		      mat `fewobs_b'=e(b)
		      mat `fewobs_V'=e(V)
		      mata: binsreg_checkdrop("`fewobs_b'", "`fewobs_V'", `nseries')
		      mat `fewobs_b'=`fewobs_b'[1,1..`nseries']
		      mat `fewobs_V'=vecdiag(`fewobs_V')
		      mat `fewobs_V'=`fewobs_V'[1,1..`nseries']
		   }
		   else {
		      error _rc
		      exit _rc
		   }
		
		   mata: `plotmatby'[|1,`dots_start'+3 \ `nbins',`dots_start'+3|]=st_matrix("`fewobs_b'")'

		   local plotnum=`plotnum'+1
		   local legendnum `legendnum' `plotnum' 
		   local col: word `counter_by' of `bycolors'
		   local sym: word `counter_by' of `bysymbols'
		   local plotcmd `plotcmd' (scatter dots_fit dots_x ///
		                 in `dots_first'/`dots_last', ///
						 mcolor(`col') msymbol(`sym') `dotsplotopt')
		   
		   if (`cintot'!=0) {
		      di as text in gr "warning: ci(0 0) is used."
		   
		      mata: `plotmatby'[|1,`ci_start'+1 \ `nbins',`ci_start'+2|]=`plotmatby'[|1,`dots_start'+1 \ `nbins',`dots_start'+2|]; ///
			        `mata_se'=sqrt(st_matrix("`fewobs_V'")'); ///
			        `plotmatby'[|1,`ci_start'+3 \ `nbins',`ci_start'+3|]=`plotmatby'[|1,`dots_start'+3 \ `nbins',`dots_start'+3|]-`mata_se'*invnormal(`alpha'); ///
			        `plotmatby'[|1,`ci_start'+4 \ `nbins',`ci_start'+4|]=`plotmatby'[|1,`dots_start'+3 \ `nbins',`dots_start'+3|]+`mata_se'*invnormal(`alpha')
			  mata: mata drop `mata_se'
			  
			  local plotnum=`plotnum'+1
			  local lty: word `counter_by' of `bylpatterns'
		      local plotcmd `plotcmd' (rcap CI_l CI_r dots_x ///
			                in `dots_first'/`dots_last', ///
							sort lcolor(`col') lpattern(`lty') `ciplotopt')
		   }
	    } 
	 
	    *********************************************
	    **** The following handles the usual case ***
	    *********************************************
		* Turn on or off?
	    if (`dotsntot'!=0&"`plot'"==""&"`fewobs'"!="T"&"`dots_fewobs'"!="T") {
	       local dotsON "T"
	    }	    
		if (`linengrid'!=0&"`plot'"==""&"`line_fewobs'"!="T"&"`fewobs'"!="T") {
	       local lineON "T"
	    }
	    if (`polyregngrid'!=0&"`plot'"==""&"`polyreg_fewobs'"!="T") {
	       local polyON "T"
	    }
	    if (`cintot'!=0&"`plot'"==""&"`ci_fewobs'"!="T"&"`fewobs'"!="T") {
	       local ciON "T"
	    }
		if (`cbngrid'!=0&"`plot'"==""&"`cb_fewobs'"!="T"&"`fewobs'"!="T") {
	       local cbON "T"
	    }
		
	    ************************
	    ****** Dots ************
	    ************************
		tempname xmean

	    if ("`dotsON'"=="T") {
		   local dots_first=`byfirst'
		   local dots_last=`byfirst'+`dots_nr'-1
		   
		   * fitting
		   tempname dots_b dots_V
		   if ((`dots_p'==`ci_p'&`dots_s'==`ci_s'&"`ciON'"=="T")| ///
		       (`dots_p'==`cb_p'&`dots_s'==`cb_s'&"`cbON'"=="T")) {
			  binsreg_fit `y_var' `x_var' `w_var' `conds' `wt', deriv(`deriv') ///
	                    p(`dots_p') s(`dots_s') type(dots)  ///
			            xcat(`xcat') kmat(`kmat') dotsmean(`dotsngrid_mean') /// 
			            xname(`xsub') yname(`ysub') catname(`xcatsub') edge(`binedges') ///
						byvalue(`byval') usereg `sorted'
		   }
		   else {
		      binsreg_fit `y_var' `x_var' `w_var' `conds' `wt', deriv(`deriv') ///
	                    p(`dots_p') s(`dots_s') type(dots)  ///
			            xcat(`xcat') kmat(`kmat') dotsmean(`dotsngrid_mean') /// 
			            xname(`xsub') yname(`ysub') catname(`xcatsub') edge(`binedges') ///
						byvalue(`byval') `sorted'
		   }
		   
		   mat `dots_b'=e(bmat)
		   mat `dots_V'=e(Vmat)
		   if (`dotsngrid_mean'!=0) mat `xmean'=e(xmat)
		   
		   * prediction
		   if (`dotsngrid_mean'==0) {
		      mata: `plotmatby'[|1,`dots_start' \ `dots_nr',`dots_end'|] = ///
			                  binsreg_plotmat("`dots_b'", "`dots_V'", ., "`kmat'", ///
		                                       `nbins', `dots_p', `dots_s', `deriv', ///
							                   "dots", `dotsngrid')
		   }
		   else {
		      mata: `plotmatby'[|1,`dots_start' \ `dots_nr',`dots_end'|] = ///
			                  binsreg_plotmat("`dots_b'", "`dots_V'", ., "`kmat'", ///
		                                        `nbins', `dots_p', `dots_s', `deriv', ///
							                    "dots", `dotsngrid', "`xmean'")
		   }		  
		  
		   * dots
		   local plotnum=`plotnum'+1
		   local legendnum `legendnum' `plotnum'
		   local col: word `counter_by' of `bycolors'
		   local sym: word `counter_by' of `bysymbols'
	       local plotcmd `plotcmd' (scatter dots_fit dots_x ///
		                 in `dots_first'/`dots_last', ///
						 mcolor(`col') msymbol(`sym') `dotsplotopt')
	    }

	    **********************************************
	    ********************* Line *******************
	    **********************************************
	    if ("`lineON'"=="T") {
		   local line_first=`byfirst'
		   local line_last=`byfirst'-1+`line_nr'

		   * fitting
		   tempname line_b line_V
		   capture confirm matrix `dots_b' `dots_V'
		   if (`line_p'==`dots_p'& `line_s'==`dots_s' & _rc==0) {
		      matrix `line_b'=`dots_b'
		      matrix `line_V'=`dots_V'
		   }
		   else {
		      if ((`line_p'==`ci_p'&`line_s'==`ci_s'&"`ciON'"=="T")| ///
		          (`line_p'==`cb_p'&`line_s'==`cb_s'&"`cbON'"=="T")) {
				 binsreg_fit `y_var' `x_var' `w_var' `conds' `wt', deriv(`deriv') ///
	                      p(`line_p') s(`line_s') type(line)  ///
			              xcat(`xcat') kmat(`kmat') dotsmean(0) /// 
			              xname(`xsub') yname(`ysub') catname(`xcatsub') edge(`binedges') ///
						  byvalue(`byval') usereg `sorted'
			  }
			  else {
		         binsreg_fit `y_var' `x_var' `w_var' `conds' `wt', deriv(`deriv') ///
	                      p(`line_p') s(`line_s') type(line)  ///
			              xcat(`xcat') kmat(`kmat') dotsmean(0) /// 
			              xname(`xsub') yname(`ysub') catname(`xcatsub') edge(`binedges') /// 
						  byvalue(`byval') `sorted'
			  }
		      mat `line_b'=e(bmat)
		      mat `line_V'=e(Vmat)
		   }

		   * prediction
		   mata: `plotmatby'[|1,`line_start' \ `line_nr',`line_end'|] = ///
			              binsreg_plotmat("`line_b'", "`line_V'", ., "`kmat'", ///
		                      `nbins', `line_p', `line_s', `deriv', ///
							  "line", `linengrid')
		   
		   * line
		   local plotnum=`plotnum'+1
		   local col: word `counter_by' of `bycolors'
		   local lty: word `counter_by' of `bylpatterns'
	       local plotcmd `plotcmd' (line line_fit line_x ///
		                 in `line_first'/`line_last', sort cmissing(n) ///
						 lcolor(`col') lpattern(`lty') `lineplotopt')
	
	    }
	 
	    ***********************************
	    ******* Polynomial fit ************
	    ***********************************
	    if ("`polyON'"=="T") {
	       local poly_first=`byfirst'
		   local poly_last=`byfirst'-1+`poly_nr'

	       mata:`plotmatby'[|1,`poly_start' \ `poly_nr',`poly_start'+2|]=binsreg_grids("`kmat'",`polyregngrid')
		   
	       local poly_series ""
	       forval i=0/`polyreg' {
		      tempvar x_var_`i'
			  qui gen `x_var_`i''=`x_var'^`i' `conds'
	          local poly_series `poly_series' `x_var_`i''
		   }
		 
		   capture reg `y_var' `poly_series' `w_var' `conds' `wt', nocon
		   * store results
		   tempname poly_b poly_V
	       if (_rc==0) {
	 	      matrix `poly_b'=e(b)
			  matrix `poly_b'=`poly_b'[1, `=`deriv'+1'..`=`polyreg'+1']
			  matrix `poly_V'=e(V)
	       }
	       else {
	          error  _rc
	   	      exit _rc
           }
		 
		   * Data for derivative
		   mata: `Xm'=J(`poly_nr',0,.)   
	       forval i=`deriv'/`polyreg' {
		      mata: `Xm'=(`Xm', ///
				      `plotmatby'[|1,`poly_start' \ `poly_nr',`poly_start'|]:^(`i'-`deriv')* ///
					  factorial(`i')/factorial(`i'-`deriv'))	
	       }
		   mata:`plotmatby'[|1,`poly_start'+3 \ `poly_nr',`poly_start'+3|]=`Xm'*st_matrix("`poly_b'")'
		   
		   mata: mata drop `Xm'
		   
		   local plotnum=`plotnum'+1
		   local col: word `counter_by' of `bycolors'
		   local lty: word `counter_by' of `bylpatterns'
		   local plotcmd `plotcmd' (line poly_fit poly_x ///
		                 in `poly_first'/`poly_last', ///
						 sort lcolor(`col') lpattern(`lty') `polyregplotopt')
		
		   * add CI for global poly?
		   if (`polyregcingrid'!=0) {
	          local polyci_first=`byfirst'
			  local polyci_last=`byfirst'-1+`polyci_nr'
			  matrix `poly_V'=`poly_V'[`=`deriv'+1'..`=`polyreg'+1',`=`deriv'+1'..`=`polyreg'+1']
	          
	          mata: `plotmatby'[|1,`polyci_start' \ `polyci_nr',`polyci_start'+2|]=binsreg_grids("`kmat'", `polyregcingrid')
		   
		      mata: `Xm'=J(`polyci_nr',0,.)
	          forval i=`deriv'/`polyreg' {
		         mata:`Xm'=(`Xm', ///
				      `plotmatby'[|1,`polyci_start' \ `polyci_nr',`polyci_start'|]:^(`i'-`deriv')* ///
					  factorial(`i')/factorial(`i'-`deriv'))	
	          }
			  mata:`mata_fit'=`Xm'*st_matrix("`poly_b'")'; ///
				   `mata_se'=sqrt(rowsum((`Xm':*(st_matrix("`poly_V'")*`Xm'')'))); ///
				   `plotmatby'[|1,`polyci_start'+3 \ `polyci_nr',`polyci_start'+3|]=`mata_fit'-`mata_se'*invnormal(`alpha'); ///
				   `plotmatby'[|1,`polyci_start'+4 \ `polyci_nr',`polyci_start'+4|]=`mata_fit'+`mata_se'*invnormal(`alpha')
		      
		      mata: mata drop `Xm' `mata_fit' `mata_se'
		
		      * poly ci
			  local plotnum=`plotnum'+1
	          local plotcmd `plotcmd' (rcap polyCI_l polyCI_r polyCI_x ///
			                in `polyci_first'/`polyci_last', ///
							sort lcolor(`col') lpattern(`lty') `ciplotopt')
		   }
	    }
	  
	 
	    **********************************
	    ******* Confidence Interval ******
	    **********************************
	    if ("`ciON'"=="T") {
	       local ci_first=`byfirst'
		   local ci_last=`byfirst'-1+`ci_nr'
		   
		   * fitting
		   tempname ci_b ci_V
		   capture confirm matrix `line_b' `line_V'
		   if (`ci_p'==`line_p'& `ci_s'==`line_s' & _rc==0) {
		         matrix `ci_b'=`line_b'
		         matrix `ci_V'=`line_V'
		   }
		   else {
		      capture confirm matrix `dots_b' `dots_V'
		      if (`ci_p'==`dots_p'& `ci_s'==`dots_s' & _rc==0) {
		         matrix `ci_b'=`dots_b'
		         matrix `ci_V'=`dots_V'
		      }
		   }
		   
		   capture confirm matrix `ci_b' `ci_V' `xmean' 
		   if (_rc!=0) {
			    binsreg_fit `y_var' `x_var' `w_var' `conds' `wt', deriv(`deriv') ///
	                    p(`ci_p') s(`ci_s') type(ci)  ///
			            xcat(`xcat') kmat(`kmat') dotsmean(`cingrid_mean') /// 
			            xname(`xsub') yname(`ysub') catname(`xcatsub') edge(`binedges') ///
						byvalue(`byval') `sorted'
				 
		        mat `ci_b'=e(bmat)
		        mat `ci_V'=e(Vmat)
		        mat `xmean'=e(xmat)
		   }
		   
		   * prediction
		   if (`cingrid_mean'==0) {
		      mata: `plotmatby'[|1,`ci_start' \ `ci_nr',`ci_end'|] = ///
			            binsreg_plotmat("`ci_b'", "`ci_V'", ///
						      `=invnormal(`alpha')', "`kmat'", ///
		                      `nbins', `ci_p', `ci_s', `deriv', "ci", ///
							  `cingrid')
		   }
		   else {
		      
		      mata: `plotmatby'[|1,`ci_start' \ `ci_nr',`ci_end'|] = ///
			            binsreg_plotmat("`ci_b'", "`ci_V'", ///
						   `=invnormal(`alpha')', "`kmat'", ///
		                   `nbins', `ci_p', `ci_s', `deriv', "ci", ///
						   `cingrid', "`xmean'")
		   }
		   
		   * ci
		   local plotnum=`plotnum'+1
		   local col: word `counter_by' of `bycolors'
		   local lty: word `counter_by' of `bylpatterns'
	       local plotcmd `plotcmd' (rcap CI_l CI_r CI_x ///
		                 in `ci_first'/`ci_last', ///
						 sort lcolor(`col') lpattern(`lty') `ciplotopt')

	    }
	 
	    *******************************
	    ***** Confidence Band *********
	    *******************************
	    tempname cval
	    scalar `cval'=.
	    if ("`cbON'"=="T") {
		   * Prepare grid for plotting
	       local cb_first=`byfirst'
		   local cb_last=`byfirst'-1+`cb_nr'
		
		   * fitting
		   tempname cb_b cb_V
		   capture confirm matrix `ci_b' `ci_V'
		   if (`cb_p'==`ci_p'& `cb_s'==`ci_s' & _rc==0) {
		      matrix `cb_b'=`ci_b'
		      matrix `cb_V'=`ci_V'
		   } 
		   else {
		      capture confirm matrix `line_b' `line_V'
		      if (`cb_p'==`line_p'& `cb_s'==`line_s' & _rc==0) {
		         matrix `cb_b'=`line_b'
		         matrix `cb_V'=`line_V'
		      }
			  else {
			     capture confirm matrix `dots_b' `dots_V'
			     if (`cb_p'==`dots_p'& `cb_s'==`dots_s' & _rc==0) {
		            matrix `cb_b'=`dots_b'
		            matrix `cb_V'=`dots_V'
		         }
				 else {
				    binsreg_fit `y_var' `x_var' `w_var' `conds' `wt', deriv(`deriv') ///
	                    p(`cb_p') s(`cb_s') type(cb)  ///
			            xcat(`xcat') kmat(`kmat') dotsmean(0) /// 
			            xname(`xsub') yname(`ysub') catname(`xcatsub') edge(`binedges') ///
						byvalue(`byval') `sorted'
					mat `cb_b'=e(bmat)
		            mat `cb_V'=e(Vmat)
				 }
			  }
		   }
		
		   * Compute critical values
		   * Prepare grid for simulation
		   local uni_last=`simsngrid'*`nbins'+`nbins'-1
		   local nseries=(`cb_p'-`cb_s'+1)*(`nbins'-1)+`cb_p'+1
	       
		   tempname cb_basis
		   mata: `cb_basis'=binsreg_grids("`kmat'", `simsngrid'); ///
		         `cb_basis'=binsreg_spdes(`cb_basis'[,1], "`kmat'", `cb_basis'[,3], `cb_p', `deriv', `cb_s'); ///
		         `Xm'=binsreg_pred(`cb_basis', st_matrix("`cb_b'")[|1 \ `nseries'|]', ///
			                       st_matrix("`cb_V'")[|1,1 \ `nseries',`nseries'|], "all"); ///
			      binsreg_pval(`cb_basis', `Xm'[,2], "`cb_V'", ".", `nsims', `nseries', "two", `=`level'/100', ".", "`cval'")		  
		   mata: mata drop `cb_basis' `Xm'
	          
	       * prediction
		   mata: `plotmatby'[|1,`cb_start' \ `cb_nr',`cb_end'|] =   ///
			                binsreg_plotmat("`cb_b'", "`cb_V'",      ///
						                    `=`cval'', "`kmat'",       ///
		                                    `nbins', `cb_p', `cb_s', `deriv', ///
							                "cb", `cbngrid')	  
							  
		   * cb
		   local plotnum=`plotnum'+1
		   local col: word `counter_by' of `bycolors'
	       local plotcmd `plotcmd' (rarea CB_l CB_r CB_x ///
		                 in `cb_first'/`cb_last', sort cmissing(n) ///
		                 lcolor(`col'%0) fcolor(`col'%20) `cbplotopt')
	    }
		mat `cvallist'=(nullmat(`cvallist') \ `cval')
		
		mata: `plotmat'=(`plotmat' \ `plotmatby')
		
		*********************************
	    **** display ********************
	    *********************************	 
	    di ""
	    * Plotting
	    if ("`plot'"=="") {
		   if (`counter_by'==1) {
		      di in smcl in gr "Binscatter plot"
	          di in smcl in gr "Bin selection method: `binselectmethod'"
	          di in smcl in gr "Placement: `placement'"
	          di in smcl in gr "Derivative: `deriv'"
		      if (`"`savedata'"'!=`""') {
		         di in smcl in gr `"Output file: `savedata'.dta"'
		      }
		   }
		   di ""
		   if ("`by'"!="") {
		      di in smcl in gr "Group: `byvarname' = " in yellow "`byvalname'"
		   }
	       di in smcl in gr "{hline 30}{c TT}{hline 15}"
	       di in smcl in gr "{lalign 1:# of observations}"      _col(30) " {c |} " _col(32) as result %7.0f `N'
	       di in smcl in gr "{lalign 1:# of distinct values}"   _col(30) " {c |} " _col(32) as result %7.0f `Ndist'
	       di in smcl in gr "{lalign 1:# of clusters}"          _col(30) " {c |} " _col(32) as result %7.0f `Nclust'	 
	       di in smcl in gr "{hline 30}{c +}{hline 15}"
	       di in smcl in gr "{lalign 1:Bin selection:}"                   _col(30) " {c |} "
		   if ("`binselectmethod'"=="User-specified") {
	          di in smcl in gr "{ralign 29:Degree of polynomial}"         _col(30) " {c |} " _col(39) as result %7.0f "."
	          di in smcl in gr "{ralign 29:# of smoothness constraints}"  _col(30) " {c |} " _col(39) as result %7.0f "."
	       }
		   else {
	          di in smcl in gr "{ralign 29:Degree of polynomial}"         _col(30) " {c |} " _col(32) as result %7.0f `dots_p'
	          di in smcl in gr "{ralign 29:# of smoothness constraints}"  _col(30) " {c |} " _col(32) as result %7.0f `dots_s'
	       }
		   di in smcl in gr "{ralign 29:# of bins}"                       _col(30) " {c |} " _col(32) as result %7.0f `nbins'
	       di in smcl in gr "{hline 30}{c BT}{hline 15}"
		   di ""
		   di in smcl in gr "{hline 9}{c TT}{hline 30}"
		   di in smcl _col(10) "{c |}" in gr _col(17) "p" _col(25) "s" _col(33) "df"
		   di in smcl in gr "{hline 9}{c +}{hline 30}"
		   if (`dotsntot'!=0) {
		      local dots_df=(`dots_p'-`dots_s'+1)*(`nbins'-1)+`dots_p'+1
		      di in smcl in gr "{lalign 1: dots}"  _col(10) "{c |}" in gr _col(17) "`dots_p'" _col(25) "`dots_s'" _col(33) "`dots_df'"
		   }
		   if ("`lineON'"=="T") {
		      local line_df=(`line_p'-`line_s'+1)*(`nbins'-1)+`line_p'+1
		      di in smcl in gr "{lalign 1: line}"  _col(10) "{c |}" in gr _col(17) "`line_p'" _col(25) "`line_s'" _col(33) "`line_df'"
		   }
		   if (`cintot'!=0) {
		      local ci_df=(`ci_p'-`ci_s'+1)*(`nbins'-1)+`ci_p'+1
		      di in smcl in gr "{lalign 1: CI}"  _col(10) "{c |}" in gr _col(17) "`ci_p'" _col(25) "`ci_s'" _col(33) "`ci_df'"
		   }
		   if ("`cbON'"=="T") {
		      local cb_df=(`cb_p'-`cb_s'+1)*(`nbins'-1)+`cb_p'+1
		      di in smcl in gr "{lalign 1: CB}"  _col(10) "{c |}" in gr _col(17) "`cb_p'" _col(25) "`cb_s'" _col(33) "`cb_df'"
		   }
		   if ("`polyON'"=="T") {
		      local poly_df=`polyreg'+1
		      di in smcl in gr "{lalign 1: polyreg}"  _col(10) "{c |}" in gr _col(17) "`polyreg'" _col(25) "NA" _col(33) "`poly_df'"
		   }
		   di in smcl in gr "{hline 9}{c BT}{hline 30}"
	    }
		
		
		mata: mata drop `plotmatby'
		local ++counter_by
	 }
	 mata: mata drop `xsub' `ysub' `binedges'
	 if (`bynum'>1) mata: mata drop `byindex'
	 capture mata: mata drop `xcatsub'
	 ****************** END loop ****************************************
	 ********************************************************************
	 
	 *******************************************
	 ******* Test everything all at once *******
	 *******************************************
	 if ((`ntestshape'!=0|`"`testmodelparfit'"'!=`""'|"`testmodelpoly'"!="")&"`test_fewobs'"!="T"&"`fewobs'"!="T") {
	    local testON "T"
	 }
	 if ("`testON'"=="T") {
	     qui binsregtest `y_var' `x_var' `w_var' `wt', deriv(`deriv')  ///
		                 testmodelparfit(`testmodelparfit') testmodelpoly(`testmodelpoly') ///
				   	     testmodel(`tmod_p' `tmod_s') /// 
		                 testshapel(`testshapel') testshaper(`testshaper') testshape2(`testshape2') ///
		                 testshape(`tsha_p' `tsha_s')  ///
		                 nbins(`nbins') binspos(`binspos')  ///
		                 nsims(`nsims') simsgrid(`simsngrid') simsseed(`simsseed') `vce' ///
						 masspoints(`masspoints') dfcheck(`dfcheck_n1' `dfcheck_n2') numdist(`Ndist') numclust(`Nclust')
	     if ("`testshapel'"!="") {
		     tempname stat_shapeL pval_shapeL
		     local testvalueL `=e(testvalueL)'
	         matrix `stat_shapeL'=e(stat_shapeL)
	         matrix `pval_shapeL'=e(pval_shapeL)
	     }
	     if ("`testshaper'"!="") {
		     tempname stat_shapeR pval_shapeR
			 local testvalueR `=e(testvalueR)'
	         matrix `stat_shapeR'=e(stat_shapeR)
	         matrix `pval_shapeR'=e(pval_shapeR)  
	     }
	     if ("`testshape2'"!="") {
		     tempname stat_shape2 pval_shape2
		     local testvalue2 `=e(testvalue2)'
	 	     matrix `stat_shape2'=e(stat_shape2)
	         matrix `pval_shape2'=e(pval_shape2)
	     }
	     if ("`testmodelpoly'"!="") {
		     tempname stat_poly pval_poly
	         local testpolyp=e(testpolyp)
	         scalar `stat_poly'=e(stat_poly)
		     scalar `pval_poly'=e(pval_poly)
	     }
	     if (`"`testmodelparfit'"'!=`""') {
		     tempname stat_model pval_model
	         local testvarlist `=e(testvarlist)'
	         matrix `stat_model'=e(stat_model)
			 matrix `pval_model'=e(pval_model)
	     }
	 }
	 
	 ******************************************************
	 ****** Display test results **************************
	 if ("`testON'"=="T") {
	 
	 if (`ntestshape'!=0|`"`testmodelparfit'"'!=`""'|"`testmodelpoly'"!="") {
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
	       di in smcl in gr "{ralign 29:Degree of polynomial}"      _col(30) " {c |} " _col(32) as result %7.0f `dots_p'
	       di in smcl in gr "{ralign 29:# of smoothness constraints}"  _col(30) " {c |} " _col(32) as result %7.0f `dots_s'
	    }
		di in smcl in gr "{ralign 29:# of bins}"                 _col(30) " {c |} " _col(32) as result %7.0f `nbins'
	    di in smcl in gr "{hline 30}{c BT}{hline 15}"
     }
	 
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
	 
	 if ("`testmodelpoly'"!=""|`"`testmodelparfit'"'!=`""') {
	    di ""
	    di in smcl in gr "Model specification Tests:"
	    di in smcl in gr "Degree: `tmod_p'" _col(15) "# of smoothness constraints: `tmod_s'"
	 }
     if ("`testmodelpoly'"!="") {
	    di ""
		di in smcl in gr "{hline 19}{c TT}{hline 30}"
		di in smcl in gr "H0: mu =" _col(20) in gr ///
		                 "{c |}" _col(22) "sup |T|"  _col(40) "p value"
		di in smcl in gr "{hline 19}{c +}{hline 30}"
	    di in smcl in gr "poly. degree  " as result `testpolyp' _col(20) in gr "{c |}" ///
		                    _col(22) as result %7.3f `stat_poly' ///
		                    _col(40) as result %7.3f `pval_poly'
		di in smcl in gr "{hline 19}{c BT}{hline 30}"
	 }
	 if (`"`testmodelparfit'"'!=`""') {
	    di ""
		di in smcl in gr `"Input file: `testmodelparfit'.dta"'
	    di in smcl in gr "{hline 19}{c TT}{hline 30}"
	    di in smcl in gr "H0: mu ="  _col(20) in gr ///
		                 "{c |}" _col(22) "sup |T|"  _col(40) "p value"
		di in smcl in gr "{hline 19}{c +}{hline 30}"
		local nfitval: word count `testvarlist'
		forval i=1/`nfitval' {
	       local val: word `i' of `testvarlist'
		   local stat=`stat_model'[`i',1]
		   local pval=`pval_model'[`i',1]
	       di in smcl in yellow "{rcenter 19:`val'}" _col(20) in gr "{c |}" ///
		                    _col(22) as result %7.3f `stat' ///
		                    _col(40) as result %7.3f `pval'
	    }
		di in smcl in gr "{hline 19}{c BT}{hline 30}"
	 }
	 
	 }
	 
	 *******************************************
	 *************** Plotting ******************
	 *******************************************
	 clear
	 if ("`plotcmd'"!="") {
	    * put data back to STATA
		mata: st_local("nr", strofreal(rows(`plotmat'))) 
		qui set obs `nr'
		
		* MAKE SURE the orderings match
		qui gen group=. in 1
		if (`dotsntot'!=0) {
		   qui gen dots_x=. in 1
		   qui gen dots_isknot=. in 1
		   qui gen dots_binid=. in 1 
		   qui gen dots_fit=. in 1
		}
	    if (`linengrid'!=0&"`fullfewobs'"=="") {
		   qui gen line_x=. in 1
		   qui gen line_isknot=. in 1
		   qui gen line_binid=. in 1 
		   qui gen line_fit=. in 1
		}
        if (`polyregngrid'!=0) {
		   qui gen poly_x=. in 1
		   qui gen poly_isknot=. in 1
		   qui gen poly_binid=. in 1
		   qui gen poly_fit=. in 1
		   if (`polyregcingrid'!=0) {
		      qui gen polyCI_x=. in 1
			  qui gen polyCI_isknot=. in 1
			  qui gen polyCI_binid=. in 1
			  qui gen polyCI_l=. in 1
			  qui gen polyCI_r=. in 1
		   }
		}
		if (`cintot'!=0) {
           qui gen CI_x=. in 1
		   qui gen CI_isknot=. in 1
		   qui gen CI_binid=. in 1
		   qui gen CI_l=. in 1
		   qui gen CI_r=. in 1 
		}
		if (`cbngrid'!=0&"`fullfewobs'"=="") {
           qui gen CB_x=. in 1
		   qui gen CB_isknot=. in 1
		   qui gen CB_binid=. in 1
		   qui gen CB_l=. in 1
		   qui gen CB_r=. in 1
		}
			
		mata: st_store(.,.,`plotmat')
		
		* Legend
		local plot_legend legend(order(
	    if ("`by'"!=""&`dotsntot'!=0) {
		  forval i=1/`bynum' {
			 local byvalname: word `i' of `byvalnamelist'
	         local plot_legend `plot_legend' `: word `i' of `legendnum'' "`byvarname'=`byvalname'"		  
	      }
		  local plot_legend `plot_legend' ))
	    }
		else {
		   local plot_legend legend(off)
		}
		
		* Plot it
	    local graphcmd twoway `plotcmd', xtitle(`x_varname') ytitle(`y_varname') `plot_legend' `options' 
	    `graphcmd'
	 }
	 mata: mata drop `plotmat' `xvec' `yvec' `byvec' `cluvec'
	 
	 
	 * Save graph data ?
	 * In the normal case
     if (`"`savedata'"'!=`""'&`"`plotcmd'"'!=`""') {		
		* Add labels
		if ("`by'"!="") {
		    if ("`bystring'"=="T") {
			   label val group `bylabel'
			   decode group, gen(`byvarname')
			}
		    else {
			   qui gen `byvarname'=group
			   if ("`bylabel'"!="") label val `byvarname' `bylabel'
			}
			label var `byvarname' "Group"
			qui drop group
			order `byvarname'
		}
		else qui drop group
		
		capture confirm variable dots_x dots_binid dots_isknot dots_fit
		if (_rc==0) {
			label var dots_x "Dots: grid"
			label var dots_binid "Dots: indicator of bins"
			label var dots_isknot "Dots: indicator of inner knot"
			label var dots_fit "Dots: fitted values"
		}
		capture confirm variable line_x line_binid line_isknot line_fit
		if (_rc==0) {
			label var line_x "Line: grid"
			label var line_binid "Line: indicator of bins"
			label var line_isknot "Line: indicator of inner knot"
			label var line_fit "Line: fitted values"
		}
		capture confirm variable poly_x poly_binid poly_isknot poly_fit
		if (_rc==0) {
			label var poly_x "Poly: grid"
			label var poly_binid "Poly: indicator of bins"
			label var poly_isknot "Poly: indicator of inner knot"
			label var poly_fit "Poly: fitted values"
		}
		capture confirm variable polyCI_x polyCI_binid polyCI_isknot polyCI_l polyCI_r
		if (_rc==0) {
			label var polyCI_x "Poly confidence interval: grid"
			label var polyCI_binid "Poly confidence interval: indicator of bins"
			label var polyCI_isknot "Poly confidence interval: indicator of inner knot"
			label var polyCI_l "Poly confidence interval: left boundary"
			label var polyCI_r "Poly confidence interval: right boundary"
		}
		capture confirm variable CI_x CI_binid CI_isknot CI_l CI_r
		if (_rc==0) {
			label var CI_x "Confidence interval: grid"
			label var CI_binid "Confidence interval: indicator of bins"
			label var CI_isknot "Confidence interval: indicator of inner knot"
			label var CI_l "Confidence interval: left boundary"
			label var CI_r "Confidence interval: right boundary"
		}
		capture confirm variable CB_x CB_binid CB_isknot CB_l CB_r
		if (_rc==0) {
			label var CB_x "Confidence band: grid"			
		    label var CB_binid "Confidence band: indicator of bins"
		    label var CB_isknot "Confidence band: indicator of inner knot"
			label var CB_l "Confidence band: left boundary"
			label var CB_r "Confidence band: right boundary"
		}
	    qui save `"`savedata'"', `replace'
	 }
	 ***************************************************************************
	 	 
	 *********************************
	 ********** Return ***************
	 *********************************
	 ereturn clear
	 * # of observations
	 ereturn scalar N=`Ntotal'
	 * by group:
	 ereturn matrix N_by=`Nlist'    
	 ereturn matrix Ndist_by=`Ndistlist'
	 ereturn matrix Nclust_by=`Nclustlist'
	 ereturn matrix nbins_by=`nbinslist'
	 *ereturn matrix knot=`kmat'
	 ereturn matrix cval_by=`cvallist'
	 * Options
	 ereturn scalar level=`level'
	 ereturn scalar dots_p=`dots_p'
	 ereturn scalar dots_s=`dots_s'
	 ereturn scalar line_p=`line_p'
	 ereturn scalar line_s=`line_s'
	 ereturn scalar ci_p=`ci_p'
	 ereturn scalar ci_s=`ci_s'
	 ereturn scalar cb_p=`cb_p'
	 ereturn scalar cb_s=`cb_s'
	 ereturn scalar testshape_p=`tsha_p'
	 ereturn scalar testshape_s=`tsha_s'
	 ereturn scalar testmodel_p=`tmod_p'
	 ereturn scalar testmodel_s=`tmod_s'
	 
	 if ("`testON'"=="T") {
	    if ("`testshapel'"!="") {
		   ereturn local testvalueL `testvalueL'
	       ereturn matrix stat_shapeL=`stat_shapeL'
	       ereturn matrix pval_shapeL=`pval_shapeL'
	    }
	    if ("`testshaper'"!="") {
		   ereturn local testvalueR `testvalueR'
	       ereturn matrix stat_shapeR=`stat_shapeR'
	       ereturn matrix pval_shapeR=`pval_shapeR'  
	    }
	    if ("`testshape2'"!="") {
		   ereturn local testvalue2 `testvalue2'
	 	   ereturn matrix stat_shape2=`stat_shape2'
	       ereturn matrix pval_shape2=`pval_shape2'
	    }
	    if ("`testmodelpoly'"!="") {
	       ereturn scalar testpolyp=`testpolyp'
	       ereturn scalar stat_poly=`stat_poly'
		   ereturn scalar pval_poly=`pval_poly'
	    }
	    if (`"`testmodelparfit'"'!=`""') {
	       ereturn local testvarlist `testvarlist'
	       ereturn matrix pval_model=`pval_model'
	       ereturn matrix stat_model=`stat_model'
	    }
	 }
end

* Helper commands
* Estimation
program define binsreg_fit, eclass
     version 13
     syntax varlist(min=2 numeric ts fv) [if] [in] [fw aw pw] [, deriv(integer 0) ///
	        p(integer 0) s(integer 0) type(string) vce(passthru)  ///
			xcat(varname numeric) kmat(name) dotsmean(integer 0) ///        /* xmean: report x-mean? */
			xname(name) yname(name) catname(name) edge(name) ///
			byvalue(string) ///
			usereg sorted]                                                 /* usereg: force the command to use reg; sored: sorted data? */
	 
	 marksample touse
	 
	 if ("`weight'"!="") local wt [`weight'`exp']
	 
	 tokenize `varlist'
	 local y_var `1'
	 local x_var `2'
	 macro shift 2
	 local w_var "`*'"
	 if ("`w_var'"==""&`p'==0&("`type'"=="dots"|"`type'"=="line")&"`usereg'"=="") {
	    local ymeanON "T"
	 }
	 local nbins=rowsof(`kmat')-1
	 
	 tempname matxmean temp_b temp_V
	 mat `matxmean'=.
	 mat `temp_b'=.
	 mat `temp_V'=.
	 if (`dotsmean'!=0|"`ymeanON'"=="T") {
	    if ("`sorted'"==""|"`weight'"!="") {
	       preserve
	       if (`dotsmean'!=0&"`ymeanON'"=="T") {
	          collapse (mean) `y_var' (mean) `x_var' if `touse' `wt', by(`xcat') fast
		      mkmat `xcat' `x_var', matrix(`matxmean')
		      mkmat `y_var', matrix(`temp_b')
		      mat `temp_b'=`temp_b''               /* row vector */
		   }
		   else if (`dotsmean'!=0&"`ymeanON'"!="T") {
		      collapse (mean) `x_var' if `touse' `wt', by(`xcat') fast
		      mkmat `xcat' `x_var', matrix(`matxmean')
		   }
		   else {
		      collapse (mean) `y_var' if `touse' `wt', by(`xcat') fast
		      mkmat `y_var', matrix(`temp_b')
		      mat `temp_b'=`temp_b'' 
		   }
		   restore
		}
		else {
		   tempname output
		   if (`dotsmean'!=0&"`ymeanON'"=="T") {
		      mata: `output'=binsreg_mean((`xname',`yname'), `catname', `nbins', `edge'); ///
			        st_matrix("`temp_b'", `output'[.,3]'); ///
					st_matrix("`matxmean'", `output'[.,1..2])
		   }
		   else if (`dotsmean'!=0&"`ymeanON'"!="T") {
		      mata: `output'=binsreg_mean(`xname', `catname', `nbins', `edge'); ///
					st_matrix("`matxmean'", `output')
		   }
		   else {
		      mata: `output'=binsreg_mean(`yname', `catname', `nbins', `edge'); ///
			        st_matrix("`temp_b'", `output'[.,2]')
		   }
		   mata: mata drop `output'
		}
	 }
	 
	 * Regression?
	 if ("`ymeanON'"!="T") {
	    if (`p'==0) {
		   capture reg `y_var' ibn.`xcat' `w_var' if `touse' `wt', nocon `vce'
		   matrix `temp_b'=e(b)
		   matrix `temp_V'=e(V)
		}
		else {
	       local nseries=(`p'-`s'+1)*(`nbins'-1)+`p'+1
	       local series ""
	       forvalues i=1/`nseries' {
	          tempvar sp`i'
	          local series `series' `sp`i''
		      qui gen `sp`i''=. in 1
	       }
	
	       if ("`byvalue'"=="noby") {
		   	  mata: binsreg_st_spdes(`xname', "`series'", "`kmat'", `catname', `p', 0, `s')
		   }
		   else {
		      mata: binsreg_st_spdes(`xname', "`series'", "`kmat'", `catname', `p', 0, `s', "`touse'")
	       }
		   capture reg `y_var' `series' `w_var' if `touse' `wt', nocon `vce'
	       * store results
	       if (_rc==0) {
		      matrix `temp_b'=e(b)
		      matrix `temp_V'=e(V)
			  mata: binsreg_checkdrop("`temp_b'", "`temp_V'", `nseries')
	       }
	       else {
	          error  _rc
	   	      exit _rc
           }
	    }
	 }
	 
	 ereturn clear
	 ereturn matrix bmat=`temp_b'
	 ereturn matrix Vmat=`temp_V'
	 ereturn matrix xmat=`matxmean'         /* xcat, xbar */
end

mata:

  // Prediction for plotting
  real matrix binsreg_plotmat(string scalar eb, string scalar eV, real scalar cval, ///
                              string scalar knotname, real scalar J, ///
                              real scalar p, real scalar s, real scalar deriv, ///
                              string scalar type, real scalar ngrid, | string scalar muxmat) 
  {
    real matrix bmat, vmat, knot, xmean, eval, out, fit, se, semat, Xm, result
	real scalar nseries
	
	nseries=(p-s+1)*(J-1)+p+1
	bmat=st_matrix(eb)[|1\nseries|]'
	if (type=="ci"|type=="cb") {
	   vmat=st_matrix(eV)[|1,1\nseries,nseries|]
	}
    
    // Prepare evaluation points
    eval=J(0,3,.)
    if (args()==11) {
	   xmean=st_matrix(muxmat)
       eval=(eval \ (xmean[,2], J(J, 1, 0), xmean[,1])) 
    }
    if (ngrid!=0) {
	   eval=(eval \ binsreg_grids(knotname, ngrid))
    }
	
	fit=J(0,1,.)
	se=J(0,1,.)
	if (p==0) {
	   if (args()==11) {
	     fit=(fit \ bmat)
	   }
	   if (ngrid!=0) {
	     fit=(fit \ (bmat#(J(ngrid,1,1)\.)))
		 fit=fit[|1 \ (rows(fit)-1)|]
	   }
	   if (type=="ci"|type=="cb") {
	      semat=sqrt(diagonal(vmat))
		  if (args()==11) {
	         se=(se \ semat)
	      }
		  if (ngrid!=0) {
		     se=(se \ (semat#(J(ngrid,1,1)\.)))
		     se=se[|1 \ (rows(se)-1)|]
		  }
	   }
	   if (type=="dots"|type=="line") {
	      out=(eval, fit)
	   }
	   else {
	      out=(eval, fit-se*cval, fit+se*cval)
	   }
	}
	else {
	   Xm=binsreg_spdes(eval[,1], knotname, eval[,3], p, deriv, s)
	   if (type=="dots"|type=="line") {
	     fit=binsreg_pred(Xm, bmat, ., "xb")[,1]
		 out=(eval, fit)
	   }
	   else {
	   	 result=binsreg_pred(Xm, bmat, vmat, "all")
		 out=(eval, result[,1]-cval*result[,2], result[,1]+cval*result[,2])
	   }
	}

	if (type=="dots"|(type=="line"&(s==0|s-deriv<=0))) {
	   out[selectindex(out[,2]:==1),4]=J(sum(out[,2]),1,.)
	}
	if (type=="ci"|(type=="cb"&(s==0|s-deriv<=0))) {
	   out[selectindex(out[,2]:==1),4..5]=J(sum(out[,2]),2,.)
	}
	
	return(out)
  }

  // Mean by bin, with sorted data
  real matrix binsreg_mean(real matrix x, real vector xcat, real scalar nbins, ///
                           real vector edge, ///
					       | real vector by, real scalar byval)
  {
     real matrix subgroup, xsub, binid, binedges, out 
	 real scalar j, ncol
	 
	 xsub=x
	 if (args()>4) {
	    subgroup = by:==byval
		xsub=select(xsub, subgroup)
	 }
	 ncol = cols(xsub)
	
	 if (edge[1]==.) {
	    if (nbins==1) {
	 	   binedges = 0 \ rows(xsub)
	    } 
	    else {
	       binid = xcat
		   if (args()>4) {
		     binid = select(binid, subgroup)
		   }
		   binedges = 0\selectindex(binid[|1 \ rows(binid)-1|]-binid[|2 \ rows(binid)|])\rows(xsub)
	    }
	 }
	 else {
	    binedges = edge
	 }
	 
	 out=J(nbins,ncol,.)
	
	 for (j=2;j<=nbins+1;j++) {
	 	out[j-1,.] = mean(xsub[|binedges[j-1]+1, 1 \ binedges[j], ncol|])
	 }
	 out=(range(1,nbins,1), out)
	 
	 return(out)
   } 

  
end


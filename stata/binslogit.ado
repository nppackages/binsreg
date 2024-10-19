*! version 1.5 19-OCT-2024 

capture program drop binslogit
program define binslogit, eclass
     version 13
	 
	 syntax varlist(min=2 numeric fv ts) [if] [in] [fw pw] [, deriv(integer 0) at(string asis) nolink ///
	        logitopt(string asis) ///
	        dots(string) dotsgrid(string) dotsplotopt(string asis) ///
			line(string) linegrid(integer 20) lineplotopt(string asis) ///
			ci(string) cigrid(string) ciplotopt(string asis) /// 
			cb(string) cbgrid(integer 20) cbplotopt(string asis) ///
			polyreg(string) polyreggrid(integer 20) polyregcigrid(integer 0) polyregplotopt(string asis) ///
			by(varname) bycolors(string asis) bysymbols(string asis) bylpatterns(string asis) ///
			nbins(string) binspos(string) binsmethod(string) nbinsrot(string) ///
			pselect(numlist integer >=0) sselect(numlist integer >=0) ///
			samebinsby randcut(numlist max=1 >=0 <=1) ///
			nsims(integer 500) simsgrid(integer 20) simsseed(numlist integer max=1 >=0) ///
			dfcheck(numlist integer max=2 >=0) masspoints(string) usegtools(string) ///
			vce(passthru) level(real 95) asyvar(string) ///
			noplot savedata(string) replace ///
			plotxrange(numlist asc max=2) plotyrange(numlist asc max=2) *]
	 
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
	 * report the results for the cond. mean model?
	 if ("`link'"!="") local transform "F"
	 else              local transform "T"
	 
	 * default vce, clustered?
	 if ("`vce'"=="") local vce "vce(robust)"
	 local vcetemp: subinstr local vce "vce(" "", all
     local vcetemp: subinstr local vcetemp ")" "", all
	 tokenize "`vcetemp'"
	 if ("`1'"=="cl"|"`1'"=="clu"|"`1'"=="clus"|"`1'"=="clust"| /// 
		 "`1'"=="cluste"|"`1'"=="cluster") {
		local clusterON "T"           /* Mark cluster is specified */
		local clustervar `2'
	 }
	 if ("`vce'"=="vce(oim)"|"`vce'"=="vce(opg)") local vce_select "vce(ols)"
	 else                                         local vce_select "`vce'"
	 
	 if ("`asyvar'"=="") local asyvar "off"
	 if ("`binsmethod'"=="rot") local binsmethod "ROT"
	 if ("`binsmethod'"=="dpi") local binsmethod "DPI"
	 if ("`binsmethod'"=="")    local binsmethod "DPI"
	 if ("`binspos'"=="es") local binspos "ES"
	 if ("`binspos'"=="qs") local binspos "QS"
	 if ("`binspos'"=="")   local binspos "QS"

	 
	 * analyze options related to degrees *************
	 if ("`dots'"!="T"&"`dots'"!="F"&"`dots'"!="") {
	    numlist "`dots'", integer max(2) range(>=0)
		local dots=r(numlist)
	 }
	 if ("`line'"!="T"&"`line'"!="F"&"`line'"!="") {
	    numlist "`line'", integer max(2) range(>=0)
		local line=r(numlist)
	 }
	 if ("`ci'"!="T"&"`ci'"!="F"&"`ci'"!="") {
	    numlist "`ci'", integer max(2) range(>=0)
		local ci=r(numlist)
	 }
	 if ("`cb'"!="T"&"`cb'"!="F"&"`cb'"!="") {
	    numlist "`cb'", integer max(2) range(>=0)
		local cb=r(numlist)
	 }
	 
	 
	 if ("`dots'"=="F") {                   /* shut down dots */
	    local dots ""
	    local dotsgrid 0
	 }
	 if ("`line'"=="F") local line ""
	 if ("`ci'"=="F")   local ci ""
	 if ("`cb'"=="F")   local cb ""
	 
	 
	 ***************************************************************
	 * 4 cases: select J, select p, user specified both, and error
	 local selection ""
	 
	 * analyze nbins
	 if ("`nbins'"=="T") local nbins=0
	 local len_nbins=0
	 if ("`nbins'"!=""&"`nbins'"!="F") {
	    numlist "`nbins'", integer sort
	    local nbins=r(numlist)
	    local len_nbins: word count `nbins'
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
		if ("`nbins'"!=""|"`pselect'"!=""|"`sselect'"!="") {
		   di as error "nbins(), pselect() or sselect() incorrectly specified."
		   exit
		}
	 }

	 
	 * 1st case: select J
	 if (("`nbins'"=="0"|`len_nbins'>1|"`nbins'"=="")&("`binspos'"=="ES"|"`binspos'"=="QS")) local selection "J"
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
		if ("`dots'"!=""&"`dots'"!="T"&"`dots'"!="F") {                         /* respect user-specified dots */
		   local plist: word 1 of `dots'
		   local slist: word 2 of `dots'
		   if ("`slist'"=="") local slist `plist'
		}
		if ("`dots'"==""|"`dots'"=="T") local dots `plist' `slist'              /* selection is based on dots */
		if ("`line'"=="T") local line `plist' `slist'
		if ("`ci'"=="T") local ci `=`plist'+1' `=`slist'+1'
		if ("`cb'"=="T") local cb `=`plist'+1' `=`slist'+1'
		local len_p=1
	    local len_s=1
     }                                                                          /* e.g., binsreg y x, nbins(a b) or nbins(T) or pselect(a) nbins(T) */
	 
	 
	 * 2nd case: select P (at least for one object) 
	 if ("`selection'"!="J" & ("`dots'"==""|"`dots'"=="T"|"`line'"=="T"|"`ci'"=="T"|"`cb'"=="T")) {
	    local pselectOK "T"                                                     /* p selection CAN be turned on as long as one of the four is T */     
	 }

	 if ("`pselectOK'"=="T" & `len_nbins'==1 & (`len_p'>1|`len_s'>1)) {
	    local selection "P"
	 }                                                                         /* e.g., binsreg y x, pselect(a b) or pselect() dots(T) */
	 
	 * 3rd case: completely user-specified J and p
	 if ((`len_p'<=1&`len_s'<=1) & "`selection'"!="J") {
		local selection "NA"
		if ("`dots'"==""|"`dots'"=="T") {
		   if (`len_p'==1&`len_s'==1) local dots `plist' `slist'
	       else                       local dots `deriv' `deriv'                            /* e.g., binsreg y x or , dots(0 0) nbins(20) */
		}
		tokenize `dots'
		if ("`2'"=="") local 2 `1'
	    if ("`line'"=="T") {
		   if (`len_p'==1&`len_s'==1) local line `plist' `slist'
	       else                       local line `dots'                            
		}
	    if ("`ci'"=="T") {
		   if (`len_p'==1&`len_s'==1) local ci `=`plist'+1' `=`slist'+1'
	       else                       local ci `=`1'+1' `=`2'+1'                            
		}
	    if ("`cb'"=="T") {
		   if (`len_p'==1&`len_s'==1) local cb `=`plist'+1' `=`slist'+1'
	       else                       local cb `=`1'+1' `=`2'+1'                            
		}
	 }
	 	 
	 * exclude all other cases
	 if ("`selection'"=="") {
	    di as error "Degree, smoothness, or # of bins are not correctly specified."
        exit
	 }

	 
	 ****** Now, extract from dots, line, etc. ************ 
	 * dots
	 tokenize `dots'
	 local dots_p "`1'"
	 local dots_s "`2'"
	 if ("`dots_p'"==""|"`dots_p'"=="T") local dots_p=.
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
	 if ("`line'"=="") local linengrid=0
	 if ("`line_p'"==""|"`line_p'"=="T") local line_p=.
	 if ("`line_s'"=="") local line_s `line_p'
	 	 
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
	 if ("`ci'"=="") local cintot=0
	 if ("`ci_p'"==""|"`ci_p'"=="T") local ci_p=.
	 if ("`ci_s'"=="") local ci_s `ci_p'
	 
	 * cb
	 tokenize `cb'
	 local cb_p "`1'"
	 local cb_s "`2'"
	 local cbngrid `cbgrid'
	 if ("`cb'"=="") local cbngrid=0
	 if ("`cb_p'"==""|"`cb_p'"=="T") local cb_p=.
	 if ("`cb_s'"=="") local cb_s `cb_p'
	 
	 * Add warnings about degrees for estimation and inference
	 if ("`selection'"=="J") {
	    if ("`ci_p'"!=".") {
		   if (`ci_p'<=`dots_p') {
		      local ci_p=`dots_p'+1
			  local ci_s=`ci_p'
			  di as text "Warning: Degree for ci() has been changed. It must be greater than the degree for dots()."
		   }
		}
	    if ("`cb_p'"!=".") {
		   if (`cb_p'<=`dots_p') {
		      local cb_p=`dots_p'+1
			  local cb_s=`cb_p'
			  di as text "Warning: Degree for cb() has been changed. It must be greater than the degree for dots()."
		   }
		}
	 }
	 if ("`selection'"=="NA") {
	    if ("`ci'"!=""|"`cb'"!="") {
		   di as text "Warning: Confidence intervals/bands are valid when nbins() is much larger than IMSE-optimal choice. Compare your choice with the IMSE-optimal one obtained by binsregselect."
	    }
	 }
	 * if selection==P, compare ci_p/cb_p with P_opt later
	 
	 * poly fit
	 local polyregngrid `polyreggrid'
	 local polyregcingrid `polyregcigrid'
	 if ("`polyreg'"!="") {
	    confirm integer n `polyreg'
	 }
	 else {
	    local polyregngrid=0
	 }
	 
	 * range of x axis and y axis?
	 tokenize `plotxrange'
     local min_xr "`1'"
	 local max_xr "`2'"
	 tokenize `plotyrange'
     local min_yr "`1'"
	 local max_yr "`2'"

	 
	 * Simuls
	 local simsngrid=`simsgrid'
	 	 
	 * Record if nbins specified by users, set default
	 local nbins_full `nbins'              /* local save common nbins */
	 if ("`selection'"=="NA")   local binselectmethod "User-specified"
     else {
	    if ("`binsmethod'"=="DPI")  local binselectmethod "IMSE-optimal plug-in choice"
        if ("`binsmethod'"=="ROT")  local binselectmethod "IMSE-optimal rule-of-thumb choice"
		if ("`selection'"=="J") local binselectmethod "`binselectmethod' (select # of bins)"
		if ("`selection'"=="P") local binselectmethod "`binselectmethod' (select degree and smoothness)"
	 }
	 
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
	 
	 * evaluate at w from another dataset?
	 if (`"`at'"'!=`""'&`"`at'"'!=`"mean"'&`"`at'"'!=`"median"'&`"`at'"'!=`"0"') local atwout "user"
	 
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
	    * use gstats tab instead of tabstat/collapse
		* use gquantiles instead of _pctile
		* use gunique instead of binsreg_uniq
		* use fasterxtile instead of irecode (within binsreg_irecode)
		* shut down local checks & do not sort
	 }
	 
	 *************************
	 **** error checks *******
	 *************************
	 if (`deriv'<0) {
	    di as error "derivative incorrectly specified."
		exit
	 }
	 if (`deriv'>1&"`transform'"=="T") {
	    di as error "deriv cannot be greater than 1 if the conditional probability is requested." 
		exit
	 }
	 if (`dotsngrid'<0|`linengrid'<0|`cingrid'<0|`cbngrid'<0|`simsngrid'<0) {
	    di as error "Number of evaluation points incorrectly specified."
		exit
	 }
	 if (`level'>100|`level'<0) {
	     di as error "Confidence level incorrectly specified."
		 exit
	 }
	 if ("`dots_p'"!=".") {
		 if (`dots_p'<`dots_s') {
			 di as error "p cannot be smaller than s."
			 exit
		 }
		 if (`dots_p'<`deriv') {
			 di as error "p for dots cannot be less than deriv."
			 exit
		 }
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
     
	 if (`"`savedata'"'!=`""') {
	    if ("`replace'"=="") {
		   confirm new file `"`savedata'.dta"'
		}
		if ("`plot'"!="") {
		    di as error "Plot cannot be turned off if graph data are requested."
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

	 marksample touse
	 markout `touse' `by', strok	 
	 qui keep if `touse'
	 local nsize=_N             /* # of rows in the original dataset */
	 
	 if ("`usegtools'"==""&("`masspoints'"!="off"|"`binspos'"=="QS")) {
	    if ("`:sortedby'"!="`x_var'") {
		  di as text in gr "Sorting dataset on `x_varname'..."
		  di as text in gr "Note: This step is omitted if dataset already sorted by `x_varname'."
		  sort `x_var', stable
		}
		local sorted "sorted"
	 }
	 
	 if ("`wtype'"=="f") qui sum `x_var' `wt', meanonly
	 else 	             qui sum `x_var', meanonly

	 local xmin=r(min)
	 local xmax=r(max)
	 local Ntotal=r(N)  /* total sample size, with wt */
	 * define the support of plot
	 if ("`plotxrange'"!="") {
	    local xsc `plotxrange'
		if (wordcount("`xsc'")==1) local xsc `xsc' `xmax'
	 }
	 else local xsc `xmin' `xmax'
	 
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
		   local nbins: word count `knotlist'
		   local first: word 1 of `knotlist'
		   local last: word `nbins' of `knotlist'
		   if (`first'<=`xmin'|`last'>=`xmax') {
			  di as error "Inner knots specified out of allowed range."
			  exit
		   }
		   else {
		      local nbins=`nbins'+1
			  local nbins_full `nbins'
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
	 
	 * Discrete x?
	 if ("`fewmasspoints'"!="") local fullfewobs "T"
	 
	 * Bin selection using the whole sample if
	 if ("`fullfewobs'"==""&"`selection'"!="NA"&(("`by'"=="")|(("`by'"!="")&("`samebinsby'"!="")))) {
	    local selectfullON "T"
	 }
	 
	 if ("`selectfullON'"=="T") {
	    local Ndist=.
	    if ("`massadj'"=="T") {
	       if ("`usegtools'"=="") {
	          mata: `binedges'=binsreg_uniq(`xvec', ., 1, "Ndist")
		      mata: mata drop `binedges'
		   }
		   else {
		      qui gunique `x_var'
			  local Ndist=r(unique)
		   }
	       local eN=min(`eN', `Ndist')
	    }
	    * # of clusters
	    local Nclust=.
	    if ("`clusterON'"=="T") {
		   if ("`usegtools'"=="") {
		      mata: st_local("Nclust", strofreal(rows(uniqrows(`cluvec'))))
		   }
		   else {
		      qui gunique `clustervar'
			  local Nclust=r(unique)
		   }
		   local eN=min(`eN', `Nclust')   /* effective sample size */
	    }
	 

		* Check effective sample size
		if ("`dots_p'"==".") local dotspcheck=6
		else                 local dotspcheck=`dots_p'
		* Check effective sample size
		if ("`nbinsrot'"==""&(`eN'<=`dfcheck_n1'+`dotspcheck'+1+`qrot')) {
		    di as text in gr "Warning: Too small effective sample size for bin selection." ///
			                 _newline _skip(9) "# of mass points or clusters used and by() option ignored."
	        local by ""                 
			local byvals ""
			local fullfewobs "T"
			local binspos "QS"                /* forced to be QS */
	    }
		else {
		   local randcut1k `randcut'
		   if ("`randcut'"=="" & `Ntotal'>5000) {
	           local randcut1k=max(5000/`Ntotal', 0.01)
		       di as text in gr "Warning: To speed up computation, bin/degree selection uses a subsample of roughly max(5,000, 0.01n) observations if the sample size n>5,000. To use the full sample, set randcut(1)."
	       }
		   if ("`selection'"=="J") {
		      qui binsregselect `y_var' `x_var' `w_var' `wt', deriv(`deriv') bins(`dots_p' `dots_s') nbins(`nbins_full') ///
				   			    absorb(`absorb') reghdfeopt(`reghdfeopt') ///
								binsmethod(`binsmethod') binspos(`binspos') nbinsrot(`nbinsrot') ///
								`vce' masspoints(`masspoints') dfcheck(`dfcheck_n1' `dfcheck_n2') ///
								numdist(`Ndist') numclust(`Nclust') randcut(`randcut1k') usegtools(`sel_gtools')
			  if (e(nbinsrot_regul)==.) {
			      di as error "Bin selection fails."
				  exit
			  }
			  if ("`binsmethod'"=="ROT") {
				  local nbins=e(nbinsrot_regul)
				  mat `mat_imse_var_rot'=J(`bynum',1,e(imse_var_rot))
	              mat `mat_imse_bsq_rot'=J(`bynum',1,e(imse_bsq_rot))
			  }
			  else if ("`binsmethod'"=="DPI") {
				 local nbins=e(nbinsdpi)
				 mat `mat_imse_var_dpi'=J(`bynum',1,e(imse_var_dpi))
	             mat `mat_imse_bsq_dpi'=J(`bynum',1,e(imse_bsq_dpi))
				 if (`nbins'==.) {
				    local nbins=e(nbinsrot_regul)
 				    mat `mat_imse_var_rot'=J(`bynum',1,e(imse_var_rot))
	                mat `mat_imse_bsq_rot'=J(`bynum',1,e(imse_bsq_rot))
					di as text in gr "Warning: DPI selection fails. ROT choice used."
				 }
			  }
		   }
		   else if ("`selection'"=="P") {
		   	  qui binsregselect `y_var' `x_var' `w_var' `wt', deriv(`deriv') nbins(`nbins_full') ///
				   			    absorb(`absorb') reghdfeopt(`reghdfeopt') ///
								pselect(`plist') sselect(`slist') ///
								binsmethod(`binsmethod') binspos(`binspos') nbinsrot(`nbinsrot') ///
								`vce' masspoints(`masspoints') dfcheck(`dfcheck_n1' `dfcheck_n2') ///
								numdist(`Ndist') numclust(`Nclust') randcut(`randcut1k') usegtools(`sel_gtools')
			  if (e(prot_regul)==.) {
			      di as error "bin selection fails."
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
			  if ("`dots'"=="T"|"`dots'"=="") {
				 local dots_p=`binsp'
				 local dots_s=`binss'
			  }
			  if ("`line'"=="T") {
			     local line_p=`binsp'
				 local line_s=`binss'
			  }
			  if ("`ci'"!="T"&"`ci'"!="") {
			     if (`ci_p'<=`binsp') {
				   local ci_p=`binsp'+1
				   local ci_s=`ci_p'
				   di as text "Warning: Degree for ci() has been changed. It must be greater than the IMSE-optimal degree."
				 }
			  }
			  if ("`ci'"=="T") {
			     local ci_p=`binsp'+1
				 local ci_s=`binss'+1
			  }
			  if ("`cb'"!="T"&"`cb'"!="") {
			     if (`cb_p'<=`binsp') {
				   local cb_p=`binsp'+1
				   local cb_s=`cb_p'
				   di as text "Warning: Degree for cb() has been changed. It must be greater than the IMSE-optimal degree."
				 }
			  }
			  if ("`cb'"=="T") {
			     local cb_p=`binsp'+1
				 local cb_s=`binss'+1
			  }
		   }
		}
	 }
	 
	 if (("`selectfullON'"=="T"|("`selection'"=="NA"&"`samebinsby'"!=""))&"`fullfewobs'"=="") {
		 * Save in a knot list
         local knotlistON "T"
		 local nbins_full=`nbins'
		 if ("`binspos'"=="ES") {
	        local stepsize=(`xmax'-`xmin')/`nbins'
	        forvalues i=1/`=`nbins'+1' {
		       mat `fullkmat'=(nullmat(`fullkmat') \ `=`xmin'+`stepsize'*(`i'-1)')
		    }
	     }
	     else if ("`binspos'"=="QS") {
			 if (`nbins'==1)  mat `fullkmat'=(`xmin' \ `xmax')
		     else {		
	           binsreg_pctile `x_var' `wt', nq(`nbins') `usegtools'
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
	 if ("`simsseed'"!="") set seed `simsseed'
	 
	 * alpha quantile (for two-sided CI)
	 local alpha=(100-(100-`level')/2)/100
	 	 
	 ***************************************************************************
	 *************** Preparation before loop************************************
	 ***************************************************************************
	 
	 ********** Prepare vars for plotting ********************
	 * names for mata objects storing graph data
	 * plotmat: final output (defined outside); 
	 * plotmatby: output for each group
	 tempname plotmat plotmatby xsub ysub byindex xcatsub
	 tempname Xm Xm0 mata_fit mata_se          /* temp name for mata obj */
	 
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
	 mata: `plotmat'=J(0,`ncolplot'+3,.)
	 
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
	    local conds ""
		if ("`by'"!="") {
		    local conds "if `by'==`byval'"     /* with "if" */
	        if ("`bylabel'"=="") local byvalname=`byval'
		    else {
			   local byvalname `: label `bylabel' `byval''
			}
			local byvalnamelist `" `byvalnamelist' `"`byvalname'"' "'
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
		   if ("`usegtools'"=="") {
		      mata: `binedges'=binsreg_uniq(`xsub', ., 1, "Ndist")
		      mata: mata drop `binedges'
		   }
		   else {
		      qui gunique `x_var' `conds'
			  local Ndist=r(unique)
		   }
		   local eN=min(`eN', `Ndist')
		   mat `Ndistlist'[`counter_by',1]=`Ndist'
	    }
		
	    * # of clusters
	    local Nclust=.
	    if ("`clusterON'"=="T") {
		   if (`bynum'==1) {
		   	  if ("`usegtools'"=="") {
			     mata: st_local("Nclust", strofreal(rows(uniqrows(`cluvec'))))
		      }
			  else {
			     qui gunique `clustervar'
				 local Nclust=r(unique)
			  }
		   }
		   else {
		      if ("`usegtools'"=="") {
		         mata: st_local("Nclust", strofreal(rows(uniqrows(select(`cluvec', `byindex')))))
		      }
			  else {
			     qui gunique `clustervar' `conds'
				 local Nclust=r(unique)
			  }
		   }
		   local eN=min(`eN', `Nclust')   /* effective SUBsample size */
		   mat `Nclustlist'[`counter_by',1]=`Nclust'
	    }
	    
		*********************************************************
		************** Prepare bins, within loop ****************
		*********************************************************
		if ("`pos'"!="user") local pos `binspos'              /* initialize pos */
		* Selection?
	    if ("`selection'"!="NA"&"`knotlistON'"!="T"&"`fullfewobs'"=="") {
		   * Check effective sample size
		   if ("`dots_p'"==".") local dotspcheck=6
		   else                 local dotspcheck=`dots_p'
		   if ("`nbinsrot'"==""&(`eN'<=`dfcheck_n1'+`dotspcheck'+1+`qrot')) {
		      di as text in gr "Warning: Too small effective sample size for bin selection." ///
			                   _newline _skip(9) "# of mass points or clusters used."
	          local fewobs "T"
			  local nbins=`eN'
			  local pos "QS"                /* forced to be QS */
	       }
		   else {
		   	  local randcut1k `randcut'
		      if ("`randcut'"=="" & `N'>5000) {
			     local randcut1k=max(5000/`N', 0.01)
				 di as text in gr "Warning: To speed up computation, bin/degree selection uses a subsample of roughly max(5,000, 0.01n) observations if the sample size n>5,000. To use the full sample, set randcut(1)."
		      }
		      if ("`selection'"=="J") {
			     qui binsregselect `y_var' `x_var' `w_var' `conds' `wt', deriv(`deriv') ///
				                   bins(`dots_p' `dots_s') nbins(`nbins_full') ///
		                           absorb(`absorb') reghdfeopt(`reghdfeopt') /// 
					   			   binsmethod(`binsmethod') binspos(`pos') nbinsrot(`nbinsrot') ///
							       `vce' masspoints(`masspoints') dfcheck(`dfcheck_n1' `dfcheck_n2') ///
							       numdist(`Ndist') numclust(`Nclust') randcut(`randcut1k') usegtools(`sel_gtools')
			     if (e(nbinsrot_regul)==.) {
			        di as error "Bin selection fails."
			        exit
			     }
			     if ("`binsmethod'"=="ROT") {
		            local nbins=e(nbinsrot_regul)
				    mat `mat_imse_bsq_rot'[`counter_by',1]=e(imse_bsq_rot)
				    mat `mat_imse_var_rot'[`counter_by',1]=e(imse_var_rot)
		         }
		         else if ("`binsmethod'"=="DPI") {
			        local nbins=e(nbinsdpi)
					mat `mat_imse_bsq_dpi'[`counter_by',1]=e(imse_bsq_dpi)
				    mat `mat_imse_var_dpi'[`counter_by',1]=e(imse_var_dpi)
			        if (`nbins'==.) {
			           local nbins=e(nbinsrot_regul)
					   mat `mat_imse_bsq_rot'[`counter_by',1]=e(imse_bsq_rot)
				       mat `mat_imse_var_rot'[`counter_by',1]=e(imse_var_rot)
			           di as text in gr "Warning: DPI selection fails. ROT choice used."
			        }
		         }
			  }
		      else if ("`selection'"=="P") {
		   	     qui binsregselect `y_var' `x_var' `w_var' `wt', deriv(`deriv') nbins(`nbins_full') ///
				     			    absorb(`absorb') reghdfeopt(`reghdfeopt') ///
					  			    pselect(`plist') sselect(`slist') ///
								    binsmethod(`binsmethod') binspos(`binspos') nbinsrot(`nbinsrot') ///
								    `vce' masspoints(`masspoints') dfcheck(`dfcheck_n1' `dfcheck_n2') ///
								    numdist(`Ndist') numclust(`Nclust') randcut(`randcut1k') usegtools(`sel_gtools')
			     if (e(prot_regul)==.) {
			        di as error "Bin selection fails."
				    exit
			     }
			     if ("`binsmethod'"=="ROT") {
				    local binsp=e(prot_regul)
				    local binss=e(srot_regul)
					mat `mat_imse_bsq_rot'[`counter_by',1]=e(imse_bsq_rot)
				    mat `mat_imse_var_rot'[`counter_by',1]=e(imse_var_rot)
			     }
			     else if ("`binsmethod'"=="DPI") {
				    local binsp=e(pdpi)
				    local binss=e(sdpi)
					mat `mat_imse_bsq_dpi'[`counter_by',1]=e(imse_bsq_dpi)
				    mat `mat_imse_var_dpi'[`counter_by',1]=e(imse_var_dpi)
				    if (`binsp'==.) {
				       local binsp=e(prot_regul)
					   local binss=e(srot_regul)
					   mat `mat_imse_bsq_rot'[`counter_by',1]=e(imse_bsq_rot)
				       mat `mat_imse_var_rot'[`counter_by',1]=e(imse_var_rot)
					   di as text in gr "Warning: DPI selection fails. ROT choice used."
				    }
			     }
				 if ("`dots'"=="T"|"`dots'"=="") {
				    local dots_p=`binsp'
					local dots_s=`binss'
				 }
				 if ("`line'"=="T") {
				    local line_p=`binsp'
					local line_s=`binss'
				 }
				 if ("`ci'"!="T"&"`ci'"!="") {
					if (`ci_p'<=`binsp') {
					   local ci_p=`binsp'+1
					   local ci_s=`ci_p'
				       di as text "Warning: Degree for ci() has been changed. It must be greater than the IMSE-optimal degree."
				    }
			     }
				 if ("`ci'"=="T") {
					local ci_p=`binsp'+1
					local ci_s=`binss'+1
				 }
				 if ("`cb'"!="T"&"`cb'"!="") {
			        if (`cb_p'<=`binsp') {
				      local cb_p=`binsp'+1
					  local cb_s=`cb_p'
				      di as text "Warning: Degree for cb() has been changed. It must be greater than the IMSE-optimal degree."
				    }
			     }
				 if ("`cb'"=="T") {
				    local cb_p=`binsp'+1
					local cb_s=`binss'+1
				 }
		      }		   
		   }
	    }
		
		if ("`selection'"=="NA"|"`knotlistON'"=="T") local nbins=`nbins_full'    /* add the universal nbins */
		*if ("`knotlistON'"=="T") local nbins=`nbins_full'
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
		      di as text in gr "Warning: Too small effective sample size for dots. # of mass points or clusters used."
	       }
		   if ("`line_p'"!=".") {
	          if ((`nbins'-1)*(`line_p'-`line_s'+1)+`line_p'+1+`dfcheck_n2'>=`eN') {
	             local line_fewobs "T"
		         di as text in gr "Warning: Too small effective sample size for line."
	          }
		   }
		   if ("`ci_p'"!=".") {
	          if ((`nbins'-1)*(`ci_p'-`ci_s'+1)+`ci_p'+1+`dfcheck_n2'>=`eN') {
	             local ci_fewobs "T"
		         di as text in gr "Warning: Too small effective sample size for CI."
	          }
		   }
		   if ("`cb_p'"!=".") {
	          if ((`nbins'-1)*(`cb_p'-`cb_s'+1)+`cb_p'+1+`dfcheck_n2'>=`eN') {
	             local cb_fewobs "T"
		         di as text in gr "Warning: Too small effective sample size for CB."
	          }
		   }
	    }
		
	    if ("`polyreg'"!="") {
	       if (`polyreg'+1>=`eN') {
	          local polyreg_fewobs "T"
		      di as text in gr "Warning: Too small effective sample size for polynomial fit."
	       }
	    }
 
	    * Generate category variable for data and save knot in matrix
		tempname kmat
	    
		if ("`knotlistON'"=="T") {
		   mat `kmat'=`fullkmat'
		   if ("`fewobs'"=="T"&"`eN'"!="`Ndist'") {
		      if (`nbins'==1)  mat `kmat'=(`xmin' \ `xmax')
		      else {		
	            binsreg_pctile `x_var' `conds' `wt', nq(`nbins') `usegtools'
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
	              binsreg_pctile `x_var' `conds' `wt', nq(`nbins') `usegtools'
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
	          di as text in gr "Warning: Repeated knots. Some bins dropped."
		      local nbins=rowsof(`kmat')-1
	       }
		   
		   binsreg_irecode `x_var' `conds', knotmat(`kmat') bin(`xcat') ///
		                                   `usegtools' nbins(`nbins') pos(`pos') knotliston(`knotlistON')
		   
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
		   mata: st_local("Ncat", strofreal(rows(uniqrows(`xcatsub'))))
		   if (`nbins'==`Ncat') {
		      mata: `binedges'=binsreg_uniq(`xsub', `xcatsub', `nbins', "uniqmin")
		   }
		   else {
		      local uniqmin=0
			  di as text in gr "Warning: There are empty bins. Specify a smaller number in nbins()."
		   }
		   
		   if ("`dots_p'"!=".") {
		      if (`uniqmin'<`dots_p'+1) {
		         local dots_fewobs "T"
		         di as text in gr "Warning: Some bins have too few distinct x-values for dots."
		      }
	       }
	       if ("`line_p'"!=".") {
		      if (`uniqmin'<`line_p'+1) {
		         local line_fewobs "T"
		         di as text in gr "Warning: Some bins have too few distinct x-values for line."
		      }
	       }
	       if ("`ci_p'"!=".") {
		      if (`uniqmin'<`ci_p'+1) {
		         local ci_fewobs "T"
		         di as text in gr "Warning: Some bins have too few distinct x-values for CI."
		      }
	       }
	       if ("`cb_p'"!=".") {
		      if (`uniqmin'<`cb_p'+1) {
		         local cb_fewobs "T"
		         di as text in gr "Warning: Some bins have too few distinct x-values for CB."
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
		mata: `plotmatby'=J(`byrange',`ncolplot'+3,.)
		if ("`byval'"!="noby") {
		   mata: `plotmatby'[.,1]=J(`byrange',1,`byval')
		}
 
		************************************************
	    **** START: prepare data for plotting***********
	    ************************************************
		local plotcmdby ""
		
		********************************
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

		
		*************************************************
	    ********** dots and ci for few obs. case ********
	    *************************************************
	    if (`dotsntot'!=0&"`plot'"==""&"`fewobs'"=="T") {
		   di as text in gr "Warning: dots(0 0) is used."
		   if (`deriv'>0) di as text in gr "Warning: deriv(0 0) is used."
		   
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
		      mata: `plotmatby'[|1,`dots_start' \ `nbins',`dots_start'|]=st_matrix("`grid'"); ///
			        `plotmatby'[|1,`dots_start'+1 \ `nbins',`dots_start'+1|]=J(`nbins',1,0)
		   }
		   
	       local nseries=`nbins'
	       capture logit `y_var' ibn.`xcat' `w_var' `conds' `wt', nocon `vce' `logitopt'
		   tempname fewobs_b fewobs_V
		   if (_rc==0) {
		      mat `fewobs_b'=e(b)
		      mat `fewobs_V'=e(V)
		      mata: binsreg_checkdrop("`fewobs_b'", "`fewobs_V'", `nseries')
			  if (`nwvar'>0) {
			     mat `fewobs_b'=`fewobs_b'[1,1..`nseries']+(`fewobs_b'[1,`=`nseries'+1'..`=`nseries'+`nwvar'']*`wval'')*J(1,`nseries',1)
		      }
			  else {
			     mat `fewobs_b'=`fewobs_b'[1,1..`nseries']
		      }
		   }
		   else {
		      error _rc
		      exit _rc
		   }
		
		   if ("`transform'"=="T") {
		      mata: `plotmatby'[|1,`dots_start'+3 \ `nbins',`dots_start'+3|]=logistic(st_matrix("`fewobs_b'"))'
           }
		   else {
		      mata: `plotmatby'[|1,`dots_start'+3 \ `nbins',`dots_start'+3|]=st_matrix("`fewobs_b'")'
		   }
		   
		   local plotnum=`plotnum'+1
		   local legendnum `legendnum' `plotnum' 
		   local col: word `counter_by' of `bycolors'
		   local sym: word `counter_by' of `bysymbols'
		   local plotcond ""
		   if ("`plotxrange'"!=""|"`plotyrange'"!="") {
		      local plotcond `plotcond' if 
			  if ("`plotxrange'"!="") {
		         local plotcond `plotcond' dots_x>=`min_xr'
			     if ("`max_xr'"!="") local plotcond `plotcond' &dots_x<=`max_xr' 
		      }
		      if ("`plotyrange'"!="") {
		         if ("`plotxrange'"=="") local plotcond `plotcond' dots_fit>=`min_yr'
				 else                    local plotcond `plotcond' &dots_fit>=`min_yr'
			     if ("`max_yr'"!="") local plotcond `plotcond' &dots_fit<=`max_yr' 
		      }
		   }

		   local plotcmdby `plotcmdby' (scatter dots_fit dots_x ///
		                   `plotcond' in `dots_first'/`dots_last', ///
						   mcolor(`col') msymbol(`sym') `dotsplotopt')
		   
		   if (`cintot'!=0) {
		      di as text in gr "Warning: ci(0 0) is used."
		   
		   	  if (`nwvar'>0) {
			     mata: `mata_se'=(I(`nseries'), J(`nseries',1,1)#st_matrix("`wval'"))
			  }
			  else {
			     mata: `mata_se'=I(`nseries')
			  }

		      mata: `plotmatby'[|1,`ci_start'+1 \ `nbins',`ci_start'+2|]=`plotmatby'[|1,`dots_start'+1 \ `nbins',`dots_start'+2|]; ///
			        `mata_se'=sqrt(rowsum((`mata_se'*st_matrix("`fewobs_V'")):*`mata_se'))
			  if ("`transform'"=="T") {
			     mata: `mata_se'=`mata_se':*(logisticden(st_matrix("`fewobs_b'"))')
			  }
			  mata: `plotmatby'[|1,`ci_start'+3 \ `nbins',`ci_start'+3|]=`plotmatby'[|1,`dots_start'+3 \ `nbins',`dots_start'+3|]-`mata_se'*invnormal(`alpha'); ///
			        `plotmatby'[|1,`ci_start'+4 \ `nbins',`ci_start'+4|]=`plotmatby'[|1,`dots_start'+3 \ `nbins',`dots_start'+3|]+`mata_se'*invnormal(`alpha')
			  mata: mata drop `mata_se'
			  
			  local plotnum=`plotnum'+1
			  local lty: word `counter_by' of `bylpatterns'
		      local plotcmdby `plotcmdby' (rcap CI_l CI_r dots_x ///
			                  `plotcond' in `dots_first'/`dots_last', ///
							  sort lcolor(`col') lpattern(`lty') `ciplotopt')
		   }
	    } 
	 
	    *********************************************
	    **** The following handles the usual case ***
	    *********************************************
		* Turn on or off?
		local dotsON ""
		local lineON ""
		local polyON ""
		local ciON   ""
		local cbON   ""
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
		   if (("`dots_p'"=="`ci_p'"&"`dots_s'"=="`ci_s'"&"`ciON'"=="T")| ///
		       ("`dots_p'"=="`cb_p'"&"`dots_s'"=="`cb_s'"&"`cbON'"=="T")) {
			  binslogit_fit `y_var' `x_var' `w_var' `conds' `wt', deriv(`deriv') ///
	                    p(`dots_p') s(`dots_s') type(dots) `vce' ///
			            xcat(`xcat') kmat(`kmat') dotsmean(`dotsngrid_mean') /// 
			            xname(`xsub') yname(`ysub') catname(`xcatsub') edge(`binedges') ///
						usereg `sorted' `usegtools' logitopt(`logitopt')
		   }
		   else {
		      binslogit_fit `y_var' `x_var' `w_var' `conds' `wt', deriv(`deriv') ///
	                    p(`dots_p') s(`dots_s') type(dots) `vce' ///
			            xcat(`xcat') kmat(`kmat') dotsmean(`dotsngrid_mean') /// 
			            xname(`xsub') yname(`ysub') catname(`xcatsub') edge(`binedges') ///
						`sorted' `usegtools' logitopt(`logitopt')
		   }
		   
		   mat `dots_b'=e(bmat)
		   mat `dots_V'=e(Vmat)
		   if (`dotsngrid_mean'!=0) mat `xmean'=e(xmat)
		   
		   * prediction
		   if (`dotsngrid_mean'==0) {
		      mata: `plotmatby'[|1,`dots_start' \ `dots_nr',`dots_end'|] = ///
			                  binslogit_plotmat("`dots_b'", "`dots_V'", ., "`kmat'", ///
		                                       `nbins', `dots_p', `dots_s', `deriv', ///
							                   "dots", `dotsngrid', "`wval'", `nwvar', ///
											   "`transform'", "`asyvar'")
		   }
		   else {
		      mata: `plotmatby'[|1,`dots_start' \ `dots_nr',`dots_end'|] = ///
			                  binslogit_plotmat("`dots_b'", "`dots_V'", ., "`kmat'", ///
		                                        `nbins', `dots_p', `dots_s', `deriv', ///
							                    "dots", `dotsngrid', "`wval'", `nwvar', ///
												"`transform'", "`asyvar'", "`xmean'")
		   }		  
		  
		   * dots
		   local plotnum=`plotnum'+1
		   if ("`cbON'"=="T") local legendnum `legendnum' `=`plotnum'+1'
		   else {
		      local legendnum `legendnum' `plotnum'
		   }
		   local col: word `counter_by' of `bycolors'
		   local sym: word `counter_by' of `bysymbols'
		   local plotcond ""
		   if ("`plotxrange'"!=""|"`plotyrange'"!="") {
		      local plotcond if
		      if ("`plotxrange'"!="") {
		         local plotcond `plotcond' dots_x>=`min_xr'
			     if ("`max_xr'"!="") local plotcond `plotcond' &dots_x<=`max_xr' 
		      }
			  if ("`plotyrange'"!="") {
		         if ("`plotxrange'"=="") local plotcond `plotcond' dots_fit>=`min_yr'
				 else                    local plotcond `plotcond' &dots_fit>=`min_yr'
			     if ("`max_yr'"!="") local plotcond `plotcond' &dots_fit<=`max_yr' 
		      }
		   }
		   
	       local plotcmdby `plotcmdby' (scatter dots_fit dots_x ///
		                   `plotcond' in `dots_first'/`dots_last', ///
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
		   if ("`line_p'"=="`dots_p'"& "`line_s'"=="`dots_s'" & _rc==0) {
		      matrix `line_b'=`dots_b'
		      matrix `line_V'=`dots_V'
		   }
		   else {
		      if (("`line_p'"=="`ci_p'"&"`line_s'"=="`ci_s'"&"`ciON'"=="T")| ///
		          ("`line_p'"=="`cb_p'"&"`line_s'"=="`cb_s'"&"`cbON'"=="T")) {
				 binslogit_fit `y_var' `x_var' `w_var' `conds' `wt', deriv(`deriv') ///
	                      p(`line_p') s(`line_s') type(line) `vce' ///
			              xcat(`xcat') kmat(`kmat') dotsmean(0) /// 
			              xname(`xsub') yname(`ysub') catname(`xcatsub') edge(`binedges') ///
						  usereg `sorted' `usegtools' logitopt(`logitopt')
			  }
			  else {
		         binslogit_fit `y_var' `x_var' `w_var' `conds' `wt', deriv(`deriv') ///
	                      p(`line_p') s(`line_s') type(line) `vce' ///
			              xcat(`xcat') kmat(`kmat') dotsmean(0) /// 
			              xname(`xsub') yname(`ysub') catname(`xcatsub') edge(`binedges') /// 
						  `sorted' `usegtools' logitopt(`logitopt')
			  }
		      mat `line_b'=e(bmat)
		      mat `line_V'=e(Vmat)
		   }

		   * prediction
		   mata: `plotmatby'[|1,`line_start' \ `line_nr',`line_end'|] = ///
			              binslogit_plotmat("`line_b'", "`line_V'", ., "`kmat'", ///
		                      `nbins', `line_p', `line_s', `deriv', ///
							  "line", `linengrid', "`wval'", `nwvar', "`transform'", "`asyvar'")
		   
		   * line
		   local plotnum=`plotnum'+1
		   local col: word `counter_by' of `bycolors'
		   local lty: word `counter_by' of `bylpatterns'
		   local plotcond ""
		   if ("`plotxrange'"!=""|"`plotyrange'"!="") {
		      local plotcond if
		      if ("`plotxrange'"!="") {
		         local plotcond `plotcond' line_x>=`min_xr'
			     if ("`max_xr'"!="") local plotcond `plotcond' &line_x<=`max_xr' 
		      }
			  if ("`plotyrange'"!="") {
		         if ("`plotxrange'"=="") local plotcond `plotcond' line_fit>=`min_yr'
				 else                    local plotcond `plotcond' &line_fit>=`min_yr'
			     if ("`max_yr'"!="") local plotcond `plotcond' &(line_fit<=`max_yr'|line_fit==.) 
		      }
		   }
		   
	       local plotcmdby `plotcmdby' (line line_fit line_x ///
		                   `plotcond' in `line_first'/`line_last', sort cmissing(n) ///
						   lcolor(`col') lpattern(`lty') `lineplotopt')
	
	    }
	 
	    ***********************************
	    ******* Polynomial fit ************
	    ***********************************
	    if ("`polyON'"=="T") {
		   if (`nwvar'>0) {
		      di as text "Note: When additional covariates w are included, the polynomial fit may not always be close to the binscatter fit."
		   }

	       local poly_first=`byfirst'
		   local poly_last=`byfirst'-1+`poly_nr'

	       mata:`plotmatby'[|1,`poly_start' \ `poly_nr',`poly_start'+2|]=binsreg_grids("`kmat'",`polyregngrid')
		   
	       local poly_series ""
	       forval i=0/`polyreg' {
		      tempvar x_var_`i'
			  qui gen `x_var_`i''=`x_var'^`i' `conds'
	          local poly_series `poly_series' `x_var_`i''
		   }
		 
		   capture logit `y_var' `poly_series' `w_var' `conds' `wt', nocon `vce' `logitopt'
		   * store results
		   tempname poly_b poly_V poly_adjw
	       if (_rc==0) {
	 	      matrix `poly_b'=e(b)
			  matrix `poly_V'=e(V)
	       }
	       else {
	          error  _rc
	   	      exit _rc
           }
		 
		   * Data for derivative
		   mata: `Xm'=J(`poly_nr',0,.); `Xm0'=J(`poly_nr',0,.)   
	       forval i=`deriv'/`polyreg' {
		      mata: `Xm'=(`Xm', ///
				      `plotmatby'[|1,`poly_start' \ `poly_nr',`poly_start'|]:^(`i'-`deriv')* ///
					  factorial(`i')/factorial(`i'-`deriv'))	
	       }
		   mata: `Xm'=(J(`poly_nr', `deriv',0), `Xm')
		   if (`nwvar'>0) {
		      if (`deriv'==0) mata: `Xm'=(`Xm', J(`poly_nr',1,1)#st_matrix("`wval'"))
			  else            mata: `Xm'=(`Xm', J(`poly_nr',`nwvar',0))
		   }
		   
		   if ("`transform'"=="T") {
		      if (`deriv'==0) {
			     mata:`plotmatby'[|1,`poly_start'+3 \ `poly_nr',`poly_start'+3|]=logistic(`Xm'*st_matrix("`poly_b'")')
			  }
			  else if (`deriv'==1) {
			     forval i=0/`polyreg' {
		            mata: `Xm0'=(`Xm0', `plotmatby'[|1,`poly_start' \ `poly_nr',`poly_start'|]:^`i')	
	             }
				 if (`nwvar'>0) mata: `Xm0'=(`Xm0', J(`poly_nr',1,1)#st_matrix("`wval'"))
			     mata:`plotmatby'[|1,`poly_start'+3 \ `poly_nr',`poly_start'+3|]=logisticden(`Xm0'*st_matrix("`poly_b'")'):* ///
				                                                                 (`Xm'*st_matrix("`poly_b'")')
			  }
		   }
		   else {
		      mata:`plotmatby'[|1,`poly_start'+3 \ `poly_nr',`poly_start'+3|]=`Xm'*st_matrix("`poly_b'")'
		   }
		   
		   mata: mata drop `Xm' `Xm0'
		   
		   local plotnum=`plotnum'+1
		   local col: word `counter_by' of `bycolors'
		   local lty: word `counter_by' of `bylpatterns'
		   local plotcond ""
		   if ("`plotxrange'"!=""|"`plotyrange'"!="") {
		      local plotcond if
		      if ("`plotxrange'"!="") {
		         local plotcond `plotcond' poly_x>=`min_xr'
			     if ("`max_xr'"!="") local plotcond `plotcond' &poly_x<=`max_xr' 
		      }
			  if ("`plotyrange'"!="") {
		         if ("`plotxrange'"=="") local plotcond `plotcond' poly_fit>=`min_yr'
				 else                    local plotcond `plotcond' &poly_fit>=`min_yr'
			     if ("`max_yr'"!="") local plotcond `plotcond' &poly_fit<=`max_yr' 
		      }
		   }

		   local plotcmdby `plotcmdby' (line poly_fit poly_x ///
		                   `plotcond' in `poly_first'/`poly_last', ///
						   sort lcolor(`col') lpattern(`lty') `polyregplotopt')
		
		   * add CI for global poly?
		   if (`polyregcingrid'!=0) {
	          local polyci_first=`byfirst'
			  local polyci_last=`byfirst'-1+`polyci_nr'
	          
	          mata: `plotmatby'[|1,`polyci_start' \ `polyci_nr',`polyci_start'+2|]=binsreg_grids("`kmat'", `polyregcingrid')
		   
		      mata: `Xm'=J(`polyci_nr',0,.); `Xm0'=J(`polyci_nr',0,.)
	          forval i=`deriv'/`polyreg' {
		         mata:`Xm'=(`Xm', ///
				      `plotmatby'[|1,`polyci_start' \ `polyci_nr',`polyci_start'|]:^(`i'-`deriv')* ///
					  factorial(`i')/factorial(`i'-`deriv'))	
	          }
 		      mata: `Xm'=(J(`polyci_nr', `deriv',0), `Xm')
		      if (`nwvar'>0) {
		         if (`deriv'==0) mata: `Xm'=(`Xm', J(`polyci_nr',1,1)#st_matrix("`wval'"))
			     else            mata: `Xm'=(`Xm', J(`polyci_nr',`nwvar',0))
		      }
			  
			  if ("`transform'"=="T") {
		         if (`deriv'==0) {
				    mata: `mata_fit'=logistic(`Xm'*st_matrix("`poly_b'")')
					mata: `mata_se'=logisticden(`Xm'*st_matrix("`poly_b'")'):* ///
					                sqrt(rowsum((`Xm'*st_matrix("`poly_V'")):*`Xm'))
				 }
				 else if (`deriv'==1) {
				    forval i=0/`polyreg' {
		               mata: `Xm0'=(`Xm0', `plotmatby'[|1,`polyci_start' \ `polyci_nr',`polyci_start'|]:^`i')	
	                }
					if (`nwvar'>0) mata: `Xm0'=(`Xm0', J(`polyci_nr',1,1)#st_matrix("`wval'"))
			        mata:`mata_fit'=logisticden(`Xm0'*st_matrix("`poly_b'")'):* ///
				                    (`Xm'*st_matrix("`poly_b'")')
									
					tempname tempobj
					mata: `tempobj'=`Xm0'*st_matrix("`poly_b'")'; ///
					      `tempobj'=logisticden(`tempobj'):*(1:-2*logistic(`tempobj')):*(`Xm'*st_matrix("`poly_b'")'):*`Xm0' + ///
					                logisticden(`tempobj'):*`Xm'; ///
						  `mata_se'=sqrt(rowsum((`tempobj'*st_matrix("`poly_V'")):*`tempobj'))
					mata: mata drop `tempobj'
				 }
			  }
			  else {
		         mata: `mata_fit'=`Xm'*st_matrix("`poly_b'")';	///			 
			           `mata_se'=sqrt(rowsum((`Xm'*st_matrix("`poly_V'")):*`Xm')) 
			  }

			  mata:`plotmatby'[|1,`polyci_start'+3 \ `polyci_nr',`polyci_start'+3|]=`mata_fit'-`mata_se'*invnormal(`alpha'); ///
				   `plotmatby'[|1,`polyci_start'+4 \ `polyci_nr',`polyci_start'+4|]=`mata_fit'+`mata_se'*invnormal(`alpha'); ///
				   `plotmatby'[selectindex(`plotmatby'[,`=`polyci_start'+1']:==1),(`=`polyci_start'+3',`=`polyci_start'+4')]=J(`=`nbins'-1',2,.)
		      
		      mata: mata drop `Xm' `Xm0' `mata_fit' `mata_se'
		
		      * poly ci
			  local plotnum=`plotnum'+1
			  local plotcond ""
		      if ("`plotxrange'"!=""|"`plotyrange'"!="") {
		         local plotcond if
		         if ("`plotxrange'"!="") {
		            local plotcond `plotcond' polyCI_x>=`min_xr'
			        if ("`max_xr'"!="") local plotcond `plotcond' &polyCI_x<=`max_xr' 
		         }
			     if ("`plotyrange'"!="") {
		            if ("`plotxrange'"=="") local plotcond `plotcond' polyCI_l>=`min_yr'
				    else                    local plotcond `plotcond' &polyCI_l>=`min_yr'
			        if ("`max_yr'"!="") local plotcond `plotcond' &polyCI_r<=`max_yr' 
		         }
		      }

	          local plotcmdby `plotcmdby' (rcap polyCI_l polyCI_r polyCI_x ///
			                  `plotcond' in `polyci_first'/`polyci_last', ///
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
		   if ("`ci_p'"=="`line_p'"& "`ci_s'"=="`line_s'" & _rc==0) {
		         matrix `ci_b'=`line_b'
		         matrix `ci_V'=`line_V'
		   }
		   else {
		      capture confirm matrix `dots_b' `dots_V'
		      if ("`ci_p'"=="`dots_p'"& "`ci_s'"=="`dots_s'" & _rc==0) {
		         matrix `ci_b'=`dots_b'
		         matrix `ci_V'=`dots_V'
		      }
		   }
		   
		   capture confirm matrix `ci_b' `ci_V' `xmean' 
		   if (_rc!=0) {
			    binslogit_fit `y_var' `x_var' `w_var' `conds' `wt', deriv(`deriv') ///
	                    p(`ci_p') s(`ci_s') type(ci) `vce' ///
			            xcat(`xcat') kmat(`kmat') dotsmean(`cingrid_mean') /// 
			            xname(`xsub') yname(`ysub') catname(`xcatsub') edge(`binedges') ///
						`sorted' `usegtools' logitopt(`logitopt')
				 
		        mat `ci_b'=e(bmat)
		        mat `ci_V'=e(Vmat)
		        mat `xmean'=e(xmat)
		   }
		   
		   * prediction
		   if (`cingrid_mean'==0) {
		      mata: `plotmatby'[|1,`ci_start' \ `ci_nr',`ci_end'|] = ///
			            binslogit_plotmat("`ci_b'", "`ci_V'", ///
						      `=invnormal(`alpha')', "`kmat'", ///
		                      `nbins', `ci_p', `ci_s', `deriv', "ci", ///
							  `cingrid', "`wval'", `nwvar', "`transform'", "`asyvar'")
		   }
		   else {
		      mata: `plotmatby'[|1,`ci_start' \ `ci_nr',`ci_end'|] = ///
			            binslogit_plotmat("`ci_b'", "`ci_V'", ///
						   `=invnormal(`alpha')', "`kmat'", ///
		                   `nbins', `ci_p', `ci_s', `deriv', "ci", ///
						   `cingrid', "`wval'", `nwvar', ///
						   "`transform'", "`asyvar'", "`xmean'")
		   }
		   
		   * ci
		   local plotnum=`plotnum'+1
		   local col: word `counter_by' of `bycolors'
		   local lty: word `counter_by' of `bylpatterns'
	       local plotcond ""
		   if ("`plotxrange'"!=""|"`plotyrange'"!="") {
		      local plotcond if
		      if ("`plotxrange'"!="") {
	             local plotcond `plotcond' CI_x>=`min_xr'
			     if ("`max_xr'"!="") local plotcond `plotcond' &CI_x<=`max_xr' 
		      }
			  if ("`plotyrange'"!="") {
		         if ("`plotxrange'"=="") local plotcond `plotcond' CI_l>=`min_yr'
				 else                    local plotcond `plotcond' &CI_l>=`min_yr'
			     if ("`max_yr'"!="") local plotcond `plotcond' &CI_r<=`max_yr' 
		      }
		   }

		   local plotcmdby `plotcmdby' (rcap CI_l CI_r CI_x ///
		                   `plotcond' in `ci_first'/`ci_last', ///
						   sort lcolor(`col') lpattern(`lty') `ciplotopt')

	    }
	 
	    *******************************
	    ***** Confidence Band *********
	    *******************************
	    tempname cval
	    scalar `cval'=.
	    if ("`cbON'"=="T") {
		   if (`nsims'<2000|`simsgrid'<50) {
	         di as text "Note: Setting at least nsims(2000) and simsgrid(50) is recommended to obtain the final results."
	       }
		   * Prepare grid for plotting
	       local cb_first=`byfirst'
		   local cb_last=`byfirst'-1+`cb_nr'
		
		   * fitting
		   tempname cb_b cb_V
		   capture confirm matrix `ci_b' `ci_V'
		   if ("`cb_p'"=="`ci_p'"& "`cb_s'"=="`ci_s'" & _rc==0) {
		      matrix `cb_b'=`ci_b'
		      matrix `cb_V'=`ci_V'
		   } 
		   else {
		      capture confirm matrix `line_b' `line_V'
		      if ("`cb_p'"=="`line_p'"& "`cb_s'"=="`line_s'" & _rc==0) {
		         matrix `cb_b'=`line_b'
		         matrix `cb_V'=`line_V'
		      }
			  else {
			     capture confirm matrix `dots_b' `dots_V'
			     if ("`cb_p'"=="`dots_p'"& "`cb_s'"=="`dots_s'" & _rc==0) {
		            matrix `cb_b'=`dots_b'
		            matrix `cb_V'=`dots_V'
		         }
				 else {
				    binslogit_fit `y_var' `x_var' `w_var' `conds' `wt', deriv(`deriv') ///
	                    p(`cb_p') s(`cb_s') type(cb) `vce' ///
			            xcat(`xcat') kmat(`kmat') dotsmean(0) /// 
			            xname(`xsub') yname(`ysub') catname(`xcatsub') edge(`binedges') ///
						`sorted' `usegtools' logitopt(`logitopt')
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
			      binsreg_pval(`cb_basis', `Xm'[,2], "`cb_V'", ".", `nsims', `nseries', "two", `=`level'/100', ".", "`cval'", "inf")		  
		   mata: mata drop `cb_basis' `Xm'
	          
	       * prediction
		   mata: `plotmatby'[|1,`cb_start' \ `cb_nr',`cb_end'|] =   ///
			                binslogit_plotmat("`cb_b'", "`cb_V'",      ///
						                      `=`cval'', "`kmat'",       ///
		                                      `nbins', `cb_p', `cb_s', `deriv', ///
							                  "cb", `cbngrid', "`wval'", `nwvar', ///
											  "`transform'", "`asyvar'")	  
							  
		   * cb
		   local plotnum=`plotnum'+1
		   local col: word `counter_by' of `bycolors'
	       local plotcond ""
		   if ("`plotxrange'"!=""|"`plotyrange'"!="") {
		      local plotcond if
		      if ("`plotxrange'"!="") {
	             local plotcond `plotcond' CB_x>=`min_xr'
			     if ("`max_xr'"!="") local plotcond `plotcond' &CB_x<=`max_xr' 
		      }
			  if ("`plotyrange'"!="") {
		         if ("`plotxrange'"=="") local plotcond `plotcond' CB_l>=`min_yr'
				 else                    local plotcond `plotcond' &CB_l>=`min_yr'
			     if ("`max_yr'"!="") local plotcond `plotcond' &(CB_r<=`max_yr'|CB_r==.) 
		      }
		   }

		   local plotcmdby (rarea CB_l CB_r CB_x ///
		                   `plotcond' in `cb_first'/`cb_last', sort cmissing(n) ///
		                   lcolor(none%0) fcolor(`col'%50) fintensity(50) `cbplotopt') `plotcmdby'
	    }
		mat `cvallist'=(nullmat(`cvallist') \ `cval')
		
		local plotcmd `plotcmd' `plotcmdby'
		
		**************************
		* generate breakpoints
		local kmat_end1=2
		local kmat_end2=rowsof(`kmat')
		if ("`kmat_end2'"=="`nbins'") {
		   local kmat_end1=1
		   local kmat_end2=`nbins'
        }
		mata: `plotmatby'[|1,`ncolplot'+1 \ `nbins',`ncolplot'+3|] = ///
		                 (range(1,`nbins',1), /// 
						  st_matrix("`kmat'")[|1, 1 \ `nbins', 1|], ///
						  st_matrix("`kmat'")[|`kmat_end1', 1 \ `kmat_end2', 1|])
		* append graph data
		mata: `plotmat'=(`plotmat' \ `plotmatby')
		
		*********************************
	    **** display ********************
	    *********************************	 
	    di ""
	    * Plotting
	    if ("`plot'"=="") {
		   if (`counter_by'==1) {
		      di in smcl in gr "Binscatter plot, logit model"
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
	       di in smcl in gr "{lalign 1:Bin/Degree selection:}"                   _col(30) " {c |} "
		   if ("`selection'"=="P") {
		   di in smcl in gr "{ralign 29:Degree of polynomial}"         _col(30) " {c |} " _col(32) as result %7.0f `binsp'
	       di in smcl in gr "{ralign 29:# of smoothness constraints}"  _col(30) " {c |} " _col(32) as result %7.0f `binss'
		   }
		   else {
	       di in smcl in gr "{ralign 29:Degree of polynomial}"         _col(30) " {c |} " _col(32) as result %7.0f `dots_p'
	       di in smcl in gr "{ralign 29:# of smoothness constraints}"  _col(30) " {c |} " _col(32) as result %7.0f `dots_s'
		   }

		   di in smcl in gr "{ralign 29:# of bins}"                       _col(30) " {c |} " _col(32) as result %7.0f `nbins'
	       if ("`binselectmethod'"!="User-specified") {
		      if ("`binsmethod'"=="ROT") {
		         di in smcl in gr "{ralign 29:imse, bias^2}"                 _col(30) " {c |} " _col(32) as result %7.3f `=`mat_imse_bsq_rot'[`counter_by',1]'
	             di in smcl in gr "{ralign 29:imse, var.}"                   _col(30) " {c |} " _col(32) as result %7.3f `=`mat_imse_var_rot'[`counter_by',1]'
			  }
			  else if ("`binsmethod'"=="DPI") {
			     di in smcl in gr "{ralign 29:imse, bias^2}"                 _col(30) " {c |} " _col(32) as result %7.3f `=`mat_imse_bsq_dpi'[`counter_by',1]'
	             di in smcl in gr "{ralign 29:imse, var.}"                   _col(30) " {c |} " _col(32) as result %7.3f `=`mat_imse_var_dpi'[`counter_by',1]'
			  }
		   }
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
	 *************** Plotting ******************
	 *******************************************
	 drop _all
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
		qui gen binid=. in 1
		qui gen lef_ep=. in 1
		qui gen rig_ep=. in 1

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
	    local graphcmd twoway `plotcmd', xtitle(`x_varname') ytitle(`y_varname') xscale(range(`xsc')) `plot_legend' `options' 
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
		label var binid "Bin: ID"
		label var lef_ep "Bin: left endpoint"
		label var rig_ep "Bin: right endpoint"
	    qui save `"`savedata'"', `replace'
	 }
	 ***************************************************************************
	 	 
	 *********************************
	 ********** Return ***************
	 *********************************
	 ereturn clear
	 * # of observations
	 ereturn scalar N=`Ntotal'
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
	 * by group:
	 *ereturn matrix knot=`kmat'
	 ereturn matrix cval_by=`cvallist'
	 ereturn matrix nbins_by=`nbinslist'
	 ereturn matrix Nclust_by=`Nclustlist'
	 ereturn matrix Ndist_by=`Ndistlist'
	 ereturn matrix N_by=`Nlist'

	 ereturn matrix imse_var_rot=`mat_imse_var_rot'
	 ereturn matrix imse_bsq_rot=`mat_imse_bsq_rot'
	 ereturn matrix	imse_var_dpi=`mat_imse_var_dpi'
	 ereturn matrix imse_bsq_dpi=`mat_imse_bsq_dpi'
end

* Helper commands
* Estimation
program define binslogit_fit, eclass
     version 13
     syntax varlist(min=2 numeric ts fv) [if] [in] [fw aw pw] [, deriv(integer 0) ///
	        p(integer 0) s(integer 0) type(string) vce(passthru)  ///
			xcat(varname numeric) kmat(name) dotsmean(integer 0) ///        /* xmean: report x-mean? */
			xname(name) yname(name) catname(name) edge(name) ///
			usereg sorted usegtools logitopt(string asis)]                  /* usereg: force the command to use reg; sored: sorted data? */
	 
	 preserve
	 marksample touse
	 qui keep if `touse'
	 
	 if ("`weight'"!="") local wt [`weight'`exp']
	 
	 tokenize `varlist'
	 local y_var `1'
	 local x_var `2'
	 macro shift 2
	 local w_var "`*'"
	 local nbins=rowsof(`kmat')-1
	 
	 tempname matxmean temp_b temp_V
	 mat `matxmean'=.
	 mat `temp_b'=.
	 mat `temp_V'=.
	 
	 if (`dotsmean'!=0) {
	    if ("`sorted'"==""|"`weight'"!=""|"`usegtools'"!="") {
	       if ("`usegtools'"=="") {
		      tempfile tmpfile
			  qui save `tmpfile', replace
			  
		      collapse (mean) `x_var' `wt', by(`xcat') fast
		      mkmat `xcat' `x_var', matrix(`matxmean')
			  
			  use `tmpfile', clear
		   }
		   else {
		      tempname obj
		      qui gstats tabstat `x_var' `wt', stats(mean) by(`xcat') matasave("`obj'")
			  mata: st_matrix("`matxmean'", (`obj'.getnum(.,1), `obj'.getOutputVar("`x_var'")))	
		      mata: mata drop `obj'
		   }
		}
		else {
		   tempname output
		   mata: `output'=binsreg_stat(`xname', `catname', `nbins', `edge', "mean", -1); ///
				 st_matrix("`matxmean'", `output')
		   mata: mata drop `output'
		}
	 }
	 
	 * Regression?
	 if (`p'==0) {
		 capture logit `y_var' ibn.`xcat' `w_var' `wt', nocon `vce' `logitopt'
		 if (_rc==0) {
		    matrix `temp_b'=e(b)
		    matrix `temp_V'=e(V)
		 }
		 else {
		    error _rc
			exit _rc
		 }
	 }
	 else {
	     local nseries=(`p'-`s'+1)*(`nbins'-1)+`p'+1
	     local series ""
	     forvalues i=1/`nseries' {
	        tempvar sp`i'
	        local series `series' `sp`i''
		    qui gen `sp`i''=. in 1
	     }
	
		 mata: binsreg_st_spdes(`xname', "`series'", "`kmat'", `catname', `p', 0, `s')
		 
		 capture logit `y_var' `series' `w_var' `wt', nocon `vce' `logitopt'
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
	 
	 
	 ereturn clear
	 ereturn matrix bmat=`temp_b'
	 ereturn matrix Vmat=`temp_V'
	 ereturn matrix xmat=`matxmean'         /* xcat, xbar */
end

mata:

  // Prediction for plotting
  real matrix binslogit_plotmat(string scalar eb, string scalar eV, real scalar cval, ///
                                string scalar knotname, real scalar J, ///
                                real scalar p, real scalar s, real scalar deriv, ///
                                string scalar type, real scalar ngrid, string scalar muwmat, ///
							    real scalar nw, string scalar transform, string scalar avar, | string scalar muxmat) 
  {
    real matrix coef, bmat, rmat, vmat, knot, xmean, wval, eval, out, fit, fit0, se, semat, Xm, Xm0, result
	real scalar nseries
	
	nseries=(p-s+1)*(J-1)+p+1
	coef=st_matrix(eb)'
	bmat=coef[|1\nseries|]
	if (nw>0) rmat=coef[|(nseries+1)\rows(coef)|]
	
	if (type=="ci"|type=="cb") {
	   vfull=st_matrix(eV)
	   vmat=vfull[|1,1\nseries,nseries|]
	}
    
    // Prepare evaluation points
    eval=J(0,3,.)
    if (args()==15) {
	   xmean=st_matrix(muxmat)
       eval=(eval \ (xmean[,2], J(J, 1, 0), xmean[,1])) 
    }
    if (ngrid!=0) {
	   eval=(eval \ binsreg_grids(knotname, ngrid))
    }
	
	// adjust w variables
	if (nw>0) {
	   wvec=st_matrix(muwmat)
       wval=wvec*rmat
	}
	else wval=0

	fit=J(0,1,.)
	se=J(0,1,.)
	if (p==0) {
	   if (args()==15) fit=(fit \ bmat)
	   if (ngrid!=0) {
	     fit=(fit \ (bmat#(J(ngrid,1,1)\.)))
		 fit=fit[|1 \ (rows(fit)-1)|]
	   }
	   if (type=="ci"|type=="cb") {
	      if (avar=="on") semat=sqrt(diagonal(vmat))
		  else {
		     if (nw>0) {
			    Xm=(I(nseries), J(nseries,1,1)#wvec)
				semat=sqrt(rowsum((Xm*vfull):*Xm))
			 }
			 else semat=sqrt(diagonal(vmat))
		  }
		  if (args()==15) se=(se \ semat)
		  if (ngrid!=0) {
		     se=(se \ (semat#(J(ngrid,1,1)\.)))
		     se=se[|1 \ (rows(se)-1)|]
		  }
	   }
	   if (type=="dots"|type=="line") {
	      if (transform=="T") out=(eval, logistic(fit:+wval))
		  else                out=(eval, fit:+wval)
	   }
	   else {
	      if (transform=="T") out=(eval, logistic(fit:+wval)-(logisticden(fit:+wval):*se)*cval, ///
		                                 logistic(fit:+wval)+(logisticden(fit:+wval):*se)*cval)
		  else                out=(eval, (fit:+wval)-se*cval, (fit:+wval)+se*cval)
	   }
	}
	else {
	   Xm=binsreg_spdes(eval[,1], knotname, eval[,3], p, deriv, s)
	   if (type=="dots"|type=="line") {
	      if (transform=="T") {
		     fit=binsreg_pred(Xm, bmat, ., "xb")[,1]
		     if (deriv==0) {
			    fit=logistic(fit:+wval)
			 }
			 if (deriv==1) {
			    Xm0=binsreg_spdes(eval[,1], knotname, eval[,3], p, 0, s)
			    fit0=binsreg_pred(Xm0, bmat, ., "xb")[,1]
			    fit=logisticden(fit0:+wval):*fit
			 }
			 out=(eval, fit)
		  }
		  else {
	         fit=binsreg_pred(Xm, bmat, ., "xb")[,1]
			 if (deriv==0) out=(eval, fit:+wval)
			 else          out=(eval, fit)
		  }
	   }
	   else {
	      if (avar=="on") {
	   	     result=binsreg_pred(Xm, bmat, vmat, "all")
	         if (transform=="T") {
		  	    Xm0=binsreg_spdes(eval[,1], knotname, eval[,3], p, 0, s)
			    fit0=binsreg_pred(Xm0, bmat, ., "xb")[,1]
		        result[,2]=logisticden(fit0:+wval):*result[,2]
			 
		        if (deriv==0) {
			       result[,1]=logistic(result[,1]:+wval)
			    }
			    else if (deriv==1) {
			       result[,1]=logisticden(fit0:+wval):*result[,1]
			    }
			 
			    out=(eval, result[,1]-cval*result[,2], result[,1]+cval*result[,2])
		     }
		     else {
		        if (deriv==0) out=(eval, (result[,1]:+wval)-cval*result[,2], (result[,1]:+wval)+cval*result[,2])
			    else          out=(eval, result[,1]-cval*result[,2], result[,1]+cval*result[,2])
	         }
		  }
		  else {
		     result=binsreg_pred(Xm, bmat, vmat, "all")
		     if (transform=="T") {				
				if (deriv==0) {
				   if (nw>0) Xm=(Xm, J(rows(Xm),1,1)#wvec)
				   result[,2]=logisticden(result[,1]:+wval):*sqrt(rowsum((Xm*vfull):*Xm))
				   result[,1]=logistic(result[,1]:+wval)
				}
				if (deriv==1) {
				   Xm0=binsreg_spdes(eval[,1], knotname, eval[,3], p, 0, s)
				   if (nw>0) {
				      Xm0=(Xm0, J(rows(Xm0),1,1)#wvec)
					  Xm=(Xm, J(rows(Xm),nw,0))
				   }
			       fit0=binsreg_pred(Xm0, coef, ., "xb")[,1]
				   Xm=logisticden(fit0):*(1:-2*logistic(fit0)):*result[,1]:*Xm0 + ///
				      logisticden(fit0):*Xm
				   result[,2]=sqrt(rowsum((Xm*vfull):*Xm))
				   result[,1]=logisticden(fit0):*result[,1]
				}
				out=(eval, result[,1]-cval*result[,2], result[,1]+cval*result[,2])
			 }
			 else {
			    if (nw>0) {
				   if (deriv==0) Xm=(Xm, J(rows(Xm),1,1)#wvec)
				   else          Xm=(Xm, J(rows(Xm),nw,0))
				}
				result=binsreg_pred(Xm, coef, vfull, "all")
				out=(eval, result[,1]-cval*result[,2], result[,1]+cval*result[,2])
			 }
		  }
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

  
end


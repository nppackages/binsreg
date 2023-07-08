*! version 1.3 03-Jul-2023

capture program drop binsregselect
program define binsregselect, eclass
     version 13
	 
	 syntax varlist(min=2 numeric ts fv) [if] [in] [fw aw pw] [, deriv(integer 0) ///
	        absorb(string asis) reghdfeopt(string asis) ///
			bins(numlist integer max=2 >=0) pselect(numlist integer >=0) sselect(numlist integer >=0) ///
			binspos(string) nbins(string) ///
			binsmethod(string) nbinsrot(string) ///
			simsgrid(integer 20) savegrid(string asis) replace ///
			dfcheck(numlist integer max=2 >=0) masspoints(string) usegtools(string) ///
			vce(passthru) useeffn(string) randcut(numlist max=1 >=0 <=1) ///
			norotnorm numdist(string) numclust(string)]   
			/* last line only for internal use */
			
	 set more off
	 
	 **************************************
	 ******** Regularization constant  ****
	 **************************************
	 local qrot=2
	 local rot_lb=1
	 local den_alp=0.975
	 
	 **************************************
	 * Create weight local
     if ("`weight'"!="") {
	    local wt [`weight'`exp']
		local wtype=substr("`weight'",1,1)
	 }
	 	 
	 **********************
	 ** Extract options ***
	 **********************	 
	 * default vce
	 if ("`vce'"=="") local vce "vce(robust)"
	 
	 * binning
	 * indictors: selectJ (F means select p)
	 local selectJ ""
	 *if ("`nbins'"=="F") local nbins ""
	 if ("`nbins'"=="T"|"`nbins'"=="") local nbins=0     /* default: select J */
	 local len_nbins=0
	 if ("`nbins'"!="") {
	    numlist "`nbins'", integer range(>=0) sort
	    local nbins=r(numlist)
	    local len_nbins: word count `nbins'
	 }
	 
	 if ("`nbins'"=="0"|`len_nbins'>1|"`bins'"!="") local selectJ "T"   
	 
	 * analyze bin- and order-related options
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
	 	 
	 if ("`bins'"!="") {
	    tokenize `bins'        /* overwrite pselect and sselect */
	    local plist "`1'"
	    local slist "`2'"
		if ("`plist'"=="") local plist=`deriv'
		if ("`slist'"=="") local slist=`plist'
     }
	 
	 local len_p: word count `plist'
	 local len_s: word count `slist'
	 
	 
	 if ((`len_p'==1&`len_s'==0)|(`len_p'==0&`len_s'==1)|(`len_p'==1&`len_s'==1)) {
	    local selectJ "T"
	 }
	 
	 if ("`selectJ'"=="T") {
	 	if (`len_p'>1|`len_s'>1) {
		   di as error "Only one p and one s are allowed."
		   exit
		}
	 	if ("`plist'"=="") local plist=`deriv'
		if ("`slist'"=="") local slist=`plist'
     }
	 
	 local len_p: word count `plist'
	 local len_s: word count `slist'
	 
	 if ((`len_p'>1|`len_s'>1) & "`selectJ'"!="T") {
	    local selectJ "F"   /* select p and s */
	 }
	 
	 if ("`selectJ'"=="") {
	    di as error "Degree, smoothness, or # of bins are not correctly specified."
        exit
	 }
	
	 * find all compatible pairs of p and s
	 tempname deg_mat
	 if ("`selectJ'"=="F") {
	    if (`len_p'>0 & `len_s'==0) {
		   mat `deg_mat'=J(`len_p', 2, .)
		   forval i=1/`len_p' {
		      local el : word `i' of `plist'
		      mat `deg_mat'[`i',1]=`el'
			  mat `deg_mat'[`i',2]=`el'
	       }
		}
	    if (`len_p'==0 & `len_s'>0) {
		   mat `deg_mat'=J(`len_s', 2, .)
		   forval i=1/`len_s' {
		      local el : word `i' of `slist'
		      mat `deg_mat'[`i',1]=`el'
			  mat `deg_mat'[`i',2]=`el'
	       }
		}
		if (`len_p'>0 & `len_s'>0) {
		   mat `deg_mat'=J(`=`len_p'*`len_s'',2,.)
		   local ncom=0
		   forval i=1/`len_p' {
		      local el_p : word `i' of `plist'
			  forval j=1/`len_s' {
		         local el_s : word `j' of `slist'
		         if (`el_p'>=`el_s') {
				    mat `deg_mat'[`=`ncom'+1',1]=`el_p'
					mat `deg_mat'[`=`ncom'+1',2]=`el_s'
					local ++ncom
			     }
	          }
	       }
		   if (`ncom'>0) mat `deg_mat'=`deg_mat'[1..`ncom', 1..2]
		   else {
		      di as error "degree and smoothness incompatible"
			  exit
		   }
		}
	 }
	 else {
	    mat `deg_mat'=(`plist', `slist')
	 }
	 
	 * take submatrix with p>=deriv
	 local ncom=0
	 local index ""
	 tempname m_deg                       /* degree matrix to be used */
	 forval i=1/`=rowsof(`deg_mat')' {
	    if (`deg_mat'[`i',1]>=`deriv') {
		   local ++ncom
		   mat `m_deg'=(nullmat(`m_deg') \ `deg_mat'[`i',1..2])
		}
	 }
	 if (`ncom'==0) {
	    di as error "Degree and smoothness incorrectly specified."
		exit
	 }
	 
	 if ("`binspos'"=="") local binspos "QS"
	 if ("`binspos'"=="es") local binspos "ES"
	 if ("`binspos'"=="qs") local binspos "QS"
	 if ("`binsmethod'"=="") local binsmethod "DPI"
	 if ("`binsmethod'"=="rot") local binsmethod "ROT"
	 if ("`binsmethod'"=="dpi") local binsmethod "DPI"
	 if ("`dfcheck'"=="") local dfcheck 20 30
	 
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
	    di as text in gr "Warning: masspoints(veryfew) not allowed for bin selection."
		local localcheck "F"
		local rot_fewobs "T"
		local dpi_fewobs "T"
	 }
	 
	 tokenize `dfcheck'
	 local dfcheck_n1 "`1'"
	 local dfcheck_n2 "`2'"
	 
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
	    * use gstats tab instead of tabstat/collapse
		* use gquantiles instead of _pctile
		* use gunique instead of binsreg_uniq
		* use fasterxtile instead of irecode (within binsreg_irecode)
		* shut down local checks & do not sort
	 }
	 
	 * cluster var?	 
	 local vcetemp: subinstr local vce "vce(" "", all
	 local vcetemp: subinstr local vcetemp ")" "", all
	 tokenize "`vcetemp'", parse(", ")
	 if ("`1'"=="cl"|"`1'"=="clu"|"`1'"=="clus"|"`1'"=="clust"|"`1'"=="cluste"|"`1'"=="cluster") {
	    if ("`3'"==""|"`3'"==",") local clusterON "T"
		local clustervar `2'  /* only keep the 1st cluster var */
     }
	
	 * use reghdfe?
	 if ("`absorb'"!="") {
	    capture which reghdfe
		if (_rc) {
		   di as error "reghdfe not installed."
		   exit
		}
	 }
	 
	 *****************************************************
	 * Error checks
	 if (`deriv'<0) {
	    di as error "deriv() incorrectly specified."
	 	exit
	 }
	 if (`simsgrid'<0) {
	    di as error "simsgrid() incorrectly specified."
		exit
	 }
	 if (`"`savegrid'"'!=`""'&"`replace'"=="") {
	    confirm new file `"`savegrid'.dta"'
	 }
	 if ("`nbinsrot'"!="") {
	    confirm integer n `nbinsrot'
	 }
	 
	 * Mark sample
	 preserve
	 
	 * Parse varlist into y_var, x_var and w_var
	 tokenize `varlist'
	 fvrevar `1', tsonly
	 local y_var "`r(varlist)'"
	 fvrevar `2', tsonly
	 local x_var "`r(varlist)'"
	 fvrevar `2', list
	 local x_varname "`r(varlist)'"
	 
	 macro shift 2
	 local w_var "`*'"
	 fvrevar `w_var', list
	 local w_varname "`r(varlist)'"
	 fvrevar `w_var', tsonly
	 local w_var "`r(varlist)'"       /* so time series operator respected */
	 
	 marksample touse      /* now renew the marker to account for missing values */
	 qui keep if `touse'
	 local eN=_N
	 local samplesize=_N
	 
	 if ("`usegtools'"==""&("`masspoints'"!="off"|"`binspos'"=="QS")) {
	    if ("`:sortedby'"!="`x_var'") sort `x_var', stable
	 }
	 
	 * Normalize support
	 tempvar z_var
	 if ("`wtype'"=="f") qui sum `x_var' `wt', meanonly
	 else                qui sum `x_var', meanonly
	 
	 local N=r(N)      /* sample size, with wt */
	 local xmin=r(min)
	 local xmax=r(max)
	 
	 tempname xvec zvec Xm binedges
		 
	 * Extract effective sample size
	 local Ndist=.
	 if ("`massadj'"=="T") {
	    if ("`numdist'"!=""&"`numdist'"!=".") {
		   local Ndist=`numdist'
		}
		else {
		   if ("`usegtools'"=="") {
	          mata: `binedges'=binsreg_uniq(st_data(.,"`x_var'"), ., 1, "Ndist")
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
		local eN=min(`eN', `Nclust')
	 }
	 
	 * Take a subsample?
	 if ("`randcut'"!="") {
	    mata: `xvec'=st_data(.,"`x_var'")
		qui keep if runiform()<=`randcut'
		local eN_sub=_N
		
		local Ndist_sub=.
		if ("`massadj'"=="T") {
		   if ("`usegtools'"=="") {
	          mata: `binedges'=binsreg_uniq(st_data(.,"`x_var'"), ., 1, "Ndist_sub")
		      mata: mata drop `binedges'
		   }
		   else {
		      qui gunique `x_var'
			  local Ndist_sub=r(unique)
		   }
		   local eN_sub=min(`eN_sub', `Ndist_sub')
	    }
	    
		local Nclust_sub=.
	    if ("`clusterON'"=="T") {
		   if ("`usegtools'"=="") {
		       mata: st_local("Nclust_sub", strofreal(rows(uniqrows(st_data(.,"`clustervar'")))))
		   }
		   else {
		       qui gunique `clustervar'
			   local Nclust_sub=r(unique)
		   }
		   local eN_sub=min(`eN_sub', `Nclust_sub')
	    }
	 }
	 else {
		local eN_sub=`eN'
		local Ndist_sub=`Ndist'
		local Nclust_sub=`Nclust'
	 }
	 
	 sum `x_var', meanonly
	 gen `z_var'=(`x_var'-`=r(min)')/(`=r(max)'-`=r(min)')
	 mata: `zvec'=st_data(., "`z_var'")   /* normalized x, subsample */

	 
	 * Define matrices here to save results
	 tempname mat_imse_bsq_rot mat_imse_var_rot mat_imse_bsq_dpi mat_imse_var_dpi ///
	          mat_J_rot_unreg mat_J_rot_reg mat_J_rot_uniq mat_J_dpi mat_J_dpi_uniq
	 mat `mat_imse_bsq_rot'=J(`ncom',1,.)
	 mat `mat_imse_var_rot'=J(`ncom',1,.)
	 mat `mat_imse_bsq_dpi'=J(`ncom',1,.)
	 mat `mat_imse_var_dpi'=J(`ncom',1,.) 
	 mat `mat_J_rot_unreg'=J(`ncom',1,.)
	 mat `mat_J_rot_reg'=J(`ncom',1,.)
	 mat `mat_J_rot_uniq'=J(`ncom',1,.) 
	 mat `mat_J_dpi'=J(`ncom',1,.)
	 mat `mat_J_dpi_uniq'=J(`ncom',1,.)
	 
	 
	 ****** START loop here **************
	 forval num=1/`ncom' {
	    
		* extract p and s from the matrix
		local p=`m_deg'[`num', 1]
		local s=`m_deg'[`num', 2]
	    
		* prepare locals for reporting
	    local imse_bsq_rot=.
	    local imse_var_rot=.
	    local imse_bsq_dpi=.
	    local imse_var_dpi=.
	 
		***************************
		******* ROT choice ********
		***************************
		tempname vcons bcons coef
		tempvar resid1 resid2   /* only used by reghdfe */
		local J_rot_reg=.
		local J_rot_unreg=.
		if ("`nbinsrot'"!="") local J_rot_reg=`nbinsrot'

		* Initial checking of sample size (for ROT)
		if (`J_rot_reg'==.&"`rot_fewobs'"=="") {
			if (`eN_sub'<=`dfcheck_n1'+`p'+1+`qrot') {
				local rot_fewobs "T"
				di as text in gr "Warning: Too small effective sample size for bin selection."
			}
		}

		if ("`rot_fewobs'"!="T"&`J_rot_reg'==.) {
			* Power series
			local series_rot ""
			forvalues i=1/`=`p'+`qrot'' {
				tempvar z_var_`i'
				qui gen `z_var_`i''=`z_var'^`i'
				local series_rot `series_rot' `z_var_`i''
			}

			* Variance Component
			if ("`absorb'"=="") capture reg `y_var' `series_rot' `w_var' `wt'
			else                capture reghdfe `y_var' `series_rot' `w_var' `wt', absorb(`absorb') resid(`resid1')

			if (_rc==0) {
				mat `coef'=e(b)
				mat `coef'=`coef'[1,`=`p'+1'..`=`p'+`qrot'']
			}
			else {
				error _rc
				exit _rc
			}

			tempvar pred_y y_var_2 pred_y2 s2
			if ("`absorb'"=="") predict `pred_y', xb
			else                predict `pred_y', xbd

			qui gen `y_var_2'=`y_var'^2                                              // move it outside
			if ("`absorb'"=="") {
				capture reg `y_var_2' `series_rot' `w_var' `wt'
				if (_rc) {
					error _rc
					exit _rc
				}
				predict `pred_y2', xb
			}
			else {
				capture reghdfe `y_var_2' `series_rot' `w_var' `wt', absorb(`absorb') resid(`resid2')
				if (_rc) {
					error _rc
					exit _rc
				}

				predict `pred_y2', xbd
			}

			qui gen `s2'=`pred_y2'-`pred_y'^2         /* sigma^2(x) var */

			* Normal density
			if ("`rotnorm'"=="") {
				if ("`wtype'"!="p") qui sum `z_var' `wt'
				else                qui sum `z_var' [aw`exp']
				local zbar=r(mean)
				local zsd=r(sd)

				tempvar fz
				* trim density from below
				local cutval=normalden(invnormal(`den_alp')*`zsd', 0, `zsd')
				qui gen `fz'=max(normalden(`z_var', `zbar', `zsd'), `cutval')	 
				if ("`binspos'"=="ES") qui replace `s2'=`s2'/`fz' 
				else                   qui replace `s2'=`s2'*(`fz'^(2*`deriv'))
			}

			if ("`wt'"!="") qui sum `s2' [aw`exp'], meanonly
			else            qui sum `s2', meanonly
			local sig2=r(mean)    
			mata: imse_v_cons(`p', `deriv', "`vcons'")
			local imse_v=`sig2'*`vcons'                  /* variance constant */

			* Bias component
			* gen data for derivative
			tempvar pred_deriv
			mata: `Xm'=J(rows(`zvec'),0,.) 

			forval i=`=`p'+1'/`=`p'+`qrot'' {
				mata:`Xm'=(`Xm',`zvec':^(`i'-`p'-1)*factorial(`i')/factorial(`i'-`p'-1))
			}

			mata: `Xm'=`Xm'*st_matrix("`coef'")'; ///
				  st_store(.,st_addvar("float", "`pred_deriv'"), `Xm':^2)

			mata: mata drop `Xm'

			if ("`rotnorm'"=="") {
				if ("`binspos'"=="QS") {
					qui replace `pred_deriv'=`pred_deriv'/(`fz'^(2*`p'+2-2*`deriv'))
				}
			}
			if ("`wt'"!="") qui sum `pred_deriv' [aw`exp'], meanonly
			else            qui sum `pred_deriv', meanonly
			local mean_deriv=r(mean)

			mata: imse_b_cons(`p', `deriv', "`bcons'")
			local imse_b=`mean_deriv'*`bcons'          /* bias constant */

			* ROT J
			local J_rot_unreg=ceil((`imse_b'*2*(`p'+1-`deriv')/               ///
				(`imse_v'*(1+2*`deriv')))^(1/(2*`p'+2+1))* ///
				`eN_sub'^(1/(2*`p'+2+1)))
			local J_rot_reg=max(`J_rot_unreg', ///
				  ceil((2*(`p'+1-`deriv')/(1+2*`deriv')*`rot_lb'*`eN_sub')^(1/(2*`p'+2+1)))) 

			local imse_bsq_rot=`imse_b'
			local imse_var_rot=`imse_v'
			
			mat `mat_imse_bsq_rot'[`num',1]=`imse_bsq_rot'
			mat `mat_imse_var_rot'[`num',1]=`imse_var_rot'
		}

		** Repeated knots? ***************
		local J_rot_uniq=`J_rot_reg'

		if (("`binsmethod'"=="DPI"|"`localcheck'"=="T")&"`masspoints'"!="veryfew") {
			tempvar zcat
			qui gen `zcat'=. in 1
			* Prepare bins
			tempname kmat

			if "`binspos'"=="ES" {
				local stepsize=1/`J_rot_reg'
				forvalues i=1/`=`J_rot_reg'+1' {
					mat `kmat'=(nullmat(`kmat') \ `=0+`stepsize'*(`i'-1)')
				}
			}
			else {
				if (`J_rot_reg'==1) {
					mat `kmat'=(0 \ 1)
				}
				else {
					binsreg_pctile `z_var' `wt', nq(`J_rot_reg') `usegtools'
					mat `kmat'=(0 \ r(Q) \ 1)
				}
			}

			mata: st_matrix("`kmat'", (0 \ uniqrows(st_matrix("`kmat'")[|2 \ `=`J_rot_reg'+1'|])))
			local J_rot_uniq=rowsof(`kmat')-1
			if ("`binsmethod'"=="DPI"&"`dpi_fewobs'"=="") {
				binsreg_irecode `z_var', knotmat(`kmat') bin(`zcat') ///
					`usegtools' nbins(`J_rot_uniq') pos(`binspos') knotliston(T)
			}
		}

		*********************************
		********** DPI Choice ***********
		*********************************
		local J_dpi=.	 
		* Check if DPI can be implemented
		if ("`J_rot_uniq'"!="."&"`binsmethod'"=="DPI"&"`masspoints'"!="veryfew") {
			* Compare with degree of freedom
			if ((`p'-`s'+1)*(`J_rot_uniq'-1)+`p'+2+`dfcheck_n2'>=`eN_sub') {
				di as text in gr  "Warning: Too small effective sample size for DPI selection."
				local dpi_fewobs "T"
			}

			* Check local effective size
			if ("`localcheck'"=="T"&"`dpi_fewobs'"!="T") {
			    mata: st_local("Ncat", strofreal(rows(uniqrows(st_data(.,"`zcat'")))))
				if (`J_rot_uniq'==`Ncat') {
		          mata: `binedges'=binsreg_uniq(`zvec', st_data(.,"`zcat'"), `J_rot_uniq', "uniqmin")
				  mata: mata drop `binedges'
		        }
				else {
		          local uniqmin=0
			      di as text in gr "Warning: There are empty bins."
		        } 
				
				if (`uniqmin'<`p'+1) {
					local dpi_fewobs "T"
					di as text in gr  "Warning: Some bins have too few distinct x-values for DPI selection."
				}
			}
		}
		else local dpi_fewobs "T"

		if ("`binsmethod'"=="DPI"&"`dpi_fewobs'"!="T") {	
			* Update vce condition
			if ("`massadj'"=="T") {
				if ("`absorb'"=="") {
					if ("`clusterON'"=="") {
						local vce "vce(cluster `z_var')"
					}
					else {
						if (`Nclust_sub'>`Ndist_sub') {
							local vce "vce(cluster `z_var')"
							di as text in gr "Warning: # of mass points < # of clusters. vce option overridden."
						}
					}
				}
				else {
					if ("`clustervar'"=="") local vce "vce(cluster `z_var')"
				}
			}

			**************************************
			* Start computation
			tempvar derivfit derivse biasterm biasterm_v projbias
			qui gen `derivfit'=. in 1
			qui gen `derivse'=. in 1
			qui gen `biasterm'=. in 1                     /* save bias */
			if (`deriv'>0) qui gen `biasterm_v'=. in 1    /* error of approx deriv */
			qui gen `projbias'=. in 1                     /* save proj of bias */

			**************************************
			* predict leading bias
			mata: bias("`z_var'", "`zcat'", "`kmat'", `p', 0, "`biasterm'")
			if (`deriv'>0) {
				mata: bias("`z_var'", "`zcat'", "`kmat'", `p', `deriv', "`biasterm_v'")
			}

			* Increase order from p to p+1
			* Expand basis
			local nseries=(`p'-`s'+1)*(`J_rot_uniq'-1)+`p'+2
			local series ""
			forvalues i=1/`nseries' {
				tempvar sp`i'
				local series `series' `sp`i''
				qui gen `sp`i''=. in 1
			}

			mata: binsreg_st_spdes(`zvec', "`series'", "`kmat'", st_data(.,"`zcat'"), `=`p'+1', 0, `=`s'+1')

			if ("`absorb'"=="") capture reg `y_var' `series' `w_var' `wt', nocon
			else                capture reghdfe `y_var' `series' `w_var' `wt', absorb(`absorb') `reghdfeopt'

			* store results
			tempname temp_b temp_V
			if (_rc==0) {
				matrix `temp_b'=e(b)
				matrix `temp_V'=e(V)
			}
			else {
				error  _rc
				exit _rc
			}

			* Predict (p+1)th derivative
			mata: `Xm'=binsreg_spdes(`zvec', "`kmat'", st_data(.,"`zcat'"), `=`p'+1', `=`p'+1', `=`s'+1'); ///
				st_store(.,"`derivfit'", (binsreg_pred(`Xm', (st_matrix("`temp_b'")[|1 \ `nseries'|])', ///
				st_matrix("`temp_V'")[|1,1 \ `nseries',`nseries'|], "xb"))[,1])
			mata: mata drop `Xm'

			qui replace `biasterm'=`derivfit'*`biasterm'
			if (`deriv'>0) qui replace `biasterm_v'=`derivfit'*`biasterm_v'
			drop `series'

			* Then get back degree-p spline, run OLS
			local nseries=(`p'-`s'+1)*(`J_rot_uniq'-1)+`p'+1
			local series ""
			forvalues i=1/`nseries' {
				tempvar sp`i'
				local series `series' `sp`i''
				qui gen `sp`i''=. in 1
			}

			mata: binsreg_st_spdes(`zvec', "`series'", "`kmat'", st_data(.,"`zcat'"), `p', 0, `s')   
			capture reg `biasterm' `series' `wt', nocon       /* project bias on X of degree p */
			tempname bias_b bias_V
			if (_rc==0) {
				matrix `bias_b'=e(b)
				matrix `bias_V'=e(V)
			}
			else {
				error _rc
				exit _rc
			}

			mata: `Xm'=binsreg_spdes(`zvec', "`kmat'", st_data(.,"`zcat'"), `p', `deriv', `s'); ///
				st_store(.,"`projbias'", binsreg_pred(`Xm', st_matrix("`bias_b'")', st_matrix("`bias_V'"), "xb")[,1])

			if (`deriv'==0) {
				qui replace `biasterm'=(`biasterm'-`projbias')^2
			}
			else {
				qui replace `biasterm'=(`biasterm_v'-`projbias')^2   /* still save in biasterm if deriv>0 */
			}

			if ("`wt'"!="") qui sum `biasterm' [aw`exp'], meanonly
			else            qui sum `biasterm', meanonly
			local m_bias=r(mean)
			local imse_b=`m_bias'*`J_rot_uniq'^(2*(`p'+1-`deriv'))

			* for variance purpose
			if ("`absorb'"=="") capture reg `y_var' `series' `w_var' `wt', nocon `vce'        
			else                capture reghdfe `y_var' `series' `w_var' `wt', absorb(`absorb') `vce' `reghdfeopt'

			* store results
			if (_rc==0) {
				matrix `temp_b'=e(b)
				matrix `temp_V'=e(V)
				tempname vcov
				mata: `vcov'=st_matrix("`temp_V'")
				if ("`absorb'"=="") mata: `vcov'=`vcov'[|1,1 \ `nseries',`nseries'|]
				else {
					mata: `vcov'=(`vcov'[|1,1 \ `nseries', `nseries'|], `vcov'[|1,cols(`vcov') \ `nseries', cols(`vcov')|] \ ///
						`vcov'[|cols(`vcov'), 1 \ cols(`vcov'), `nseries'|], `vcov'[cols(`vcov'), cols(`vcov')]); ///
						`Xm'=(`Xm', J(rows(`Xm'),1,1))
				}
			}
			else {
				error  _rc
				exit _rc
			}

			mata: st_store(., "`derivse'", (binsreg_pred(`Xm', ., `vcov', "se")[,2]):^2)
			mata: mata drop `vcov'

			if ("`wt'"!="") qui sum `derivse' [aw`exp'], meanonly
			else            qui sum `derivse', meanonly
			local m_se=r(mean)
			local imse_v=`m_se'/(`J_rot_uniq'^(1+2*`deriv'))	   
			mata: mata drop `Xm'   	   

			* DPI J
			local J_dpi=ceil((`imse_b'*2*(`p'+1-`deriv')/               ///
				(`imse_v'*(1+2*`deriv')))^(1/(2*`p'+2+1)))

			local imse_bsq_dpi=`imse_b'
			local imse_var_dpi=`imse_v'*`eN_sub'
			
			mat `mat_imse_bsq_dpi'[`num',1]=`imse_bsq_dpi'
			mat `mat_imse_var_dpi'[`num',1]=`imse_var_dpi'

		}
		local J_dpi_uniq=`J_dpi'



		************************************************
		* update J if useeffn or subsample specified
		if ("`useeffn'"!=""|"`randcut'"!="") {
			if ("`useeffn'"!="") local scaling=(`useeffn'/`eN')^(1/(2*`p'+2+1))
			if ("`randcut'"!="") local scaling=(`eN'/`eN_sub')^(1/(2*`p'+2+1))

			if (`J_rot_unreg'!=.) local J_rot_unreg=ceil(`J_rot_unreg'*`scaling')
			if (`J_rot_reg'!=.)   local J_rot_reg=ceil(`J_rot_reg'*`scaling')
			if (`J_rot_uniq'!=.)  local J_rot_uniq=ceil(`J_rot_uniq'*`scaling')
			if (`J_dpi'!=.)       local J_dpi=ceil(`J_dpi'*`scaling')
			if (`J_dpi_uniq'!=.)  local J_dpi_uniq=ceil(`J_dpi_uniq'*`scaling')
		}
		
		mat `mat_J_rot_unreg'[`num',1]=`J_rot_unreg'
		mat `mat_J_rot_reg'[`num',1]=`J_rot_reg'
		mat `mat_J_rot_uniq'[`num',1]=`J_rot_uniq'
		mat `mat_J_dpi'[`num',1]=`J_dpi'
		mat `mat_J_dpi_uniq'[`num',1]=`J_dpi_uniq'
	 
	 }
	 
	 ****** END loop *******
	 tempname ord_rot_unreg ord_rot_reg ord_rot_uniq ord_dpi ord_dpi_uniq ///
	          ind_rot_unreg ind_rot_reg ind_rot_uniq ind_dpi ind_dpi_uniq
	 local imse_var_dpi_upd=.
	 local imse_bsq_dpi_upd=.
	 if ("`selectJ'"=="F") {
	    * output a row vector of p and s
		foreach name in "rot_unreg" "rot_reg" "rot_uniq" "dpi" "dpi_uniq" {
		   mata: findmindex("`mat_J_`name''", "`ind_`name''", `nbins', `ncom')
		   mat `ord_`name''=`m_deg'[`ind_`name'',1..2]
		   local J_`name'=`nbins'
		}
		if (`nbins'!=`=`mat_J_dpi'[`ind_dpi',1]') {
	       qui bins_imse `y_var' `z_var' `w_var' `wt', deriv(`deriv') ///
               p(`=`ord_dpi'[1,1]') s(`=`ord_dpi'[1,2]') nbins(`J_dpi') eN_sub(`eN_sub') ///
		       binspos(`binspos') `vce' `usegtools' ///
		       zvec(`zvec') absorb(`absorb') reghdfeopt(`reghdfeopt')
		   local imse_var_dpi_upd=e(imse_var)
		   local imse_bsq_dpi_upd=e(imse_bsq)
		}
	 }
	 else {
	    if (`len_nbins'>1) {
		   tempname m_nbins
		   forval i=1/`len_nbins' {
		      local el: word `i' of `nbins' 
		      mat `m_nbins'=(nullmat(`m_nbins') \ `el')
		   }
		   * output a scalar
		   foreach name in "rot_unreg" "rot_reg" "rot_uniq" "dpi" "dpi_uniq" {
		      mata: findmindex("`m_nbins'", "`ind_`name''", `=`mat_J_`name''[1,1]', `len_nbins')
			  local J_`name'=`m_nbins'[`ind_`name'',1]
		   }
		}
		foreach name in "rot_unreg" "rot_reg" "rot_uniq" "dpi" "dpi_uniq" {
		   mat `ord_`name''=`m_deg'[1,1..2]
		}
	 }
	 mata: mata drop `zvec'
	 
	 * Reconstruct knot list
	 tempname xkmat
	 
	 if ("`binsmethod'"=="ROT"&"`rot_fewobs'"!="T") {
	    local Jselected=`J_rot_uniq'
		local pselected=`ord_rot_uniq'[1,1]
		local sselected=`ord_rot_uniq'[1,2]
	 }
	 else if ("`binsmethod'"=="DPI"&"`dpi_fewobs'"!="T") {
	    local Jselected=`J_dpi'
		local pselected=`ord_dpi'[1,1]
		local sselected=`ord_dpi'[1,2]
	 }
	 else {
	    local Jselectfail "T"
	 }
	 
	 if ("`Jselectfail'"!="T"&"`useeffn'"=="") {
	    if ("`binspos'"=="ES") {
	       local stepsize=1/`Jselected'   
	       forvalues i=1/`=`Jselected'+1' {
		      mat `xkmat'=(nullmat(`xkmat') \ `=`xmin'+`stepsize'*(`i'-1)')
		   }
	    }
	    else {
		   if (`Jselected'>1) {
	          if ("`randcut'"!="") {
			     qui set obs `samplesize'
			     mata: st_store(., "`x_var'", `xvec')
				 mata: mata drop `xvec'
			  }
			  binsreg_pctile `x_var' `wt', nq(`Jselected')
			  mat `xkmat'=(`xmin'\ r(Q) \ `xmax')
		   }
		   else mat `xkmat'=(`xmin' \ `xmax')		   
	    }
		
		* Renew if needed
		mata: st_matrix("`xkmat'", (`xmin' \ uniqrows(st_matrix("`xkmat'")[|2 \ `=`Jselected'+1'|])))
	    if (`Jselected'!=rowsof(`xkmat')-1) {
		   local Jselected=rowsof(`xkmat')-1
		}
		if ("`binsmethod'"=="DPI") {
		   local J_dpi_uniq=`Jselected'
		}
	 }
	 else mat `xkmat'=.
	 
	 
	 if ("`binsmethod'"=="DPI") local method "IMSE-optimal plug-in choice"
	 else                       local method "IMSE-optimal rule-of-thumb choice"
	 if ("`selectJ'"=="T") local method "`method' (select # of bins)"
	 else                  local method "`method' (select degree and smoothness)"
	 
	 
	 if ("`binspos'"=="ES") {
	     local placement "Evenly-spaced"
	 }
	 else {
	     local placement "Quantile-spaced"
	 }
	 
	 * Save data?
	 if (`"`savegrid'"'!=`""') {
	    if ("`Jselectfail'"!="T"&"`useeffn'"=="") {
		   clear
	       local obs=`simsgrid'*`Jselected'+`Jselected'-1
		   qui set obs `obs'
		   
		   qui gen `x_varname'=. in 1
		   label var `x_varname' "Eval. point"
	       qui gen binsreg_isknot=. in 1
		   label var binsreg_isknot "Is the eval. point an inner knot"
		   qui gen binsreg_bin=. in 1
		   label var binsreg_bin "indicator of bin"
		   mata: st_store(., (1,2,3), binsreg_grids("`xkmat'", `simsgrid'))
		   foreach var of local w_varname {
		      qui gen `var'=0
	   	   }
		   
		   qui save `"`savegrid'"', `replace'
	    }
		else {
		   di as text in gr "Warning: Grid not saved. Selection fails or useeffn() is specified."
		}
	 }

	 * Display
	 di ""
	 di in smcl in gr "Bin selection for binscatter estimates"
	 di in smcl in gr "Method: `method'"
	 di in smcl in gr "Position: `placement'"
	 if (`"`savegrid'"'!=`""') {
	    di in smcl in gr `"Output file: `savegrid'.dta"' 
	 }
	 di ""
	 di in smcl in gr "{hline 28}{c TT}{hline 10}"
	 di in smcl in gr "{ralign 27:# of observations}"    _col(28) " {c |} " _col(30) as result %7.0f `N'
	 di in smcl in gr "{ralign 27:# of distince values}" _col(28) " {c |} " _col(30) as result %7.0f `Ndist'
	 di in smcl in gr "{ralign 27:# of clusters}"        _col(28) " {c |} " _col(30) as result %7.0f `Nclust'
     if ("`useeffn'"=="") {
	 di in smcl in gr "{ralign 27:eff. sample size}"     _col(28) " {c |} " _col(30) as result %7.0f `eN'
	 }
	 else {
	 di in smcl in gr "{ralign 27:eff. sample size}"     _col(28) " {c |} " _col(30) as result %7.0f `useeffn'
     }
	 
	 foreach name in "rot_unreg" "rot_reg" "rot_uniq" "dpi" "dpi_uniq" {
		local df_`name'=.
		if (`J_`name''!=.) local df_`name'=(`ord_`name''[1,1]-`ord_`name''[1,2]+1)*(`J_`name''-1)+`ord_`name''[1,1]+1
	 }
	 
	 if ("`selectJ'"=="F") {
	    local imse_bsq_rot=`mat_imse_bsq_rot'[`ind_rot_unreg',1]
		local imse_var_rot=`mat_imse_var_rot'[`ind_rot_unreg',1]
	    if ("`imse_var_dpi_upd'"!=".") {
		   local imse_bsq_dpi=`imse_bsq_dpi_upd'
		   local imse_var_dpi=`imse_var_dpi_upd'
		}
		else {
		   local imse_bsq_dpi=`mat_imse_bsq_dpi'[`ind_dpi',1]
		   local imse_var_dpi=`mat_imse_var_dpi'[`ind_dpi',1]
	 	}
		di in smcl in gr "{hline 28}{c +}{hline 10}"
	    di in smcl in gr "{ralign 27:# of bins}"  _col(28) " {c |} " _col(30) as result %7.0f `nbins'
        di in smcl in gr "{hline 28}{c BT}{hline 10}"
	    di ""
	    di in smcl in gr "{hline 14}{c TT}{hline 8}{c TT}{hline 7}{c TT}{hline 7}{c TT}{hline 15}{c TT}{hline 14}"
		di in smcl in gr "{rcenter 13: method}" _col(13) " {c |} " "{center 7: p}" ///
		                                        _col(22) "{c |}" "{rcenter 7: s}" ///
												_col(31) "{c |}" "{center 7: df}" ///
												_col(40) "{c |}" "{center 14: imse, bias^2}" ///
												_col(56) "{c |}" "{center 14: imse, var.}"
		di in smcl in gr "{hline 14}{c +}{hline 8}{c +}{hline 7}{c +}{hline 7}{c +}{hline 15}{c +}{hline 14}"
		di in smcl in gr "{rcenter 13: ROT-POLY}"  _col(13) " {c |} " as result %4.0f `ord_rot_unreg'[1,1]  ///
		                                           _col(23) in gr " {c |} " as result %4.0f `ord_rot_unreg'[1,2]  ///
		                                           _col(32) in gr "{c |}" as result %5.0f `df_rot_unreg' ///
			                                       _col(40) in gr "{c |}  " as result %7.3f `imse_bsq_rot' ///
												   _col(56) in gr "{c |}  " as result %7.3f `imse_var_rot'
		di in smcl in gr "{rcenter 13: ROT-REGUL}" _col(13) " {c |} " as result %4.0f `ord_rot_reg'[1,1] ///
		                                           _col(23) in gr " {c |} " as result %4.0f `ord_rot_reg'[1,2]  ///
		                                           _col(32) in gr "{c |}" as result %5.0f `df_rot_reg' ///
			                                       _col(40) in gr "{c |}  " as result %7.3f . ///                     
												   _col(56) in gr "{c |}  " as result %7.3f .
		di in smcl in gr "{rcenter 13: ROT-UKNOT}" _col(13) " {c |} " as result %4.0f `ord_rot_uniq'[1,1] ///
		                                           _col(23) in gr " {c |} " as result %4.0f `ord_rot_uniq'[1,2]  ///
                                                   _col(32) in gr "{c |}" as result %5.0f `df_rot_uniq' ///
			                                       _col(40) in gr "{c |}  " as result %7.3f . ///
												   _col(56) in gr "{c |}  " as result %7.3f .
		di in smcl in gr "{rcenter 13: DPI}"       _col(13) " {c |} " as result %4.0f `ord_dpi'[1,1] ///
		                                           _col(23) in gr " {c |} " as result %4.0f `ord_dpi'[1,2]  ///
		                                           _col(32) in gr "{c |}" as result %5.0f `df_dpi' ///
			                                       _col(40) in gr "{c |}  " as result %7.3f `imse_bsq_dpi' ///
												   _col(56) in gr "{c |}  " as result %7.3f `imse_var_dpi'
		di in smcl in gr "{rcenter 13: DPI-UKNOT}" _col(13) " {c |} " as result %4.0f `ord_dpi_uniq'[1,1] ///
		                                           _col(23) in gr " {c |} " as result %4.0f `ord_dpi_uniq'[1,2]  ///
                                                   _col(32) in gr "{c |}" as result %5.0f `df_dpi_uniq' ///
			                                       _col(40) in gr "{c |}  " as result %7.3f . ///
												   _col(56) in gr "{c |}  " as result %7.3f .
		di in smcl in gr "{hline 14}{c BT}{hline 8}{c BT}{hline 7}{c +}{hline 7}{c BT}{hline 15}{c BT}{hline 14}"
		di in smcl in gr "p: degree of polynomial. s: # of smoothness constraints. df: degrees of freedom."
	 }
	 else {
	 local imse_bsq_rot=`mat_imse_bsq_rot'[1,1]
	 local imse_var_rot=`mat_imse_var_rot'[1,1]
	 local imse_bsq_dpi=`mat_imse_bsq_dpi'[1,1]
	 local imse_var_dpi=`mat_imse_var_dpi'[1,1]
	 di in smcl in gr "{hline 28}{c +}{hline 10}"
	 di in smcl in gr "{ralign 27:Degree of polynomial}"  _col(28) " {c |} " _col(30) as result %7.0f `m_deg'[1,1]
	 di in smcl in gr "{ralign 27:# of smoothness constraint}"  _col(28) " {c |} " _col(30) as result %7.0f `m_deg'[1,2] 
     di in smcl in gr "{hline 28}{c BT}{hline 10}"
	 di ""
	 di in smcl in gr "{hline 14}{c TT}{hline 12}{c TT}{hline 10}{c TT}{hline 14}{c TT}{hline 14}"
	 di in smcl in gr "{rcenter 13: method}" _col(13) " {c |} " "{center 11: # of bins}" _col(26) "{c |}" "{rcenter 10: df}" _col(39) "{c |}" "{center 14: imse, bias^2}" _col(54) "{c |}" "{center 14: imse, var.}"
	 di in smcl in gr "{hline 14}{c +}{hline 12}{c +}{hline 10}{c +}{hline 14}{c +}{hline 14}"
	 di in smcl in gr "{rcenter 13: ROT-POLY}"  _col(13) " {c |} " as result %7.0f `J_rot_unreg'        _col(28) in gr "{c |}" as result %7.0f `df_rot_unreg' ///
	                                            _col(39) in gr "{c |}  " as result %7.3f `imse_bsq_rot' _col(54) in gr "{c |}  " as result %7.3f `imse_var_rot'
	 di in smcl in gr "{rcenter 13: ROT-REGUL}" _col(13) " {c |} " as result %7.0f `J_rot_reg'          _col(28) in gr "{c |}" as result %7.0f `df_rot_reg' ///
	                                            _col(39) in gr "{c |}  " as result %7.3f .              _col(54) in gr "{c |}  " as result %7.3f .
	 di in smcl in gr "{rcenter 13: ROT-UKNOT}" _col(13) " {c |} " as result %7.0f `J_rot_uniq'         _col(28) in gr "{c |}" as result %7.0f `df_rot_uniq' ///
	                                            _col(39) in gr "{c |}  " as result %7.3f .              _col(54) in gr "{c |}  " as result %7.3f .
	 di in smcl in gr "{rcenter 13: DPI}"       _col(13) " {c |} " as result %7.0f `J_dpi'              _col(28) in gr "{c |}" as result %7.0f `df_dpi' ///
	                                            _col(39) in gr "{c |}  " as result %7.3f `imse_bsq_dpi' _col(54) in gr "{c |}  " as result %7.3f `imse_var_dpi'
	 di in smcl in gr "{rcenter 13: DPI-UKNOT}" _col(13) " {c |} " as result %7.0f `J_dpi_uniq'         _col(28) in gr "{c |}" as result %7.0f `df_dpi_uniq' ///
	                                            _col(39) in gr "{c |}  " as result %7.3f .              _col(54) in gr "{c |}  " as result %7.3f .
	 di in smcl in gr "{hline 14}{c BT}{hline 12}{c BT}{hline 10}{c BT}{hline 14}{c BT}{hline 14}"
	 di in smcl in gr "df: degrees of freedom."
	 }
	 * return
	 * notes: J_rot_uniq is obtained possibly based on the subsample; J_dpi_uniq is ALWAYS obtained based on the full sample
	 ereturn clear
	 ereturn scalar N=`N'
	 ereturn scalar Ndist=`Ndist'
	 ereturn scalar Nclust=`Nclust'
	 ereturn scalar deriv=`deriv'
	 ereturn scalar imse_bsq_rot=`imse_bsq_rot'
	 ereturn scalar imse_var_rot=`imse_var_rot'
	 ereturn scalar imse_bsq_dpi=`imse_bsq_dpi'
	 ereturn scalar imse_var_dpi=`imse_var_dpi'
	 
	 ereturn scalar nbinsrot_poly=`J_rot_unreg'
	 ereturn scalar nbinsrot_regul=`J_rot_reg'
	 ereturn scalar nbinsrot_uknot=`J_rot_uniq'
	 ereturn scalar nbinsdpi=`J_dpi'
	 ereturn scalar nbinsdpi_uknot=`J_dpi_uniq'
	 
	 ereturn scalar prot_poly=`ord_rot_unreg'[1,1]
	 ereturn scalar prot_regul=`ord_rot_reg'[1,1]
	 ereturn scalar prot_uknot=`ord_rot_uniq'[1,1]
	 ereturn scalar pdpi=`ord_dpi'[1,1]
	 ereturn scalar pdpi_uknot=`ord_dpi_uniq'[1,1]
	 
	 ereturn scalar srot_poly=`ord_rot_unreg'[1,2]
	 ereturn scalar srot_regul=`ord_rot_reg'[1,2]
	 ereturn scalar srot_uknot=`ord_rot_uniq'[1,2]
	 ereturn scalar sdpi=`ord_dpi'[1,2]
	 ereturn scalar sdpi_uknot=`ord_dpi_uniq'[1,2]

	 ereturn matrix knot=`xkmat'
	 tempname m_p m_s
	 mat `m_p'=`m_deg'[1..`ncom',1]
	 mat `m_s'=`m_deg'[1..`ncom',2]
	 ereturn matrix m_p=`m_p'
	 ereturn matrix m_s=`m_s'
	 ereturn matrix m_nbinsrot_poly=`mat_J_rot_unreg'
	 ereturn matrix m_nbinsrot_regul=`mat_J_rot_reg'
	 ereturn matrix m_nbinsrot_uknot=`mat_J_rot_uniq'
	 ereturn matrix	m_nbinsdpi=`mat_J_dpi'
	 ereturn matrix m_nbinsdpi_uknot=`mat_J_dpi_uniq'
	
	 ereturn matrix m_imse_bsq_dpi=`mat_imse_bsq_dpi'
	 ereturn matrix m_imse_var_dpi=`mat_imse_var_dpi'
	 ereturn matrix m_imse_bsq_rot=`mat_imse_bsq_rot'
	 ereturn matrix m_imse_var_rot=`mat_imse_var_rot'
	 
	 

end

* Helper command
program define bins_imse, eclass

   version 13
   syntax varlist(min=2 numeric ts fv) [if] [in] [fw aw pw] [, deriv(integer 0) ///
          p(integer 0) s(integer 0) nbins(integer 0) eN_sub(integer 0) ///
		  binspos(string) vce(passthru) usegtools ///
		  zvec(name) absorb(string asis) reghdfeopt(string asis)]
   
   preserve
   marksample touse
   qui keep if `touse'
	 
   if ("`weight'"!="") local wt [`weight'`exp']
   
   tokenize `varlist'
   local y_var `1'
   local z_var `2'
   macro shift 2
   local w_var "`*'"
   
   tempvar zcat
   qui gen `zcat'=. in 1
   * Prepare bins
   tempname kmat
	
   if "`binspos'"=="ES" {
	  local stepsize=1/`nbins'
	  forvalues i=1/`=`nbins'+1' {
	 	 mat `kmat'=(nullmat(`kmat') \ `=0+`stepsize'*(`i'-1)')
	  }
   }
   else {
	  if (`nbins'==1) {
		 mat `kmat'=(0 \ 1)
	  }
	  else {
		 binsreg_pctile `z_var' `wt', nq(`nbins') `usegtools'
		 mat `kmat'=(0 \ r(Q) \ 1)
	  }
   }
   
   binsreg_irecode `z_var', knotmat(`kmat') bin(`zcat') ///
				   `usegtools' nbins(`nbins') pos(`binspos') knotliston(T)

   
	* Start computation
	tempvar derivfit derivse biasterm biasterm_v projbias
	qui gen `derivfit'=. in 1
	qui gen `derivse'=. in 1
	qui gen `biasterm'=. in 1                     /* save bias */
	if (`deriv'>0) qui gen `biasterm_v'=. in 1    /* error of approx deriv */
	qui gen `projbias'=. in 1                     /* save proj of bias */

	**************************************
	* predict leading bias
	mata: bias("`z_var'", "`zcat'", "`kmat'", `p', 0, "`biasterm'")
	if (`deriv'>0) {
		mata: bias("`z_var'", "`zcat'", "`kmat'", `p', `deriv', "`biasterm_v'")
	}

	* Increase order from p to p+1
	* Expand basis
	local nseries=(`p'-`s'+1)*(`nbins'-1)+`p'+2
	local series ""
	forvalues i=1/`nseries' {
		tempvar sp`i'
		local series `series' `sp`i''
		qui gen `sp`i''=. in 1
	}

	mata: binsreg_st_spdes(`zvec', "`series'", "`kmat'", st_data(.,"`zcat'"), `=`p'+1', 0, `=`s'+1')

	if ("`absorb'"=="") capture reg `y_var' `series' `w_var' `wt', nocon
	else                capture reghdfe `y_var' `series' `w_var' `wt', absorb(`absorb') `reghdfeopt'

	* store results
	tempname temp_b temp_V
	if (_rc==0) {
		matrix `temp_b'=e(b)
		matrix `temp_V'=e(V)
	}
	else {
		error  _rc
		exit _rc
	}
     
	tempname Xm
	* Predict (p+1)th derivative
	mata: `Xm'=binsreg_spdes(`zvec', "`kmat'", st_data(.,"`zcat'"), `=`p'+1', `=`p'+1', `=`s'+1'); ///
		  st_store(.,"`derivfit'", (binsreg_pred(`Xm', (st_matrix("`temp_b'")[|1 \ `nseries'|])', ///
		           st_matrix("`temp_V'")[|1,1 \ `nseries',`nseries'|], "xb"))[,1])
	mata: mata drop `Xm'

	qui replace `biasterm'=`derivfit'*`biasterm'
	if (`deriv'>0) qui replace `biasterm_v'=`derivfit'*`biasterm_v'
	drop `series'

	* Then get back degree-p spline, run OLS
	local nseries=(`p'-`s'+1)*(`nbins'-1)+`p'+1
	local series ""
	forvalues i=1/`nseries' {
		tempvar sp`i'
		local series `series' `sp`i''
		qui gen `sp`i''=. in 1
	}

	mata: binsreg_st_spdes(`zvec', "`series'", "`kmat'", st_data(.,"`zcat'"), `p', 0, `s')   
	capture reg `biasterm' `series' `wt', nocon       /* project bias on X of degree p */
	tempname bias_b bias_V
	if (_rc==0) {
		matrix `bias_b'=e(b)
		matrix `bias_V'=e(V)
	}
	else {
		error _rc
		exit _rc
	}

	mata: `Xm'=binsreg_spdes(`zvec', "`kmat'", st_data(.,"`zcat'"), `p', `deriv', `s'); ///
		  st_store(.,"`projbias'", binsreg_pred(`Xm', st_matrix("`bias_b'")', st_matrix("`bias_V'"), "xb")[,1])

	if (`deriv'==0) {
		qui replace `biasterm'=(`biasterm'-`projbias')^2
	}
	else {
		qui replace `biasterm'=(`biasterm_v'-`projbias')^2   /* still save in biasterm if deriv>0 */
	}

	if ("`wt'"!="") qui sum `biasterm' [aw`exp'], meanonly
	else            qui sum `biasterm', meanonly
	local m_bias=r(mean)
	local imse_b=`m_bias'*`nbins'^(2*(`p'+1-`deriv'))

	* for variance purpose
	if ("`absorb'"=="") capture reg `y_var' `series' `w_var' `wt', nocon `vce'        
	else                capture reghdfe `y_var' `series' `w_var' `wt', absorb(`absorb') `vce' `reghdfeopt'

	* store results
	if (_rc==0) {
		matrix `temp_b'=e(b)
		matrix `temp_V'=e(V)
		tempname vcov
		mata: `vcov'=st_matrix("`temp_V'")
		if ("`absorb'"=="") mata: `vcov'=`vcov'[|1,1 \ `nseries',`nseries'|]
		else {
			mata: `vcov'=(`vcov'[|1,1 \ `nseries', `nseries'|], `vcov'[|1,cols(`vcov') \ `nseries', cols(`vcov')|] \ ///
				  `vcov'[|cols(`vcov'), 1 \ cols(`vcov'), `nseries'|], `vcov'[cols(`vcov'), cols(`vcov')]); ///
				  `Xm'=(`Xm', J(rows(`Xm'),1,1))
		}
	}
	else {
		error  _rc
		exit _rc
	}

	mata: st_store(., "`derivse'", (binsreg_pred(`Xm', ., `vcov', "se")[,2]):^2)
	mata: mata drop `vcov'

	if ("`wt'"!="") qui sum `derivse' [aw`exp'], meanonly
	else            qui sum `derivse', meanonly
	local m_se=r(mean)
	local imse_v=`m_se'/(`nbins'^(1+2*`deriv'))	   
	mata: mata drop `Xm'   	   

	* DPI J
	local J_dpi=ceil((`imse_b'*2*(`p'+1-`deriv')/               ///
		(`imse_v'*(1+2*`deriv')))^(1/(2*`p'+2+1)))

	local imse_bsq_dpi=`imse_b'
	local imse_var_dpi=`imse_v'*`eN_sub'
 
    ereturn clear
	ereturn scalar imse_bsq=`imse_bsq_dpi'
	ereturn scalar imse_var=`imse_var_dpi'
	
end


version 13
mata:
    // Constant in variance
    void imse_v_cons(real scalar degree, real scalar deriv, string scalar vcons)
	{  
	   real scalar v_cons, m
	   real matrix V, Vderiv
	   
	   m=degree+1
       if (deriv==0) {
	      v_cons=m
	   }
	   else {
	      V=J(m, m, .)
		  Vderiv=J(m, m, 0)
		  
		  for (i=1; i<=m; i++){
		      for (j=1; j<=i; j++) {
			       V[i,j]=1/(i+j-1)
			       if (i>deriv & j>deriv) {
				       Vderiv[i,j]=1/(i+j-1-2*deriv)* /*
					   */          (factorial(i-1)/factorial(i-1-deriv))* /*
					   */          (factorial(j-1)/factorial(j-1-deriv))
				   }
			  }
		  }
		  V=makesymmetric(V)
		  Vderiv=makesymmetric(Vderiv)
		  v_cons=trace(invsym(V)*Vderiv)
	   }
	   
	   // return results
	   st_numscalar(vcons,v_cons)
	}
	
	// Constant in bias
	void imse_b_cons(real scalar degree, real scalar deriv, string scalar bcons, | real scalar s)
	{  
	   real scalar b_cons, m, bernum
	   m=degree+1
       if (args()<4) {
	      b_cons=1/(2*(m-deriv)+1)/factorial(m-deriv)^2/comb(2*(m-deriv), m-deriv)^2
	   }
	   else {
	      if (degree==0) {
		     bernum=1/6
		  }
		  else if (degree==1) {
		     bernum=1/30
		  }
		  else if (degree==2) {
		     bernum=1/42
		  }
		  else  if (degree==3) {
		     bernum=1/30
		  }
		  else if (degree==4) {
		     bernum=5/66
		  }
		  else if (degree==5) {
		     bernum=691/2730
		  }
		  else if (degree==6) {
		     bernum=7/6
		  }
		  else {
		     _error("p>6 not allowed.")
		  }
	   	  b_cons=1/factorial(2*(m-deriv))*bernum
	   }
	   
	   // return results
	   st_numscalar(bcons, b_cons)
	}
  
    // Bernoulli polynomial
    real vector bernpoly(real vector x, real scalar degree)
    {
      n=rows(x)
	  if (degree==0) {
	     bernx=J(n,1,1)
	  }
	  else if (degree==1) {
	     bernx=x:-0.5
	  }
	  else if (degree==2) {
	     bernx=x:^2-x:+1/6
	  }
	  else if (degree==3) {
	     bernx=x:^3-1.5*x:^2+0.5*x
	  }
	  else if (degree==4) {
	     bernx=x:^4-2*x:^3+x:^2:-1/30
	  }
	  else if (degree==5) {
	     bernx=x:^5-2.5*x:^4+5/3*x:^3-1/6*x
	  }
	  else if (degree==6) {
	     bernx=x:^6-3*x:^5+2.5*x:^4-0.5*x:^2:+1/42
	  }
	  else {
	     _error("p is too large.")
	  }
	  return(bernx)
  }

  // Leading bias for splines
  void bias(string scalar Var, string scalar Xcat, string scalar knotname, ///
            real scalar degree, real scalar deriv, ///
			string scalar biasname, | string scalar select)
  {
    if (args()<7) {
	   X=st_data(., (Var))
	   xcat=st_data(., (Xcat))
	   st_view(bias=., ., (biasname))
	}
	else {
       X=st_data(., (Var), select)
	   xcat=st_data(., (Xcat), select)
	   st_view(bias=.,.,(biasname), select)
	}
	knot=st_matrix(knotname)
 	h=knot[|2 \ length(knot)|]-knot[|1 \ (length(knot)-1)|]
	h=h[xcat]
	if (rows(h)==1) {
	   h=h'
	}
	tl=knot[|1 \ (length(knot)-1)|]
	tl=tl[xcat]
	if (rows(tl)==1) {
	   tl=tl'
	}
	bern=bernpoly((X-tl):/h, degree+1-deriv)/factorial(degree+1-deriv):*(h:^(degree+1-deriv))
	bias[.,.]=bern
  }
  
  // find the minimum
  void findmindex(string scalar matname, string scalar outname, ///
                  real scalar J, real scalar nr)
  {
     real matrix A
	 
     A=sort((abs(st_matrix(matname):-J), (1::nr)), 1)
	 st_numscalar(outname, A[1,2])
  }
  
end


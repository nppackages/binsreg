*! version 1.3 03-Jul-2023 
 
* Generalized pctile function

program define binsreg_pctile, rclass
   version 13
   
   syntax varlist(max=1 numeric) [if] [in] [fw aw pw] [, nq(integer 2) usegtools]
   
   /* nq must be >=2, # of bins */
   /* used internally, no error checks */
   
   marksample touse
   
   if ("`weight'"!="") local wt [`weight'`exp']
   
   local nk = `nq'-1   /* # of quantiles */
   tempname A
   
   mat `A'=J(`nk',1,.)
	  
   if `nq' <= 1001 {
	  if ("`usegtools'"=="") _pctile `varlist' if `touse' `wt', nq(`nq')
	  else                   gquantiles `varlist' if `touse' `wt', _pctile nq(`nq')
	  forvalues j = 1/`nk' {
		 mat `A'[`j',1] = r(r`j')
	  }
   }
   else {
	  local plist    = ""
      local plistlen = 0
	  local entryint = 1
	  local plistsize = 0
	  local plistnextsize = 0

	  forvalues j = 1/`nk' {
		 local plistnext = 100*`j'/`nq'
		 local plistsize = `plistsize' + `plistnextsize'
		 local plistnextsize : strlen local plistnext
		 if (`plistsize'+`plistnextsize' < `c(macrolen)') & (`plistlen' <= 999) {
			if "`plist'"=="" {
			   local plist "`plistnext'"
			}
			else local plist "`plist',`plistnext'"

			local ++plistlen

			if `j' == `nk' {
				if ("`usegtools'"=="") _pctile `varlist' if `touse' `wt', p(`plist')
				else                   gquantiles `varlist' if `touse' `wt', _pctile p(`plist')
				forvalues k = 1/`plistlen' {
					mat `A'[`entryint'+`k'-1,1] = r(r`k')
				}
			}
		 }
		 else {
			if ("`usegtools'"=="") _pctile `varlist' if `touse' `wt', p(`plist')
			else                   gquantiles `varlist' if `touse' `wt', _pctile p(`plist')
			forvalues k = 1/`plistlen' {
				mat `A'[`entryint'+`k'-1,1] = r(r`k')
			}
			if (`j'<`nk') {
				local plist "`plistnext'"
				local entryint = `entryint' + `plistlen'
				local plistlen = 1
				local plistsize = 0
			}
			else {
				if ("`usegtools'"=="") _pctile `varlist' if `touse' `wt', p(`plistnext')
				else                   gquantiles `varlist' if `touse' `wt', _pctile p(`plistnext')
				mat `A'[`nk',1]=r(r1)
			}
		 }
	  }
   }
  
  return clear
  return matrix Q=`A'
   
end

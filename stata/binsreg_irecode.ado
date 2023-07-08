*! version 1.3 03-Jul-2023 
 
* Generalized irecode function

program define binsreg_irecode
   version 13
   
   syntax varlist(max=1 numeric) [if] [in], knotmat(name) BINid(varname numeric) ///
                                            [usegtools nbins(string) pos(string) knotliston(string)] 
   
   /* knot is a FULL knot matrix with boundaries */
   /* used internally, no error checks */
   
   marksample touse
   
   if ("`usegtools'"==""|"`pos'"!="QS") {
      confirm variable `binid'
	  confirm matrix `knotmat'
   
      local n=rowsof(`knotmat')
      if (`n'==2) qui replace `binid'=1 if `touse'
      else if (`n'==3) qui replace `binid'=1+irecode(`varlist', `knotmat'[2,1]) if `touse'
      else {
         local J = `n' - 2
         local knots = `knotmat'[2,1]
         local assignedbins = 1
	     local knotstsize : strlen local knots
	     local knotsnextsize = 0

         forvalues j = 2/`J' {
	        local knotsnext = `knotmat'[`j'+1,1]
	        local knotstsize = `knotstsize' + `knotsnextsize'
	        local knotsnextsize : strlen local knotsnext

	        if (`knotstsize' + `knotsnextsize' < `c(macrolen)') & (`j' - `assignedbins' < 248){
		       local knots "`knots',`knotsnext'"
		       if `j' == `J' {
			      qui replace `binid' = `assignedbins' + irecode(`varlist',`knots') if `touse' & `binid'==.
		       }
	        }
	        else {
		       qui replace `binid' = `assignedbins' + irecode(`varlist',`knots') if `binid'==. & `touse'
		       if `j'<`J' {
			     qui replace `binid' = . if `varlist' > `knotmat'[`j'+1,1] & `touse'
			     local knots = `knotmat'[`j'+1,1]
			     local assignedbins  = `j'
			     local knotstsize = 0
		       }
			   else {
			     qui replace `binid' = `n'-1 if `varlist' > `knotmat'[`j'+1,1] & `touse'
			   }
	        }
         } 
      }
   }
   else {
      capture confirm variable `binid'
	  if (!_rc) drop `binid'
	  
      if (`nbins'==1) qui replace `binid'=1 if `touse' 
      else {
	     tempvar cat
         if ("`knotliston'"=="T") fasterxtile `binid'=`varlist', nq(`nbins')
		 else                     fasterxtile `binid'=`varlist'  if `touse', nq(`nbins')
		 *qui replace `binid'=`cat' if `touse'
	  }
   }
   
   
end

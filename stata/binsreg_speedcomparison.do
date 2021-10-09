clear all 
set obs 10000000

local rep=10
gen x=runiform()
gen y=x
mat time=J(10,`rep',.)

forval i=1/`rep' {
   timer on 1
   binscatter y x
   timer off 1

   timer on 2
   binscatter2 y x
   timer off 2
   

   timer on 3
   binsreg y x, polyreg(1) nbins(20)
   timer off 3
   
   timer on 4
   binsreg y x, polyreg(1) nbins(20) masspoints(off)
   timer off 4

   timer on 5
   binsreg y x, polyreg(1) nbins(20) masspoints(off) usegtools(on)
   timer off 5

   
   preserve
   sort x

   timer on 6
   binscatter y x
   timer off 6

   timer on 7
   binscatter2 y x
   timer off 7

   timer on 8
   binsreg y x, polyreg(1) nbins(20)
   timer off 8
   
   timer on 9
   binsreg y x, polyreg(1) nbins(20) masspoints(off)
   timer off 9

   timer on 10
   binsreg y x, polyreg(1) nbins(20) masspoints(off) usegtools(on)
   timer off 10

   restore

   timer list
   forval j=1/10 {
      mat time[`j',`i']=r(t`j')
   }
   
   timer clear
}

mat ave_time=time*J(`rep',1,1)/`rep'

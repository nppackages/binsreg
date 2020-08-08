clear all
set obs 10000000

gen x=runiform()
gen y=x
mat time=J(8,1,.)

forval i=1/1 {
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
   
   preserve
   sort x

   timer on 5
   binscatter y x
   timer off 5

   timer on 6
   binscatter2 y x
   timer off 6

   timer on 7
   binsreg y x, polyreg(1) nbins(20)
   timer off 7
   
   timer on 8
   binsreg y x, polyreg(1) nbins(20) masspoints(off)
   timer off 8

   restore

   timer list
   forval j=1/8 {
      mat time[`j',`i']=r(t`j')
   }
   
   timer clear
}

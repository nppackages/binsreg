* version 1.3, 03-Jul-2023
*****************************************************************
****** This file contains necessary mata functions used in ******
******************** BINSREG Package ****************************
version 13
 mata:
 
  // Generate spline design matrix using MATA objects
  real matrix binsreg_spdes(real vector xvar, ///
                            string scalar knotname, real vector xcatvar, ///
							real scalar degree, real scalar deriv, real scalar s, ///
							| real vector by, real scalar byvalue)  
				  /* s: # of constraints */   
  { 
	real matrix knot, x, xcat, subgroup, exknot, ind_lk, lk, rk, bs, tl, tr, w, vm, rind
	real scalar k, n, i, sm, cind

	// read data into mata
	knot=st_matrix(knotname)
    k=rows(knot)
	if (args()<7) {
	   x=xvar
	   xcat=xcatvar
	}
	else {
	   subgroup = by:==byvalue
	   x=select(xvar, subgroup)
	   xcat=select(xcatvar, subgroup)
	}
	
	n=rows(x)
	sm=s-1                /* so that resulting spline is C^sm */
    
	// extended knots
	if (k>2) {
        exknot=(J(degree+1,1,knot[1])\ ///
		        knot[|2\(k-1)|]#J(degree-sm,1,1)\ /*
		*/      J(degree+1,1,knot[k]))
    }
    else {
     	exknot=(J(degree+1,1,knot[1])\ /*
		*/      J(degree+1,1,knot[k]))
    }
	
	// left and right knot subscripts for each x
	ind_lk=(degree+1):+(xcat:-1)*(degree-sm)
	lk=J(n,degree,.)
	rk=J(n,degree,.)
	
	// knot matrix: left and right endpoints of local supports
	for (i=1; i<=degree; ++i) {
	    lk[.,(degree-i+1)]=exknot[ind_lk:+(-i+1)]
		rk[.,i]=exknot[ind_lk:+i]
	}
	
	// initial nonzero b spline vector
	bs=J(n,1,1)
	
	// loop based on recursive formula
	if (degree >= 1) {
	   if (degree < deriv) {
	       bs=J(n,degree+1,0)
	   }
	   else if (degree > deriv) {	   
	       // loop: up to the degree of "degree-deriv"
		   for (i=1; i<=degree-deriv; i++) {
	           tl=lk[|1,(degree-i+1) \ n,degree|]
		       tr=rk[|1,1 \ n,i|]
		       w=(x:-tl):/(tr-tl)
		       bs=((1:-w):*bs, J(n,1,0))+(J(n,1,0), w:*bs)
	       }
		   // loop: derivative
		   if (deriv > 0) {
	          for (i=degree-deriv+1; i<=degree; i++) {
	              tl=lk[|1,(degree-i+1) \ n,degree|]
		          tr=rk[|1,1 \ n,i|]
		          w=1:/(tr-tl)
		          bs=((J(n,1,0), w:*bs)-(w:*bs, J(n,1,0)))*i
	          }
		   }
        }
		else {
		   // degree=deriv, so only need to account for derivative
	       for (i=1; i<=degree; i++) {
	           tl=lk[|1,(degree-i+1) \ n,degree|]
		       tr=rk[|1,1 \ n,i|]
		       w=1:/(tr-tl)
		       bs=((J(n,1,0), w:*bs)-(w:*bs, J(n,1,0)))*i
	       }
		}
	}
	
	// bs should be an n by degree+1 matrix containing nonzeros of
	// design matrix for each x
	
	// expand it to a full matrix: bin by bin
	// zero all
	vm=J(n,(k-2)*(degree-sm)+degree+1,0)
	
	
	// fill in nonzeros
	for (i=1; i<=k-1; i++) {
	    rind=selectindex(xcat:==i)		
	    cind=(i-1)*(degree-sm)+1
	    vm[rind,cind..(cind+degree)]=bs[rind,]
	}
	
	return(vm)
  }
  
  mata mosave binsreg_spdes(), replace

  
  // Gen spline design (transferred to STATA)
  void binsreg_st_spdes(real vector xvar, string scalar spvar, ///
                        string scalar knotname, real vector catname, ///
						real scalar degree, real scalar deriv, real scalar s, ///
						| string scalar select)  
			  /* s: # of constraints */   
  {  
	real matrix X, vm, xcat
	
	xcat = catname
	X=binsreg_spdes(xvar, knotname, xcat, degree, deriv, s)
	if (args()<8) {
	   st_view(vm=., ., (spvar))
	}
	else {
	   st_view(vm=., ., (spvar), select)
	}
	
	// put it back to stata
	vm[.,.]=X
  }
   mata mosave binsreg_st_spdes(), replace

  

  // Generate grids for plotting and simulation
  real matrix binsreg_grids(string scalar knotname, real scalar n)
  {
    real vector knot, isknot, bin, h, grid, V
	real scalar r
	
    knot=st_matrix(knotname)
	r=rows(knot)
	
	isknot=J(r-1,1,1)#(1 \ J(n,1,0))   /* indicator for knots */
	isknot=isknot[|2 \ rows(isknot)|]
	
	bin=range(1, r-1, 1)#J(n+1,1,1)
	bin=bin[|1 \ (rows(bin)-1)|]
	
	h=(knot[|2\r|]-knot[|1 \ (r-1)|])/(n+1)
	grid=(knot[1] \ h#J(n+1,1,1))                /* grids containing knots */
	
	// remove the first and last
	grid=runningsum(grid)
	grid=grid[|2\(rows(grid)-1)|]
	
	V=(grid, isknot, bin)
	return(V)
  }
  
  mata mosave binsreg_grids(), replace

 
  // General prediction function
  real matrix binsreg_pred(real matrix X, real matrix beta, real matrix cov, ///
                           string scalar type)
  {
    real matrix est, se
	
	est=J(rows(X), 1, .)
	se=J(rows(X), 1, .)
    if (type=="all"|type=="xb") {
	    est=X*beta
	}
	if (type=="all"|type=="se") {
	    se=sqrt(rowsum((X*cov):*X))
	}
	return((est, se))
  }
  mata mosave binsreg_pred(), replace

  
  // Simulation-based pval and critical value
  void binsreg_pval(real matrix X, real vector se, ///
                    string scalar covname, string scalar mtest, ///
		            real scalar rep, real scalar k, ///
		            string scalar side, real scalar alpha, ///
		            string scalar pmat, string scalar cval, string scalar metric)
  {
     real matrix cov, mt, pval, tvec, num, U, V, sv, t, p
	 real scalar i,j,w,j1, lp
	 
	 cov=st_matrix(covname)[|1,1\k,k|]
	 if (metric!="inf") {
	    lp=strtoreal(metric)
	 }
	 if (mtest!=".") {
	     mt=st_matrix(mtest)
	     pval=J(rows(mt),1,0)
	 }
	 if (side!=".") tvec=J(rep,1,.)
	 if (rank(cov)==k) {
	     num=X*cholesky(cov)
	 }
	 else {
	     svd(cov, U=., sv=., V=.)
		 pragma unused V
		 num=X*U*diag(sv:^0.5)*U'
	 }
	  
	 for (i=1; i<=rep; i++) {
	     t=(num*rnormal(k,1,0,1)):/se
		 // for p-vals
		 if (mtest!=".") {
		    for (j=1;j<=rows(mt); j++) {
		        // 1: left; 2: right; 3: two-sided
		        if (mt[j,2]==1) {
			        pval[j]=pval[j]+(max(t)>=mt[j,1])
			    }
			    else if (mt[j,2]==2) {
   			        pval[j]=pval[j]+(min(t)<=mt[j,1])
			    }
			    else {
				    if (metric=="inf") {
			           pval[j]=pval[j]+(max(abs(t))>=mt[j,1])
					}
					else {
					   pval[j]=pval[j]+(mean(abs(t):^lp)^(1/lp)>=mt[j,1])
					}
			    }
		    }
		 }
		 // for critical value
		 if (side=="two") {
		     if (metric=="inf") {
		        tvec[i]=max(abs(t))
			 }
			 else {
			    tvec[i]=mean(abs(t):^lp)^(1/lp)
			 }
	     }
		 else if (side=="left") {
			 tvec[i]=max(t)
	     }
	     else if (side=="right") {
			 tvec[i]=min(t)
	     }
	 }
	  
	 if (mtest!=".") st_matrix(pmat,pval/rep)
	 if (side!=".") {
	     p=order(tvec,1)
	     w = alpha*rep
         j = floor(w)
         w = 0.5 + 0.5*((w - j)>0)
	     j1=j+1
	     if (j<1) j=1
	     if (j1>rep) j1=rep
	     st_numscalar(cval, (1-w)*tvec[p[j]] + w*tvec[p[j1]])
	 }
   }
   mata mosave binsreg_pval(), replace


   // Check dropping vars
   void binsreg_checkdrop(string scalar betaname, string scalar covname, real scalar k, | string scalar useqreg)
   {
     real matrix beta, se
	 real scalar isdrop
	 
     beta=st_matrix(betaname)
	 beta=beta[|1\k|]'
     se=diagonal(st_matrix(covname))[|1\k|]
     isdrop=sum(((beta:==0)+(se:==0)):==2)
	 if (args()==3) {
		pragma unused useqreg
		if (isdrop>0) {
	       display("{gr:Warning: some {it:basis functions} are dropped.}")
        }
     }
	 else {
	    if (isdrop>1) {
	       display("{gr:Warning: some {it:basis functions} are dropped.}")
        }
	 }
   }
   mata mosave binsreg_checkdrop(), replace

   // Check eff. N by bin
   real matrix binsreg_uniq(real vector x, real vector xcat, real scalar nbins, ///
                     string scalar uniqmin, |real vector by, real scalar byval)
   {
     real matrix subgroup, xsub, binid, binedges, unique 
	 real scalar j
	 
	 xsub=x
	 if (args()>4) {
	    subgroup = by:==byval
		xsub=select(xsub, subgroup)
	 }
	
	 if (nbins==1) {
	 	binedges = 0 \ rows(xsub)
	 } 
	 else {
	    binid=xcat
		if (args()>4) {
		   binid = select(binid, subgroup)
		}
		binedges = 0\selectindex(binid[|1 \ rows(binid)-1|]-binid[|2 \ rows(binid)|])\rows(xsub)
	 }
	 
	 unique=J(nbins,1,.)
	
	 for (j=1;j<=nbins;j++){
	    if (binedges[j+1]-binedges[j]==1) {
		   unique[j] = 1
		}
		else {
	 	   unique[j] = sum(xsub[|binedges[j]+1 \ binedges[j+1]-1|]:!=xsub[|binedges[j]+2 \ binedges[j+1]|]) + 1
		}
	 }
	 st_local(uniqmin, strofreal(min(unique)))
	 
	 return(binedges)
   }
   
   mata mosave binsreg_uniq(), replace
   
   
   // Stats by bin, with sorted data
   real matrix binsreg_stat(real matrix x, real vector xcat, real scalar nbins, ///
                            real vector edge, string scalar stat, real scalar quantile, ///
	 				        | real vector by, real scalar byval)
   {
      real matrix subgroup, xsub, binid, binedges, out 
	  real scalar j, ncol
	 
	  xsub=x
	  if (args()>6) {
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
		    if (args()>6) {
		      binid = select(binid, subgroup)
		    }
		    binedges = 0\selectindex(binid[|1 \ rows(binid)-1|]-binid[|2 \ rows(binid)|])\rows(xsub)
	     }
	  }
	  else {
	     binedges = edge
	  }
	 
	  out=J(nbins,ncol,.)
	 
	  if (stat=="mean") {
	     for (j=2;j<=nbins+1;j++) {
	 	    out[j-1,.] = mean(xsub[|binedges[j-1]+1, 1 \ binedges[j], ncol|])
	     }
	  }
	  else {
	     for (j=2;j<=nbins+1;j++) {
	 	    out[j-1,.] = binsreg_cquantile(xsub[|binedges[j-1]+1, 1 \ binedges[j], ncol|], quantile)
	     }
	  }
	 
	  out=(range(1,nbins,1), out)
	 
	  return(out)
   }
   
   mata mosave binsreg_stat(), replace
 
   // column quantile based on sorted data
   real matrix binsreg_cquantile(real matrix X, real scalar quantile)
   {
      real matrix p, out
	  real scalar i, w, j, j1, nrow, ncol
       
	  nrow=rows(X)
	  ncol=cols(X)
	  w = quantile*nrow
      j = floor(w)
      w = 0.5 + 0.5*((w - j)>0)
	  j1=j+1
	  if (j<1) j=1
	  if (j1>nrow) j1=nrow
	  
	  out=J(1, ncol, .)
	  for (i=1;i<=ncol;i++) {
	     p=order(X[,i],1)
	     out[1,i]=(1-w)*X[p[j],i] + w*X[p[j1],i]
      }
	  
	  return(out)  
   }
   
   mata mosave binsreg_cquantile(), replace
  
end


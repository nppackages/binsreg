{smcl}
{* *! version 1.3 03-JUL-2023}{...}
{viewerjumpto "Syntax" "binsregselect##syntax"}{...}
{viewerjumpto "Description" "binsregselect##description"}{...}
{viewerjumpto "Options" "binsregselect##options"}{...}
{viewerjumpto "Examples" "binsregselect##examples"}{...}
{viewerjumpto "Stored results" "binsregselect##stored_results"}{...}
{viewerjumpto "References" "binsregselect##references"}{...}
{viewerjumpto "Authors" "binsregselect##authors"}{...}
{cmd:help binsregselect}
{hline}

{title:Title}

{p 4 8}{hi:binsregselect} {hline 2} Data-driven IMSE-Optimal Partitioning/Binning Selection for Binscatter.{p_end}


{marker syntax}{...}
{title:Syntax}

{p 4 18} {cmdab:binsregselect} {depvar} {it:indvar} [{it:othercovs}] {ifin} {weight} [{cmd:,} {opt deriv(v)}{p_end}
{p 18 18} {opt absorb(absvars)} {opt reghdfeopt(reghdfe_option)}{p_end}
{p 18 18} {opt bins(p s)} {opt binspos(position)} {opt binsmethod(method)} {opt nbinsrot(#)} {opt nbins(nbinsopt)}{p_end}
{p 18 18} {cmd:pselect(}{it:{help numlist}}{cmd:)} {cmd:sselect(}{it:{help numlist}}{cmd:)}{p_end}
{p 18 18} {opt simsgrid(#)} {opt savegrid(filename)} {opt replace}{p_end}
{p 18 18} {opt dfcheck(n1 n2)} {opt masspoints(masspointsoption)}{p_end}
{p 18 18} {cmd:vce(}{it:{help vcetype}}{cmd:)} {opt usegtools(on/off)} {opt useeffn(#)} {opt randcut(#)} ]{p_end}

{p 4 8} where {depvar} is the dependent variable, {it:indvar} is the independent variable for binning, and {it:othercovs} are other covariates to be controlled for.{p_end}

{p 4 8} The degree of the piecewise polynomial p, the number of smoothness constraints s, and the derivative order v are integers 
satisfying 0 <= s,v <= p, which can take different values in each case.{p_end}

{p 4 8} {opt fweight}s, {opt aweight}s and {opt pweight}s are allowed; see {help weight}.{p_end}

{marker description}{...}
{title:Description}

{p 4 8} {cmd:binsregselect} implements data-driven procedures for selecting the number of bins for binscatter estimation.
The selected number is optimal in minimizing the (asymptotic) integrated mean squared error (IMSE).
{p_end}


{marker options}{...}
{title:Options}

{dlgtab:Estimand}

{p 4 8} {opt deriv(v)} specifies the derivative order of the regression function for estimation, testing and plotting.
The default is {cmd:deriv(0)}, which corresponds to the function itself.{p_end}

{dlgtab:Reghdfe}

{p 4 8} {opt absorb(absvars)} specifies categorical variables (or interactions) representing the fixed effects to be absorbed.
This is equivalent to including an indicator/dummy variable for each category of each {it:absvar}. When {cmd:absorb()} is specified,
the community-contributed command {cmd:reghdfe} instead of the command {cmd:regress} is used.
{p_end}

{p 4 8} {opt reghdfeopt(reghdfe_option)} options to be passed on to {cmd:reghdfe}. Important: {cmd:absorb()} and {cmd:vce()} should not be specified within this option.
{p_end}

{p 4 8} For more information about the community-contributed command {cmd:reghdfe}, please see {browse "http://scorreia.com/software/reghdfe/":http://scorreia.com/software/reghdfe/}.
{p_end}

{dlgtab:Binning/Degree/Smoothness Selection}

{p 4 8} {opt bins(p s)} sets a piecewise polynomial of degree {it:p} with {it:s} smoothness constraints for
data-driven (IMSE-optimal) selection of the partitioning/binning scheme.
The default is {cmd:bins(0 0)}, which corresponds to piecewise constant (canonical binscatter).

{p 4 8} {opt binspos(position)} specifies the position of binning knots.
The default is {cmd:binspos(qs)}, which corresponds to quantile-spaced binning (canonical binscatter).
Other option is {cmd:es} for evenly-spaced binning.
{p_end}

{p 4 8} {opt binsmethod(method)} specifies the method for data-driven selection of the number of bins.
The default is {cmd:binsmethod(dpi)}, which corresponds to the IMSE-optimal direct plug-in rule.
The other option is: {cmd:rot} for rule of thumb implementation.
{p_end}

{p 4 8} {opt nbinsrot(#)} specifies an initial number of bins value used to construct the DPI number of bins selector.
If not specified, the data-driven ROT selector is used instead.
{p_end}

{p 4 8} {opt nbins(nbinsopt)} sets the number of bins for degree/smoothness selection. 
If {cmd:nbins(T)} is specified, the command selects the number of bins instead, 
given the specified degree and smoothness. 
If a {help numlist:numlist} with more than one number is specified, 
the command selects the number of bins within this list.
{p_end}

{p 4 8} {opt pselect(numlist)} specifies a list of numbers within which the degree of polynomial {it:p} for 
point estimation is selected.
{p_end}

{p 4 8} {opt sselect(numlist)} specifies a list of numbers within which the number of smoothness constraints {it:s}
for point estimation is selected. If not specified, for each value {it:p} supplied in the 
option {cmd:pselect()}, only the piecewise polynomial with the maximum smoothness is considered, i.e., {it:s=p}.  
{p_end}

{p 4 8} Note: To implement the degree or smoothness selection, in addition to {cmd:pselect()} 
or {cmd:sselect()}, {cmd:nbins(#)} must be specified.
{p_end}

{dlgtab:Evaluation Points Grid Generation}

{p 4 8} {opt simsgrid(#)} specifies the number of evaluation points of an evenly-spaced grid within each bin used
for evaluation of the supremum (infimum or Lp metric) operation needed to construct confidence bands and hypothesis testing procedures.
The default is {cmd:simsgrid(20)}, which corresponds to 20 evenly-spaced evaluation points within each bin for
approximating the supremum (or infimum) operator.
{p_end}

{p 4 8} {opt savegrid(filename)} specifies a filename for storing the simulation grid of evaluation points.
It contains the following variables:
{it:indvar}, which is a sequence of evaluation points used in approximation;
all control variables in {it:othercovs}, which take values of zero for prediction purpose;
{it:binsreg_isknot}, indicating  whether the evaluation point is an inner knot;
and {it:binsreg_bin}, indicating which bin the evaluation point belongs to.
{p_end}

{p 4 8} {opt replace} overwrites the existing file when saving the grid.
{p_end}

{dlgtab:Mass Points and Degrees of Freedom}

{p 4 8} {opt dfcheck(n1 n2)} sets cutoff values for minimum effective sample size checks,
which take into account the number of unique values of {it:indvar} (i.e., adjusting for the number of mass points),
number of clusters, and degrees of freedom of the different statistical models considered.
The default is {cmd:dfcheck(20 30)}. See Cattaneo, Crump, Farrell and Feng (2023c) for more details.
{p_end}

{p 4 8} {opt masspoints(masspointsoption)} specifies how mass points in {it:indvar} are handled.
By default, all mass point and degrees of freedom checks are implemented.
Available options:
{p_end}
{p 8 8} {opt masspoints(noadjust)} omits mass point checks and the corresponding effective sample size adjustments.{p_end}
{p 8 8} {opt masspoints(nolocalcheck)} omits within-bin mass point and degrees of freedom checks.{p_end}
{p 8 8} {opt masspoints(off)} sets {opt masspoints(noadjust)} and {opt masspoints(nolocalcheck)} simultaneously.{p_end}
{p 8 8} {opt masspoints(veryfew)} forces the command to proceed as if {it:indvar} has only a few number of mass points (i.e., distinct values).
In other words, forces the command to proceed as if the mass point and degrees of freedom checks were failed.{p_end}

{dlgtab:Other Options}

{p 4 8} {cmd:vce(}{it:{help vcetype}}{cmd:)} specifies the {it:vcetype} for variance estimation
used by the command {help regress##options:regress} (or {cmd:reghdfe} if {cmd:absorb()} is specified).
The default is {cmd:vce(robust)}.
{p_end}

{p 4 8}{opt usegtools(on/off)} forces the use of several commands in the community-distributed Stata package {cmd:gtools}
to speed the computation up, if {it:on} is specified.
Default is {cmd:usegtools(off)}.
{p_end}

{p 4 8} For more information about the package {cmd:gtools}, please see {browse "https://gtools.readthedocs.io/en/latest/index.html":https://gtools.readthedocs.io/en/latest/index.html}.
{p_end}

{p 4 8} {opt useeffn(#)} specifies the effective sample size {it:#} to be used when computing the (IMSE-optimal) number of bins.
This option is useful for extrapolating the optimal number of bins to larger (or smaller) datasets than the one used to compute it.
{p_end}

{p 4 8} {opt randcut(#)} specifies the upper bound on a uniformly distributed variable used to draw a subsample for bins selection.
Observations for which {cmd:runiform()<=#} are used. # must be between 0 and 1.
    
{marker examples}{...}
{title:Examples}

{p 4 8} Setup{p_end}
{p 8 8} . {stata sysuse auto}{p_end}

{p 4 8} Select IMSE-optimal number of bins using DPI-procedure{p_end}
{p 8 8} . {stata binsregselect mpg weight foreign}{p_end}


{marker stored_results}{...}
{title:Stored results}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(Ndist)}}number of distinct values{p_end}
{synopt:{cmd:e(Nclust)}}number of clusters{p_end}
{synopt:{cmd:e(deriv)}}order of derivative{p_end}
{synopt:{cmd:e(imse_bsq_rot)}}bias constant in IMSE, ROT selection{p_end}
{synopt:{cmd:e(imse_var_rot)}}variance constant in IMSE, ROT selection{p_end}
{synopt:{cmd:e(imse_bsq_dpi)}}bias constant in IMSE, DPI selection{p_end}
{synopt:{cmd:e(imse_var_dpi)}}variance constant in IMSE, DPI selection{p_end}
{synopt:{cmd:e(nbinsrot_poly)}}ROT number of bins, unregularized{p_end}
{synopt:{cmd:e(nbinsrot_regul)}}ROT number of bins, regularized or user-specified{p_end}
{synopt:{cmd:e(nbinsrot_uknot)}}ROT number of bins, unique knots{p_end}
{synopt:{cmd:e(nbinsdpi)}}DPI number of bins{p_end}
{synopt:{cmd:e(nbinsdpi_uknot)}}DPI number of bins, unique knots{p_end}
{synopt:{cmd:e(prot_poly)}}ROT degree of polynomial, unregularized{p_end}
{synopt:{cmd:e(prot_regul)}}ROT degree of polynomial, regularized or user-specified{p_end}
{synopt:{cmd:e(prot_uknot)}}ROT degree of polynomial, unique knots{p_end}
{synopt:{cmd:e(pdpi)}}DPI degree of polynomial{p_end}
{synopt:{cmd:e(pdpi_uknot)}}DPI degree of polynomial, unique knots{p_end}
{synopt:{cmd:e(srot_poly)}}ROT number of smoothness constraints, unregularized{p_end}
{synopt:{cmd:e(srot_regul)}}ROT number of smoothness constraints, regularized or user-specified{p_end}
{synopt:{cmd:e(srot_uknot)}}ROT number of smoothness constraints, unique knots{p_end}
{synopt:{cmd:e(sdpi)}}DPI number of smoothness constraints{p_end}
{synopt:{cmd:e(sdpi_uknot)}}DPI number of smoothness constraints, unique knots{p_end}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(knot)}}numlist of knots{p_end}
{synopt:{cmd:e(m_p)}}vector of degrees of polynomial{p_end}
{synopt:{cmd:e(m_s)}}vector of number of smoothness constraints{p_end}
{synopt:{cmd:e(m_nbinsrot_poly)}}ROT number of bins, unregularized, for each pair of degree and smoothness{p_end}
{synopt:{cmd:e(m_nbinsrot_regul)}}ROT number of bins, regularized or user-specified, for each pair of degree and smoothness{p_end}
{synopt:{cmd:e(m_nbinsrot_uknot)}}ROT number of bins, unique knots, for each pair of degree and smoothness{p_end}
{synopt:{cmd:e(m_nbinsdpi)}}DPI number of bins, for each pair of degree and smoothness{p_end}
{synopt:{cmd:e(m_nbinsdpi_uknot)}}DPI number of bins, unique knots, for each pair of degree and smoothness{p_end}
{synopt:{cmd:e(m_imse_bsq_rot)}}bias constant in IMSE, ROT selection, for each pair of degree and smoothness{p_end}
{synopt:{cmd:e(m_imse_var_rot)}}variance constant in IMSE, ROT selection, for each pair of degree and smoothness{p_end}
{synopt:{cmd:e(m_imse_bsq_dpi)}}bias constant in IMSE, DPI selection, for each pair of degree and smoothness{p_end}
{synopt:{cmd:e(m_imse_var_dpi)}}variance constant in IMSE, DPI selection, for each pair of degree and smoothness{p_end}

{marker references}{...}
{title:References}

{p 4 8} Cattaneo, M. D., R. K. Crump, M. H. Farrell, and Y. Feng. 2023a.
{browse "https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2023_AER.pdf":On Binscatter}.
Working Paper.
{p_end}

{p 4 8} Cattaneo, M. D., R. K. Crump, M. H. Farrell, and Y. Feng. 2023b.
{browse "https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2023_NonlinearBinscatter.pdf":Nonlinear Binscatter Methods}.
Working Paper.
{p_end}

{p 4 8} Cattaneo, M. D., R. K. Crump, M. H. Farrell, and Y. Feng. 2023c.
{browse "https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2023_Stata.pdf":Binscatter Regressions}.
Working Paper.
{p_end}


{marker authors}{...}
{title:Authors}

{p 4 8} Matias D. Cattaneo, Princeton University, Princeton, NJ.
{browse "mailto:cattaneo@princeton.edu":cattaneo@princeton.edu}.
{p_end}

{p 4 8} Richard K. Crump, Federal Reserve Band of New York, New York, NY.
{browse "mailto:richard.crump@ny.frb.org":richard.crump@ny.frb.org}.
{p_end}

{p 4 8} Max H. Farrell, UC Santa Barbara, Santa Barbara, CA.
{browse "mailto:mhfarrell@gmail.com":mhfarrell@gmail.com}.
{p_end}

{p 4 8} Yingjie Feng, Tsinghua University, Beijing, China.
{browse "mailto:fengyingjiepku@gmail.com":fengyingjiepku@gmail.com}.
{p_end}


{smcl}
{* *! version 0.3 10-JUN-2021}{...}
{viewerjumpto "Syntax" "binspwc##syntax"}{...}
{viewerjumpto "Description" "binspwc##description"}{...}
{viewerjumpto "Options" "binspwc##options"}{...}
{viewerjumpto "Examples" "binspwc##examples"}{...}
{viewerjumpto "Stored results" "binspwc##stored_results"}{...}
{viewerjumpto "References" "binspwc##references"}{...}
{viewerjumpto "Authors" "binspwc##authors"}{...}
{cmd:help binspwc}
{hline}

{title:Title}

{p 4 8}{hi:binspwc} {hline 2} Data-Driven Nonparametric Pairwise Group Comparison using Binscatter.{p_end}


{marker syntax}{...}
{title:Syntax}

{p 4 16} {cmdab:binspwc} {depvar} {it:indvar} [{it:covars}] {ifin} {weight} [ {cmd:,} {opt estmethod(cmdname)} {opt deriv(v)} {opt by(varname)}{p_end}
{p 16 16} {opt pwc(p s)} {opt testtype(type)} {opt lp(metric)}{p_end}
{p 16 16} {opt bins(p s)} {opt bynbins(numlist)} {opt binspos(position)} {opt binsmethod(method)} {opt nbinsrot(#)} {opt samebinsby}{p_end}
{p 16 16} {opt nsims(#)} {opt simsgrid(#)} {opt simsseed(seed)}{p_end}
{p 16 16} {opt dfcheck(n1 n2)} {opt masspoints(masspointsoption)}{p_end}
{p 16 16} {cmd:vce(}{it:{help vcetype}}{cmd:)} ]{p_end}

{p 4 8} where {depvar} is the dependent variable, {it:indvar} is the independent variable for binning, and {it:covars} are other covariates to be controlled for.{p_end}

{p 4 8} p, s and v are integers satisfying 0 <= s,v <= p, which can take different values in each case.{p_end}

{p 4 8} {opt fweight}s, {opt aweight}s and {opt pweight}s are allowed; see {help weight}.{p_end}

{marker description}{...}
{title:Description}

{p 4 8} {cmd:binspwc} implements binscatter-based hypothesis testing procedures for pairwise group comparison of binscatter estimators, following the results in
{browse "https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2019_Binscatter.pdf":Cattaneo, Crump, Farrell and Feng (2021a)}.
If the binning scheme is not set by the user, the companion command {help binsregselect:binsregselect} is used to implement binscatter in a data-driven (optimal) way and inference procedures are based on robust bias correction.
Binned scatter plots based on different models can be constructed using the companion commands {help binsreg:binsreg}, {help binsqreg: binsqreg}, {help binslogit:binslogit} and {help binsprobit:binsprobit}.
{p_end}

{p 4 8} A detailed introduction to this command is given in
{browse "https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2019_Stata.pdf":Cattaneo, Crump, Farrell and Feng (2021b)}.
A companion R package with the same capabilities is available (see website below).
{p_end}

{p 4 8} Companion commands: {help binsreg:binsreg} for binscatter least squares regression with robust inference procedures and plots, {help binsqreg:binsqreg} for binscatter quantile regression with robust inference procedures and plots, {help binslogit:binslogit} for binscatter logit estimation with robust inference procedures and plots, {help binsprobit:binsprobit} for binscatter probit estimation with robust inference procedures and plots, and {help binsregselect:binsregselect} data-driven (optimal) binning selection.{p_end}

{p 4 8} Related Stata and R packages are available in the following website:{p_end}

{p 8 8} {browse "https://nppackages.github.io/":https://nppackages.github.io/}{p_end}


{marker options}{...}
{title:Options}

{dlgtab:Estimand}

{p 4 8} {opt estmethod(cmdname)} specifies the binscatter model. The default is {cmd:estmethod(reg)}, which corresponds to the binscatter least squares regression. Other options are: {cmd:estmethod(qreg #)} for binscatter quantile regression where # is the quantile to be estimated, {cmd:estmethod(logit)} for binscatter logistic regression and {cmd:estmethod(probit)} for binscatter probit regression.
{p_end}  

{p 4 8} {opt deriv(v)} specifies the derivative order of the regression function for estimation, testing and plotting.
The default is {cmd:deriv(0)}, which corresponds to the function itself.
{p_end}

{p 4 8} {opt by(varname)} specifies the variable containing the group indicator to perform subgroup analysis; both numeric and string variables are supported.  When {opt by(varname)} is specified, {cmdab:binspwc} implements estimation by each subgroup separately and then conduct {it:all} pairwise comparison tests. By default, the binning structure is selected for each subgroup separately, but see the option samebinsby below for imposing a common binning structure across subgroups.

{dlgtab:Pairwise Group Comparison Testing}

{p 4 8} {opt pwc(p s)} sets a piecewise polynomial of degree {it:p} with {it:s} smoothness constraints for pairwise group comparison.
The default is {cmd:pwc(3 3)}, which corresponds to a cubic B-spline estimate of the function of interest for each group.
{p_end}

{p 4 8} {opt testtype(type)} specifies the type of pairwise comparison test. The default is {opt testtype(2)}, which corresponds to a two-sided test of the form H0: {it:mu_1(x)=mu_2(x)}. Other options are: {opt testtype(l)} for the one-sided test of the form H0: {it:mu_1(x)<=mu_2(x)} and {opt testtype(r)} for the one-sided test of the form H0: {it:mu_1(x)>=mu_2(x)}.
{p_end}

{p 4 8} {opt lp(metric)} specifies a Lp metric used for a (two-sided) test for the difference between two groups. The default is {cmd:lp(inf)}, which corresponds to the sup-norm. Other options are Lp(q) for a positive integer q.
{p_end}
 
{dlgtab:Partitioning/Binning Selection}

{p 4 8} {opt bins(p s)} sets a piecewise polynomial of degree {it:p} with {it:s} smoothness constraints for data-driven (IMSE-optimal) selection of the partitioning/binning scheme.
The default is {cmd:bins(0 0)}, which corresponds to piecewise constant (canonical binscatter).

{p 4 8} {opt bynbins(numlist)} sets a {help numlist} of numbers of bins for partitioning/binning of {it:indvar}, which is applied to the binscatter estimation for each group. The ordering of the group follows the result of the {help tabulate oneway:tabulate}. If a single number of bins is specified, it applies to the estimation for all groups. If not specified, the number of bins is selected via the companion command {help binsregselect:binsregselect} in a data-driven, optimal way whenever possible.
{p_end}

{p 4 8} {opt binspos(position)} specifies the position of binning knots.
The default is {cmd:binspos(qs)}, which corresponds to quantile-spaced binning (canonical binscatter).
Other options are: {cmd:es} for evenly-spaced binning, or a {help numlist} for manual specification of the positions of inner knots (which must be within the range of {it:indvar}).
{p_end}

{p 4 8} {opt binsmethod(method)} specifies the method for data-driven selection of the number of bins via the companion command {help binsregselect:binsregselect}.
The default is {cmd:binsmethod(dpi)}, which corresponds to the IMSE-optimal direct plug-in rule.
The other option is: {cmd:rot} for rule of thumb implementation.
{p_end}

{p 4 8} {opt nbinsrot(#)} specifies an initial number of bins value used to construct the DPI number of bins selector.
If not specified, the data-driven ROT selector is used instead.
{p_end}

{p 4 8} {opt samebinsby} forces a common partitioning/binning structure across all subgroups specified by the option {cmd:by()}.
The knots positions are selected according to the option {cmd:binspos()} and using the full sample.
If {cmd:nbins()} is not specified, then the number of bins is selected via the companion command {help binsregselect:binsregselect} and using the full sample.{p_end}

{dlgtab:Simulation}

{p 4 8} {opt nsims(#)} specifies the number of random draws for constructing confidence bands and hypothesis testing.
The default is {cmd:nsims(500)}, which corresponds to 500 draws from a standard Gaussian random vector of size [(p+1)*J - (J-1)*s].
{p_end}

{p 4 8} {opt simsgrid(#)} specifies the number of evaluation points of an evenly-spaced grid within each bin used for evaluation of the supremum (or infimum) operation needed to construct confidence bands and hypothesis testing procedures.
The default is {cmd:simsgrid(20)}, which corresponds to 20 evenly-spaced evaluation points within each bin for approximating the supremum (or infimum) operator.
{p_end}

{p 4 8} {opt simsseed(#)} sets the seed for simulations.
{p_end}

{dlgtab:Mass Points and Degrees of Freedom}

{p 4 8} {opt dfcheck(n1 n2)} sets cutoff values for minimum effective sample size checks, which take into account the number of unique values of {it:indvar} (i.e., adjusting for the number of mass points), number of clusters, and degrees of freedom of the different statistical models considered.
The default is {cmd:dfcheck(20 30)}. See Cattaneo, Crump, Farrell and Feng (2019b) for more details.
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

{p 4 8} {cmd:vce(}{it:{help vcetype}}{cmd:)} specifies the {it:vcetype} for variance estimation used by the commands {help regress##options:regress}, {help logit##options:logit} or {help qreg##qreg_options:qreg}.
The default is {cmd:vce(robust)}.
{p_end}


{marker examples}{...}
{title:Examples}

{p 4 8} Test the difference between two groups{p_end}
{p 8 8} . {stata binspwc y x w, by(t)}{p_end}


{marker stored_results}{...}
{title:Stored results}

{synoptset 17 tabbed}{...}
{p2col 5 17 21 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(p)}}degree of polynomial for bin selection{p_end}
{synopt:{cmd:e(s)}}smoothness of polynomial for bin selection{p_end}
{synopt:{cmd:e(pwc_p)}}degree of polynomial for testing{p_end}
{synopt:{cmd:e(pwc_s)}}smoothnes of polynomial for testing{p_end}
{p2col 5 17 21 2: Locals}{p_end}
{synopt:{cmd:e(byvalue)}}name of groups found in {cmd:by()}{p_end}
{p2col 5 17 21 2: Matrices}{p_end}
{synopt:{cmd:e(N_by)}}number of observations for each group{p_end}
{synopt:{cmd:e(Ndist_by)}}number of distinct values for each group{p_end}
{synopt:{cmd:e(Nclust_by)}}number of clusters for each group{p_end}
{synopt:{cmd:e(nbins_by)}}number of bins for each group{p_end}
{synopt:{cmd:e(stat)}}test statistics for all pairwise comparisons{p_end}
{synopt:{cmd:e(pval)}}p values for all pairwise comparisons{p_end}


{marker references}{...}
{title:References}

{p 4 8} Cattaneo, M. D., R. K. Crump, M. H. Farrell, and Y. Feng. 2021a.
{browse "https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2019_Binscatter.pdf":On Binscatter}.
{it:arXiv:1902.09608}.
{p_end}

{p 4 8} Cattaneo, M. D., R. K. Crump, M. H. Farrell, and Y. Feng. 2021b.
{browse "https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2019_Stata.pdf":Binscatter Regressions}.
{it:arXiv:1902.09615}.
{p_end}


{marker authors}{...}
{title:Authors}

{p 4 8} Matias D. Cattaneo, Princeton University, Princeton, NJ.
{browse "mailto:cattaneo@princeton.edu":cattaneo@princeton.edu}.
{p_end}

{p 4 8} Richard K. Crump, Federal Reserve Band of New York, New York, NY.
{browse "mailto:richard.crump@ny.frb.org":richard.crump@ny.frb.org}.
{p_end}

{p 4 8} Max H. Farrell, University of Chicago, Chicago, IL.
{browse "mailto:max.farrell@chicagobooth.edu":max.farrell@chicagobooth.edu}.
{p_end}

{p 4 8} Yingjie Feng, Princeton University, Princeton, NJ.
{browse "mailto:yingjief@princeton.edu":yingjief@princeton.edu}.
{p_end}


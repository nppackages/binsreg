{smcl}
{* *! version 0.8 12-OCT-2021}{...}
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

{p 4 12} {cmdab:binspwc} {depvar} {it:indvar} [{it:othercovs}] {ifin} {weight} {cmd:,} {opt by(varname)} [{p_end}
{p 12 12} {opt estmethod(cmdname)} {opt deriv(v)} {opt at(position)} {opt nolink}{p_end}
{p 12 12} {opt absorb(absvars)} {opt reghdfeopt(reghdfe_option)}{p_end}
{p 12 12} {opt pwc(p s)} {opt testtype(type)} {opt lp(metric)}{p_end}
{p 12 12} {opt bins(p s)} {opt bynbins(numlist)} {opt binspos(position)} {opt binsmethod(method)} {opt nbinsrot(#)} {opt samebinsby} {opt randcut(#)}{p_end}
{p 12 12} {opt nsims(#)} {opt simsgrid(#)} {opt simsseed(seed)}{p_end}
{p 12 12} {opt dfcheck(n1 n2)} {opt masspoints(masspointsoption)}{p_end}
{p 12 12} {cmd:vce(}{it:{help vcetype}}{cmd:)} {opt asyvar(on/off)} {opt usegtools(on/off)} ]{p_end}

{p 4 8} where {depvar} is the dependent variable, {it:indvar} is the independent variable for binning, and {it:othercovs} are other covariates to be controlled for.{p_end}

{p 4 8} p, s and v are integers satisfying 0 <= s,v <= p, which can take different values in each case.{p_end}

{p 4 8} {opt fweight}s, {opt aweight}s and {opt pweight}s are allowed; see {help weight}.{p_end}

{marker description}{...}
{title:Description}

{p 4 8} {cmd:binspwc} implements binscatter-based hypothesis testing procedures for pairwise group comparison of binscatter estimators, following the results in
{browse "https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2021_Binscatter.pdf":Cattaneo, Crump, Farrell and Feng (2021a)}.
If the binning scheme is not set by the user, the companion command {help binsregselect:binsregselect} is used to implement binscatter
in a data-driven (optimal) way and inference procedures are based on robust bias correction.
Binned scatter plots based on different models can be constructed using the companion commands {help binsreg:binsreg},
{help binsqreg: binsqreg}, {help binslogit:binslogit} and {help binsprobit:binsprobit}.
{p_end}

{p 4 8} A detailed introduction to this command is given in
{browse "https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2021_Stata.pdf":Cattaneo, Crump, Farrell and Feng (2021b)}.
Companion R and Python packages with the same capabilities are available (see website below).
{p_end}

{p 4 8} Companion commands: {help binsreg:binsreg} for binscatter least squares regression with robust inference procedures and plots,
{help binsqreg:binsqreg} for binscatter quantile regression with robust inference procedures and plots,
{help binslogit:binslogit} for binscatter logit estimation with robust inference procedures and plots,
{help binsprobit:binsprobit} for binscatter probit estimation with robust inference procedures and plots, and
{help binsregselect:binsregselect} for data-driven (optimal) binning selection.{p_end}

{p 4 8} Related Stata, R and Python packages are available in the following website:{p_end}

{p 8 8} {browse "https://nppackages.github.io/":https://nppackages.github.io/}{p_end}


{marker options}{...}
{title:Options}

{dlgtab:Estimand}

{p 4 8} {opt by(varname)} specifies the variable containing the group indicator to perform subgroup analysis; both numeric and string variables are supported.
When {opt by(varname)} is specified, {cmdab:binspwc} implements estimation for each subgroup separately and then conduct {it:all} pairwise comparison tests.
By default, the binning structure is selected for each subgroup separately, but see the option
{cmd:samebinsby} below for imposing a common binning structure across subgroups.
This option is required.
{p_end}

{p 4 8} {opt estmethod(cmdname)} specifies the binscatter model. The default is {cmd:estmethod(reg)},
which corresponds to the binscatter least squares regression. Other options are:
{cmd:estmethod(qreg #)} for binscatter quantile regression where # is the quantile to be estimated,
{cmd:estmethod(logit)} for binscatter logistic regression and {cmd:estmethod(probit)} for binscatter probit regression.
{p_end}

{p 4 8} {opt deriv(v)} specifies the derivative order of the regression function for estimation, testing and plotting.
The default is {cmd:deriv(0)}, which corresponds to the function itself.
{p_end}

{p 4 8} {opt at(position)} specifies the values of {it:othercovs} at which the estimated function is evaluated for plotting.
The default is {cmd:at(mean)}, which corresponds to the mean of {it:othercovs}. Other options are: {cmd:at(median)} for the
median of {it:othercovs}, {cmd:at(0)} for zeros, and {cmd:at(filename)} for particular values of {it:othercovs} saved in another file.
{p_end}

{p 4 8} Note: when {cmd:at(mean)} or {cmd:at(median)} is specified, all factor variables in {it:othercovs} (if specified)
are excluded from the evaluation (set as zero).
{p_end}

{p 4 8}{opt nolink} specifies that the function within the inverse link (logistic) function be reported instead of
the conditional probability function. This option is used only if logit or probit model is specified in {cmd:estmethod()}.
{p_end}

{dlgtab:Reghdfe}

{p 4 8} {opt absorb(absvars)} specifies categorical variables (or interactions) representing the fixed effects to be absorbed.
This is equivalent to including an indicator/dummy variable for each category of each {it:absvar}.
When {cmd:absorb()} is specified, the community-contributed command {cmd:reghdfe} instead of the command {cmd:regress} is used.
{p_end}

{p 4 8} {opt reghdfeopt(reghdfe_option)} options to be passed on to the command {cmd:reghdfe}. 
Important: {cmd:absorb()} and {cmd:vce()} should not be specified within this option.
{p_end}

{p 4 8} For more information about the community-contributed command {cmd:reghdfe},
please see {browse "http://scorreia.com/software/reghdfe/":http://scorreia.com/software/reghdfe/}.

{dlgtab:Pairwise Group Comparison Testing}

{p 4 8} {opt pwc(p s)} sets a piecewise polynomial of degree {it:p} with {it:s} smoothness constraints for pairwise group comparison.
The default is {cmd:pwc(3 3)}, which corresponds to a cubic B-spline estimate of the function of interest for each group.
{p_end}

{p 4 8} {opt testtype(type)} specifies the type of pairwise comparison test. The default is {opt testtype(2)},
which corresponds to a two-sided test of the form H0: {it:mu_1(x)=mu_2(x)}. Other options are: {opt testtype(l)}
for the one-sided test of the form H0: {it:mu_1(x)<=mu_2(x)} and {opt testtype(r)} for the one-sided test of the form H0: {it:mu_1(x)>=mu_2(x)}.
{p_end}

{p 4 8} {opt lp(metric)} specifies an Lp metric used for a (two-sided) test for the difference between two groups. The default is {cmd:lp(inf)},
which corresponds to the sup-norm. Other options are {cmd:Lp(q)} for a positive integer {cmd:q}.
{p_end}
 
{dlgtab:Partitioning/Binning Selection}

{p 4 8} {opt bins(p s)} sets a piecewise polynomial of degree {it:p} with {it:s} smoothness constraints for data-driven (IMSE-optimal)
selection of the partitioning/binning scheme.
The default is {cmd:bins(2 2)}, which corresponds to a quadratic spline estimate.

{p 4 8} {opt bynbins(numlist)} sets a {help numlist} of numbers of bins for partitioning/binning of {it:indvar},
which is applied to the binscatter estimation for each group. The ordering of the group follows
the result of {help tabulate oneway:tabulate}. If a single number of bins is specified, it applies to the estimation for all groups.
If not specified, the number of bins is selected via the companion command {help binsregselect:binsregselect} in a data-driven, optimal way whenever possible.
{p_end}

{p 4 8} {opt binspos(position)} specifies the position of binning knots.
The default is {cmd:binspos(qs)}, which corresponds to quantile-spaced binning (canonical binscatter).
Other options are: {cmd:es} for evenly-spaced binning, or a {help numlist} for manual specification of
the positions of inner knots (which must be within the range of {it:indvar}).
{p_end}

{p 4 8} {opt binsmethod(method)} specifies the method for data-driven selection of the number of bins via
the companion command {help binsregselect:binsregselect}.
The default is {cmd:binsmethod(dpi)}, which corresponds to the IMSE-optimal direct plug-in rule.
The other option is: {cmd:rot} for rule of thumb implementation.
{p_end}

{p 4 8} {opt nbinsrot(#)} specifies an initial number of bins value used to construct the DPI number of bins selector.
If not specified, the data-driven ROT selector is used instead.
{p_end}

{p 4 8} {opt samebinsby} forces a common partitioning/binning structure across all subgroups specified by the option {cmd:by()}.
The knots positions are selected according to the option {cmd:binspos()} and using the full sample.
If {cmd:nbins()} is not specified, then the number of bins is selected via the companion
command {help binsregselect:binsregselect} and using the full sample.{p_end}

{p 4 8} {opt randcut(#)} specifies the upper bound on a uniformly distributed variable used to draw a subsample for bins selection.
Observations for which {cmd:runiform()<=#} are used. # must be between 0 and 1.{p_end}

{dlgtab:Simulation}

{p 4 8} {opt nsims(#)} specifies the number of random draws for hypothesis testing.
The default is {cmd:nsims(500)}, which corresponds to 500 draws from a standard Gaussian random vector of size [(p+1)*J - (J-1)*s].
{p_end}

{p 4 8} {opt simsgrid(#)} specifies the number of evaluation points of an evenly-spaced grid
within each bin used for evaluation of the supremum (infimum or Lp metric) operation needed to
construct confidence bands and hypothesis testing procedures.
The default is {cmd:simsgrid(20)}, which corresponds to 20 evenly-spaced evaluation points
within each bin for approximating the supremum (infimum or Lp metric) operator.
{p_end}

{p 4 8} {opt simsseed(#)} sets the seed for simulations.
{p_end}

{dlgtab:Mass Points and Degrees of Freedom}

{p 4 8} {opt dfcheck(n1 n2)} sets cutoff values for minimum effective sample size checks,
which take into account the number of unique values of {it:indvar} (i.e., adjusting for the number of mass points),
number of clusters, and degrees of freedom of the different statistical models considered.
The default is {cmd:dfcheck(20 30)}. See Cattaneo, Crump, Farrell and Feng (2021b) for more details.
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

{p 4 8} {cmd:vce(}{it:{help vcetype}}{cmd:)} specifies the {it:vcetype} for variance estimation used by
the commands {help regress##options:regress},
{help logit##options:logit}, {help logit##options:logit}, {help qreg##qreg_options:qreg} or {cmd:reghdfe}.
The default is {cmd:vce(robust)}.
{p_end}

{p 4 8} {opt asyvar(on/off)} specifies the method used to compute standard errors.
If {cmd:asyvar(on)} is specified, the standard error of the nonparametric component is used and the uncertainty
related to other control variables {it:othercovs} is omitted.
Default is {cmd:asyvar(off)}, that is, the uncertainty related to {it:othercovs} is taken into account.
{p_end}

{p 4 8}{opt usegtools(on/off)} forces the use of several commands in the community-distributed Stata package
{cmd:gtools} to speed the computation up, if {it:on} is specified.
Default is {cmd:usegtools(off)}.
{p_end}

{p 4 8} For more information about the package {cmd:gtools}, please see {browse "https://gtools.readthedocs.io/en/latest/index.html":https://gtools.readthedocs.io/en/latest/index.html}.
{p_end}

{marker examples}{...}
{title:Examples}

{p 4 8} Setup{p_end}
{p 8 8} . {stata  sysuse auto}{p_end}

{p 4 8} Generate two groups{p_end}
{p 8 8} . {stata  gen group=price>5000}{p_end}

{p 4 8} Test for the difference between two groups{p_end}
{p 8 8} . {stata binspwc mpg weight foreign, by(group)}{p_end}


{marker stored_results}{...}
{title:Stored results}

{synoptset 17 tabbed}{...}
{p2col 5 17 21 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(p)}}degree of polynomial for bin selection{p_end}
{synopt:{cmd:e(s)}}smoothness of polynomial for bin selection{p_end}
{synopt:{cmd:e(pwc_p)}}degree of polynomial for testing{p_end}
{synopt:{cmd:e(pwc_s)}}smoothness of polynomial for testing{p_end}
{p2col 5 17 21 2: Macros}{p_end}
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
{browse "https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2021_Binscatter.pdf":On Binscatter}.
{it:arXiv:1902.09608}.
{p_end}

{p 4 8} Cattaneo, M. D., R. K. Crump, M. H. Farrell, and Y. Feng. 2021b.
{browse "https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2021_Stata.pdf":Binscatter Regressions}.
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

{p 4 8} Yingjie Feng, Tsinghua University, Beijing, China.
{browse "mailto:fengyingjiepku@gmail.com":fengyingjiepku@gmail.com}.
{p_end}


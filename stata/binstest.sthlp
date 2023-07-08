{smcl}
{* *! version 1.3 03-JUL-2023}{...}
{viewerjumpto "Syntax" "binstest##syntax"}{...}
{viewerjumpto "Description" "binstest##description"}{...}
{viewerjumpto "Options" "binstest##options"}{...}
{viewerjumpto "Examples" "binstest##examples"}{...}
{viewerjumpto "Stored results" "binstest##stored_results"}{...}
{viewerjumpto "References" "binstest##references"}{...}
{viewerjumpto "Authors" "binstest##authors"}{...}
{cmd:help binstest}
{hline}

{title:Title}

{p 4 8}{hi:binstest} {hline 2} Data-Driven Nonparametric Shape Restriction and Parametric Model Specification Testing using Binscatter.{p_end}


{marker syntax}{...}
{title:Syntax}

{p 4 13} {cmdab:binstest} {depvar} {it:indvar} [{it:othercovs}] {ifin} {weight} [ {cmd:,} {p_end}
{p 13 13} {opt estmethod(cmdname)} {opt deriv(v)} {opt at(position)} {opt nolink}{p_end}
{p 13 13} {opt absorb(absvars)} {opt reghdfeopt(reghdfe_option)}{p_end}
{p 13 13} {opt testmodel(testmodelopt)} {opt testmodelparfit(filename)} {opt testmodelpoly(p)}{p_end}
{p 13 13} {opt testshape(testshapeopt)} {opt testshapel(numlist)} {opt testshaper(numlist)} {opt testshape2(numlist)} {opt lp(metric)}{p_end}
{p 13 13} {opt bins(p s)} {opt nbins(nbinsopt)} {opt binspos(position)} {opt binsmethod(method)} {opt nbinsrot(#)} {opt randcut(#)}{p_end}
{p 13 13} {cmd:pselect(}{it:{help numlist}}{cmd:)} {cmd:sselect(}{it:{help numlist}}{cmd:)}{p_end}
{p 13 13} {opt nsims(#)} {opt simsgrid(#)} {opt simsseed(seed)}{p_end}
{p 13 13} {opt dfcheck(n1 n2)} {opt masspoints(masspointsoption)}{p_end}
{p 13 13} {cmd:vce(}{it:{help vcetype}}{cmd:)} {opt asyvar(on/off)} {opt estmethodopt(cmd_option)} {opt usegtools(on/off)} ]{p_end}

{p 4 8} where {depvar} is the dependent variable, {it:indvar} is the independent variable for binning, and {it:othercovs}
are other covariates to be controlled for.{p_end}

{p 4 8} The degree of the piecewise polynomial p, the number of smoothness constraints s, and the derivative order v are integers 
satisfying 0 <= s,v <= p, which can take different values in each case.{p_end}

{p 4 8} At least one test has to be specified via {opt testmodelparfit()}, {opt testmodelpoly()}, {opt testshapel()},
{opt testshaper()} and/or {opt testshape2()}.
{p_end}

{p 4 8} {opt fweight}s, {opt aweight}s and {opt pweight}s are allowed; see {help weight}.{p_end}

{marker description}{...}
{title:Description}

{p 4 8} {cmd:binstest} implements binscatter-based hypothesis testing procedures for parametric functional forms of
and nonparametric shape restrictions on the regression function estimators, following the results in
{browse "https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2023_AER.pdf":Cattaneo, Crump, Farrell and Feng (2023a)} and
{browse "https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2023_NonlinearBinscatter.pdf":Cattaneo, Crump, Farrell and Feng (2023b)}.
If the binning scheme is not set by the user, the companion command {help binsregselect:binsregselect} is used
to implement binscatter in a data-driven (optimal) way and inference procedures are based on robust bias correction.
Binned scatter plots based on different models can be constructed using the companion commands {help binsreg:binsreg},
{help binsqreg: binsqreg}, {help binslogit:binslogit} and {help binsprobit:binsprobit}.
{p_end}

{p 4 8} A detailed introduction to this command is given in
{browse "https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2023_Stata.pdf":Cattaneo, Crump, Farrell and Feng (2023c)}.
Companion R and Python packages with the same capabilities are available (see website below).
{p_end}

{p 4 8} Companion commands: {help binsreg:binsreg} for binscatter regression with robust inference procedures and plots,
{help binsqreg:binsqreg} for binscatter quantile regression with robust inference procedures and plots,
{help binslogit:binslogit} for binscatter logit estimation with robust inference procedures and plots,
{help binsprobit:binsprobit} for binscatter probit estimation with robust inference procedures and plots,
and {help binsregselect:binsregselect} for data-driven (optimal) binning selection.{p_end}

{p 4 8} Related Stata, R and Python packages are available in the following website:{p_end}

{p 8 8} {browse "https://nppackages.github.io/":https://nppackages.github.io/}{p_end}


{marker options}{...}
{title:Options}

{dlgtab:Estimand}

{p 4 8} {opt estmethod(cmdname)} specifies the binscatter model. The default is {cmd:estmethod(reg)},
which corresponds to the binscatter least squares regression. Other options are: {cmd:estmethod(qreg #)}
for binscatter quantile regression where # is the quantile to be estimated, {cmd:estmethod(logit)} for
binscatter logistic regression and {cmd:estmethod(probit)} for binscatter probit regression.
{p_end}

{p 4 8} {opt deriv(v)} specifies the derivative order of the regression function for estimation, testing and plotting.
The default is {cmd:deriv(0)}, which corresponds to the function itself.
{p_end}

{p 4 8} {opt at(position)} specifies the values of {it:othercovs} at which the estimated function is evaluated for plotting.
The default is {cmd:at(mean)}, which corresponds to the mean of {it:othercovs}. Other options are: {cmd:at(median)} for the median of {it:othercovs},
{cmd:at(0)} for zeros, and {cmd:at(filename)} for particular values of {it:othercovs} saved in another file.
{p_end}

{p 4 8} Note: When {cmd:at(mean)} or {cmd:at(median)} is specified, all factor variables in {it:othercovs} (if specified)
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

{p 4 8} For more information about the community-contributed command {cmd:reghdfe}, please see {browse "http://scorreia.com/software/reghdfe/":http://scorreia.com/software/reghdfe/}.

{dlgtab:Parametric Model Specification Testing}

{p 4 8} {opt testmodel(testmodelopt)} sets the degree of polynomial and the number of smoothness constraints for parametric model specification testing. 
If {cmd:testmodel(p s)} is specified, a piecewise polynomial of degree {it:p} with {it:s} smoothness constraints is used.
If {cmd:testmodel(T)} or {cmd:testmodel()} is specified, 
{cmd:testmodel(1 1)} is used unless the degree {it:p} or smoothness {it:s} selection
is requested via the option {cmd:pselect()} or {cmd:sselect()} (see more details in the explanation 
of {cmd:pselect()} and {cmd:sselect()}). 
The default is {cmd:testmodel()}.
{p_end}

{p 4 8} {opt testmodelparfit(filename)} specifies a dataset which contains the evaluation grid and fitted values of the model(s) to be tested against.
The file must have a variable with the same name as {it:indvar}, which contains a series of evaluation points at which
the binscatter model and the parametric model of interest are compared with each other.
Each parametric model is represented by a variable named as {it:binsreg_fit*}, which must contain the fitted values at the corresponding evaluation points.
{p_end}

{p 4 8} {opt testmodelpoly(p)} specifies the degree of a global polynomial model to be tested against.
{p_end}
 
{dlgtab:Nonparametric Shape Restriction Testing}

{p 4 8} {opt testshape(testshapeopt)} sets the degree of polynomial and the number of smoothness constraints
for nonparametric shape restriction testing. If {cmd:testshape(p s)} is specified, 
a piecewise polynomial of degree {it:p} with {it:s} smoothness constraints is used. 
If {cmd:testshape(T)} or {cmd:testshape()} is specified, 
{cmd:testshape(1 1)} is used unless the degree {it:p} or smoothness {it:s} selection
is requested via the option {cmd:pselect()} or {cmd:sselect()} (see more details in the explanation 
of {cmd:pselect()} and {cmd:sselect()}). 
The default is {cmd:testshape()}.
{p_end}

{p 4 8} {opt testshapel(numlist)} specifies a {help numlist} of null boundary values for hypothesis testing.
Each number {it:a} in the {it:numlist} corresponds to one boundary of a one-sided hypothesis test to the left of the form H0: {it:sup_x mu(x)<=a}.
{p_end}

{p 4 8} {opt testshaper(numlist)} specifies a {help numlist} of null boundary values for hypothesis testing.
Each number {it:a} in the {it:numlist} corresponds to one boundary of a one-sided hypothesis test to the right of the form H0: {it:inf_x mu(x)>=a}.
{p_end}

{p 4 8} {opt testshape2(numlist)} specifies a {help numlist} of null boundary values for hypothesis testing.
Each number {it:a} in the {it:numlist} corresponds to one boundary of a two-sided hypothesis test of the
form H0: {it:sup_x |mu(x)-a|=0}.
{p_end}

{dlgtab:Metric for Hypothesis Testing}

{p 4 8} {opt lp(metric)} specifies an Lp metric used for (two-sided) parametric model specification testing and/or shape restriction testing.
The default is {cmd:lp(inf)},
which corresponds to the sup-norm. Other options are {cmd:lp(q)} for a positive integer {cmd:q}.
{p_end}
 
{dlgtab:Binning/Degree/Smoothness Selection}

{p 4 8} {opt bins(p s)} sets a piecewise polynomial of degree {it:p} with {it:s} smoothness constraints for
data-driven (IMSE-optimal) selection of the partitioning/binning scheme.
The default is {cmd:bins(0 0)}, which corresponds to the piecewise constant.

{p 4 8} {opt nbins(nbinsopt)} sets the number of bins for partitioning/binning of {it:indvar}. 
If {cmd:nbins(T)} or {cmd:nbins()} (default) is specified, the number of bins is selected via the companion command {help binsregselect:binsregselect} 
in a data-driven, optimal way whenever possible. If a {help numlist:numlist} with more than one number is specified, 
the number of bins is selected within this list via the companion command {help binsregselect:binsregselect}.
{p_end}

{p 4 8} {opt binspos(position)} specifies the position of binning knots.
The default is {cmd:binspos(qs)}, which corresponds to quantile-spaced binning (canonical binscatter).
Other options are: {cmd:es} for evenly-spaced binning, or a {help numlist} for manual specification of the positions
of inner knots (which must be within the range of {it:indvar}).
{p_end}

{p 4 8} {opt binsmethod(method)} specifies the method for data-driven selection of the number of bins via
the companion command {help binsregselect:binsregselect}.
The default is {cmd:binsmethod(dpi)}, which corresponds to the IMSE-optimal direct plug-in rule.
The other option is: {cmd:rot} for rule of thumb implementation.
{p_end}

{p 4 8} {opt nbinsrot(#)} specifies an initial number of bins value used to construct the DPI number of bins selector.
If not specified, the data-driven ROT selector is used instead.
{p_end}

{p 4 8} {opt randcut(#)} specifies the upper bound on a uniformly distributed variable used to draw a subsample 
for bins/degree/smoothness selection.
Observations for which {cmd:runiform()<=#} are used. # must be between 0 and 1. 
By default, max(5000, 0.01n) observations are used if the samples size n>5000.
{p_end}

{p 4 8} {opt pselect(numlist)} specifies a list of numbers within which the degree of polynomial {it:p} for 
point estimation is selected. If the selected optimal degree is {it:p}, then piecewise polynomials 
of degree {it:p+1} are used to conduct testing 
for nonparametric shape restrictions or parametric model specifications.
{p_end}

{p 4 8} {opt sselect(numlist)} specifies a list of numbers within which the number of smoothness constraints {it:s}
for point estimation.  If the selected optimal smoothness is {it:s}, 
then piecewise polynomials with {it:s+1} smoothness constraints are used to conduct testing 
for nonparametric shape restrictions or parametric model specifications.
If not specified, for each value {it:p} supplied in the 
option {cmd:pselect()}, only the piecewise polynomial with the maximum smoothness is considered, i.e., {it:s=p}. 
{p_end}

{p 4 8} Note: To implement the degree or smoothness selection, in addition to {cmd:pselect()} 
or {cmd:sselect()}, {cmd:nbins(#)} must be specified.
{p_end}

{dlgtab:Simulation}

{p 4 8} {opt nsims(#)} specifies the number of random draws for hypothesis testing.
The default is {cmd:nsims(500)}, which corresponds to 500 draws from a standard Gaussian random vector of size [(p+1)*J - (J-1)*s].
Setting at least {cmd:nsims(2000)} is recommended to obtain the final results.
{p_end}

{p 4 8} {opt simsgrid(#)} specifies the number of evaluation points of an evenly-spaced grid within each bin used
for evaluation of the supremum (infimum or Lp metric) operation needed for hypothesis testing procedures.
The default is {cmd:simsgrid(20)}, which corresponds to 20 evenly-spaced evaluation points within
each bin for approximating the supremum (infimum or Lp metric) operator.
Setting at least {cmd:simsgrid(50)} is recommended to obtain the final results.
{p_end}

{p 4 8} {opt simsseed(#)} sets the seed for simulations.
{p_end}

{dlgtab:Mass Points and Degrees of Freedom}

{p 4 8} {opt dfcheck(n1 n2)} sets cutoff values for minimum effective sample size checks, which take into account the number of unique values of {it:indvar}
(i.e., adjusting for the number of mass points), number of clusters, and degrees of freedom of the different statistical models considered.
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

{p 4 8} {cmd:vce(}{it:{help vcetype}}{cmd:)} specifies the {it:vcetype} for variance estimation used by the commands {help regress##options:regress},
{help logit##options:logit}, {help probit##options:probit},
{help qreg##qreg_options:qreg} or {cmd:reghdfe}. The default is {cmd:vce(robust)}.
{p_end}

{p 4 8} {opt asyvar(on/off)} specifies the method used to compute standard errors.
If {cmd:asyvar(on)} is specified, the standard error of the nonparametric component is used and the
uncertainty related to other control variables {it:othercovs} is omitted. Default is {cmd:asyvar(off)},
that is, the uncertainty related to {it:othercovs} is taken into account.
{p_end}

{p 4 8} {opt estmethodopt(cmd_option)} options to be passed on to the estimation command specified in {cmd:estmethod()}.
For example, options that control for the optimization process can be added here.
{p_end}

{p 4 8}{opt usegtools(on/off)} forces the use of several commands in the community-distributed Stata package {cmd:gtools}
to speed the computation up, if {it:on} is specified.
Default is {cmd:usegtools(off)}.
{p_end}

{p 4 8} For more information about the package {cmd:gtools}, please see {browse "https://gtools.readthedocs.io/en/latest/index.html":https://gtools.readthedocs.io/en/latest/index.html}.
{p_end}

{marker examples}{...}
{title:Examples}

{p 4 8} Setup{p_end}
{p 8 8} . {stata sysuse auto}

{p 4 8} Test for linearity{p_end}
{p 8 8} . {stata binstest mpg weight foreign, testmodelpoly(1)}{p_end}

{p 4 8} Test for monotonicity{p_end}
{p 8 8} . {stata binstest mpg weight foreign, deriv(1) bins(1 1) testshapel(0)}{p_end}


{marker stored_results}{...}
{title:Stored results}

{synoptset 17 tabbed}{...}
{p2col 5 17 21 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(Ndist)}}number of distinct values{p_end}
{synopt:{cmd:e(Nclust)}}number of clusters{p_end}
{synopt:{cmd:e(nbins)}}number of bins{p_end}
{synopt:{cmd:e(p)}}degree of polynomial for bin selection{p_end}
{synopt:{cmd:e(s)}}smoothness of polynomial for bin selection{p_end}
{synopt:{cmd:e(testshape_p)}}degree of polynomial for testing shape restrictions{p_end}
{synopt:{cmd:e(testshape_s)}}smoothness of polynomial for testing shape restrictions{p_end}
{synopt:{cmd:e(testmodel_p)}}degree of polynomial for testing model specifications{p_end}
{synopt:{cmd:e(testmodel_s)}}smoothness of polynomial for testing model specifications{p_end}
{synopt:{cmd:e(testpolyp)}}degree of polynomial regression model{p_end}
{synopt:{cmd:e(stat_poly)}}statistic for testing global polynomial model{p_end}
{synopt:{cmd:e(pval_poly)}}p value for testing global polynomial model{p_end}
{synopt:{cmd:e(imse_var_rot)}}variance constant in IMSE, ROT selection{p_end}
{synopt:{cmd:e(imse_bsq_rot)}}bias constant in IMSE, ROT selection{p_end}
{synopt:{cmd:e(imse_var_dpi)}}variance constant in IMSE, DPI selection{p_end}
{synopt:{cmd:e(imse_bsq_dpi)}}bias constant in IMSE, DPI selection{p_end}
{p2col 5 17 21 2: Macros}{p_end}
{synopt:{cmd:e(testvarlist)}}varlist found in {cmd:testmodel()}{p_end}
{synopt:{cmd:e(testvalue2)}}values in {cmd:testshape2()}{p_end}
{synopt:{cmd:e(testvalueR)}}values in {cmd:testshaper()}{p_end}
{synopt:{cmd:e(testvalueL)}}values in {cmd:testshapel()}{p_end}
{p2col 5 17 21 2: Matrices}{p_end}
{synopt:{cmd:e(pval_model)}}p values for {cmd:testmodel()}{p_end}
{synopt:{cmd:e(stat_model)}}statistics for {cmd:testmodel()}{p_end}
{synopt:{cmd:e(pval_shape2)}}p values for {cmd:testshape2()}{p_end}
{synopt:{cmd:e(stat_shape2)}}statistics for {cmd:testshape2()}{p_end}
{synopt:{cmd:e(pval_shapeR)}}p values for {cmd:testshaper()}{p_end}
{synopt:{cmd:e(stat_shapeR)}}statistics for {cmd:testshaper()}{p_end}
{synopt:{cmd:e(pval_shapeL)}}p values for {cmd:testshapel()}{p_end}
{synopt:{cmd:e(stat_shapeL)}}statistics for {cmd:testshapel()}{p_end}

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


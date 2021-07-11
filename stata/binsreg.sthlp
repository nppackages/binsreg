{smcl}
{* *! version 0.4.1 11-JUL-2021}{...}
{viewerjumpto "Syntax" "binsreg##syntax"}{...}
{viewerjumpto "Description" "binsreg##description"}{...}
{viewerjumpto "Options" "binsreg##options"}{...}
{viewerjumpto "Examples" "binsreg##examples"}{...}
{viewerjumpto "Stored results" "binsreg##stored_results"}{...}
{viewerjumpto "References" "binsreg##references"}{...}
{viewerjumpto "Authors" "binsreg##authors"}{...}
{cmd:help binsreg}
{hline}

{title:Title}

{p 4 8}{hi:binsreg} {hline 2} Data-Driven Binscatter Least Squares Estimation with Robust Inference Procedures and Plots.{p_end}


{marker syntax}{...}
{title:Syntax}

{p 4 12} {cmdab:binsreg} {depvar} {it:indvar} [{it:covars}] {ifin} {weight} [ {cmd:,} {opt deriv(v)} {opt at(position)}{p_end}
{p 12 12} {opt absorb(absvars)} {opt reghdfeopt(reghdfe_option)}{p_end}
{p 12 12} {opt dots(p s)} {opt dotsgrid(dotsgridoption)} {opt dotsplotopt(dotsoption)}{p_end}
{p 12 12} {opt line(p s)} {opt linegrid(#)} {opt lineplotopt(lineoption)}{p_end}
{p 12 12} {opt ci(p s)} {opt cigrid(cigridoption)} {opt ciplotopt(rcapoption)}{p_end}
{p 12 12} {opt cb(p s)} {opt cbgrid(#)} {opt cbplotopt(rareaoption)}{p_end}
{p 12 12} {opt polyreg(p)} {opt polyreggrid(#)} {opt polyregcigrid(#)} {opt polyregplotopt(lineoption)}{p_end}
{p 12 12} {opth by(varname)} {cmd:bycolors(}{it:{help colorstyle}list}{cmd:)} {cmd:bysymbols(}{it:{help symbolstyle}list}{cmd:)} {cmd:bylpatterns(}{it:{help linepatternstyle}list}{cmd:)}{p_end}
{p 12 12} {opt nbins(#)} {opt binspos(position)} {opt binsmethod(method)} {opt nbinsrot(#)} {opt samebinsby} {opt randcut(#)}{p_end}
{p 12 12} {opt nsims(#)} {opt simsgrid(#)} {opt simsseed(seed)}{p_end}
{p 12 12} {opt dfcheck(n1 n2)} {opt masspoints(masspointsoption)}{p_end}
{p 12 12} {cmd:vce(}{it:{help vcetype}}{cmd:)} {opt asyvar(on/off)}{p_end}
{p 12 12} {opt level(level)} {opt usegtools(on/off)} {opt noplot} {opt savedata(filename)} {opt replace}{p_end}
{p 12 12} {opt plotxrange(min max)} {opt plotyrange(min max)} {it:{help twoway_options}} ]{p_end}

{p 4 8} where {depvar} is the dependent variable, {it:indvar} is the independent variable for binning, and {it:covars} are other covariates to be controlled for.{p_end}

{p 4 8} p, s and v are integers satisfying 0 <= s,v <= p, which can take different values in each case.{p_end}

{p 4 8} {opt fweight}s, {opt aweight}s and {opt pweight}s are allowed; see {help weight}.{p_end}

{marker description}{...}
{title:Description}

{p 4 8} {cmd:binsreg} implements binscatter least squares estimation with robust inference procedure and plots, following the results in
{browse "https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2019_Binscatter.pdf":Cattaneo, Crump, Farrell and Feng (2021a)}.
Binscatter provides a flexible way of describing the mean relationship between two variables, after possibly adjusting for other covariates, based on partitioning/binning of the independent variable of interest.
The main purpose of this command is to generate binned scatter plots with curve estimation with robust pointwise confidence intervals and uniform confidence band.
If the binning scheme is not set by the user, the companion command {help binsregselect:binsregselect} is used to implement binscatter in a data-driven (optimal) way.
Hypothesis testing about the regression function can be conducted via the companion command {help binstest:binstest}. Hypothesis testing about pairwise group comparison can be conducted via the companion command {help binspwc: binspwc}.
{p_end}

{p 4 8} A detailed introduction to this command is given in
{browse "https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2019_Stata.pdf":Cattaneo, Crump, Farrell and Feng (2021b)}.
A companion R package with the same capabilities is available (see website below).
{p_end}

{p 4 8} Companion commands: {help binstest:binstest} for hypothesis testing, and {help binsregselect:binsregselect} data-driven (optimal) binning selection.{p_end}

{p 4 8} Related Stata and R packages are available in the following website:{p_end}

{p 8 8} {browse "https://nppackages.github.io/":https://nppackages.github.io/}{p_end}


{marker options}{...}
{title:Options}

{dlgtab:Estimand}

{p 4 8} {opt deriv(v)} specifies the derivative order of the regression function for estimation and plotting.
The default is {cmd:deriv(0)}, which corresponds to the function itself.
{p_end}

{p 4 8} {opt at(position)} specifies the values of {it:covars} at which the estimated function is evaluated for plotting.
The default is {cmd:at(mean)}, which corresponds to the mean of {it:covars}. Other options are: {cmd:at(median)} for the median of {it:covars}, {cmd:at(0)} for zeros,
and {cmd:at(filename)} for particular values of {it:covars} saved in another file.

{p 4 8} Note: when {cmd:at(mean)} or {cmd:at(median)} is specified, all factor variables in {it:covars} (if specified) are excluded from the evaluation.
{p_end}

{dlgtab:Reghdfe}

{p 4 8} {opt absorb(absvars)} specifies categorical variables (or interactions) representing the fixed effects to be absorbed. This is equivalent to including an indicator/dummy variable for each category of each {it:absvar}.
When {cmd:absorb()} is specified, the community-contributed command {cmd:reghdfe} instead of the command {cmd:regress} is used.
{p_end}

{p 4 8} {opt reghdfeopt(reghdfe_option)} options to be passed on to the command {cmd:reghdfe}.
{p_end}

{p 4 8} {it:Important:} 

{p 6 9} 1. Fixed effects added via {cmd:absorb()} are not used in the evaluation of the estimated function, regardless of the option specified within {cmd:at()}.
To plot the binscatter function for a particular category of interest, save the values of {it:covars} at which the function is evaluated in another file, say,
{cmd:wval}, specify the corresponding factor variables in {it:covar} directly, and add the option {cmd:at(wval)}.
{p_end}

{p 6 9} 2. {cmd:absorb()} and {cmd:vce()} should not be specified within {cmd:reghdfeopt()}.
{p_end}

{p 6 9} 3. Make sure the package reghdfe installed has a version number greater than or equal to 5.9.0 (03jun2020). An older version may result in an error in Mata.
{p_end}

{p 4 8} For more information about the community-contributed command {cmd:reghdfe}, please see {browse "http://scorreia.com/software/reghdfe/":http://scorreia.com/software/reghdfe/}.
{p_end}

{dlgtab:Dots}

{p 4 8} {opt dots(p s)} sets a piecewise polynomial of degree {it:p} with {it:s} smoothness constraints for point estimation and plotting as "dots".
The default is {cmd:dots(0 0)}, which corresponds to piecewise constant (canonical binscatter).
{p_end}

{p 4 8} {opt dotsgrid(dotsgridoption)} specifies the number and location of dots within each bin to be plotted.
Two options are available: {it:mean} and a {it:numeric} non-negative integer.
The option {opt dotsgrid(mean)} adds the sample average of {it:indvar} within each bin to the grid of evaluation points.
The option {opt dotsgrid(#)} adds {it:#} number of evenly-spaced points to the grid of evaluation points for each bin.
Both options can be used simultaneously: for example, {opt dotsgrid(mean 5)} generates six evaluation points within each bin containing the sample mean of {it:indvar} within each bin and five evenly-spaced points.
Given this choice, the dots are point estimates evaluated over the selected grid within each bin.
The default is {opt dotsgrid(mean)}, which corresponds to one dot per bin evaluated at the sample average of {it:indvar} within each bin (canonical binscatter).
{p_end}

{p 4 8} {opt dotsplotopt(dotsoption)} standard graphs options to be passed on to the {help twoway:twoway} command to modify the appearance of the plotted dots.
{p_end}

{dlgtab:Line}
  
{p 4 8} {opt line(p s)} sets a piecewise polynomial of degree {it:p} with {it:s} smoothness constraints for plotting as a "line".
By default, the line is not included in the plot unless explicitly specified.
Recommended specification is {cmd:line(3 3)}, which adds a cubic B-spline estimate of the regression function of interest to the binned scatter plot. 
{p_end}

{p 4 8} {opt linegrid(#)} specifies the number of evaluation points of an evenly-spaced grid within each bin used for evaluation of the point estimate set by the {cmd:line(p s)} option.
The default is {cmd:linegrid(20)}, which corresponds to 20 evenly-spaced evaluation points within each bin for fitting/plotting the line.
{p_end}

{p 4 8} {opt lineplotopt(lineoption)} standard graphs options to be passed on to the {help twoway:twoway} command to modify the appearance of the plotted line.
{p_end}

{dlgtab:Confidence Intervals}

{p 4 8} {opt ci(p s)} specifies the piecewise polynomial of degree {it:p} with {it:s} smoothness constraints used for constructing confidence intervals.
By default, the confidence intervals are not included in the plot unless explicitly specified.
Recommended specification is {cmd:ci(3 3)}, which adds confidence intervals based on a cubic B-spline estimate of the regression function of interest to the binned scatter plot.
{p_end}

{p 4 8} {opt cigrid(cigridoption)} specifies the number and location of evaluation points in the grid used to construct the confidence intervals set by the {opt ci(p s)} option.
Two options are available: {it:mean} and a {it:numeric} non-negative integer.
The option {opt cigrid(mean)} adds the sample average of {it:indvar} within each bin to the grid of evaluation points.
The option {opt cigrid(#)} adds {it:#} number of evenly-spaced points to the grid of evaluation points for each bin.
Both options can be used simultaneously: for example, {opt cigrid(mean 5)} generates six evaluation points within each bin containing the sample mean of {it:indvar} within each bin and five evenly-spaced points.
The default is {opt cigrid(mean)}, which corresponds to one evaluation point set at the sample average of {it:indvar} within each bin for confidence interval construction.
{p_end}

{p 4 8} {opt ciplotopt(rcapoption)} standard graphs options to be passed on to the {help twoway:twoway} command to modify the appearance of the confidence intervals.
{p_end}

{dlgtab:Confidence Band}

{p 4 8} {opt cb(p s)} specifies the piecewise polynomial of degree {it:p} with {it:s} smoothness constraints used for constructing the confidence band.
By default, the confidence band is not included in the plot unless explicitly specified.
Recommended specification is {cmd:cb(3 3)}, which adds a confidence band based on a cubic B-spline estimate of the regression function of interest to the binned scatter plot.
{p_end}

{p 4 8} {opt cbgrid(#)} specifies the number of evaluation points of an evenly-spaced grid within each bin used for evaluation of the point estimate set by the {cmd:cb(p s)} option.
The default is {cmd:cbgrid(20)}, which corresponds to 20 evenly-spaced evaluation points within each bin for confidence band construction.
{p_end}

{p 4 8} {opt cbplotopt(rareaoption)} standard graphs options to be passed on to the {help twoway:twoway} command to modify the appearance of the confidence band.
{p_end}

{dlgtab:Global Polynomial Regression}

{p 4 8} {opt polyreg(p)} sets the degree {it:p} of a global polynomial regression model for plotting.
By default, this fit is not included in the plot unless explicitly specified.
Recommended specification is {cmd:polyreg(3)}, which adds a fourth order global polynomial fit of the regression function of interest to the binned scatter plot. 
{p_end}

{p 4 8} {opt polyreggrid(#)} specifies the number of evaluation points of an evenly-spaced grid within each bin used for evaluation of the point estimate set by the {cmd:polyreg(p)} option.
The default is {cmd:polyreggrid(20)}, which corresponds to 20 evenly-spaced evaluation points within each bin for confidence interval construction.
{p_end}

{p 4 8} {opt polyregcigrid(#)} specifies the number of evaluation points of an evenly-spaced grid within each bin used for constructing confidence intervals based on polynomial regression set by the {cmd:polyreg(p)} option.
The default is {cmd:polyregcigrid(0)}, which corresponds to not plotting confidence intervals for the global polynomial regression approximation.
{p_end}

{p 4 8} {opt polyregplotopt(lineoption)} standard graphs options to be passed on to the {help twoway:twoway} command to modify the appearance of the global polynomial regression fit.
{p_end}

{dlgtab:Subgroup Analysis}

{p 4 8} {opt by(varname)} specifies the variable containing the group indicator to perform subgroup analysis; both numeric and string variables are supported.
When {opt by(varname)} is specified, {cmdab:binsreg} implements estimation and inference by each subgroup separately, but produces a common binned scatter plot.
By default, the binning structure is selected for each subgroup separately, but see the option {cmd:samebinsby} below for imposing a common binning structure across subgroups.
{p_end}

{p 4 8} {cmd:bycolors(}{it:{help colorstyle}list}{cmd:)} specifies an ordered list of colors for plotting each subgroup series defined by the option {opt by()}.
{p_end}

{p 4 8} {cmd:bysymbols(}{it:{help symbolstyle}list}{cmd:)} specifies an ordered list of symbols for plotting each subgroup series defined by the option {opt by()}.
{p_end}

{p 4 8} {cmd:bylpatterns(}{it:{help linepatternstyle}list}{cmd:)} specifies an ordered list of line patterns for plotting each subgroup series defined by the option {opt by()}.
{p_end}

{dlgtab:Partitioning/Binning Selection}

{p 4 8} {opt nbins(#)} sets the number of bins for partitioning/binning of {it:indvar}.
If not specified, the number of bins is selected via the companion command {help binsregselect:binsregselect} in a data-driven, optimal way whenever possible.
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

{p 4 8} {opt samebinsby} forces a common partitioning/binning structure across all subgroups specified by the option {cmd:by()}. The knots positions are selected according to the option {cmd:binspos()} and using the full sample.
If {cmd:nbins()} is not specified, then the number of bins is selected via the companion command {help binsregselect:binsregselect} and using the full sample.
{p_end}

{p 4 8} {opt randcut(#)} specifies the upper bound on a uniformly distributed variable used to draw a subsample for bins selection.
Observations for which {cmd:runiform()<=#} are used. # must be between 0 and 1.
{p_end}

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

{dlgtab:Standard Error}

{p 4 8} {cmd:vce(}{it:{help vcetype}}{cmd:)} specifies the {it:vcetype} for variance estimation used by the command {help regress##options:regress}.
The default is {cmd:vce(robust)}.
{p_end}

{p 4 8} {opt asyvar(on/off)} specifies the method used to compute standard errors.
If {cmd:asyvar(on)} is specified, the standard error of the nonparametric component is used and the uncertainty related to other control variables {it:covars} is omitted.
Default is {cmd:asyvar(off)}, that is, the uncertainty related to {it:covars} is taken into account.
{p_end}

{dlgtab:Other Options}

{p 4 8} {opt level(#)} sets the nominal confidence level for confidence interval and confidence band estimation.
{p_end}

{p 4 8}{opt usegtools(on/off)} forces the use of several commands in the community-distributed Stata package {cmd:gtools} to speed the computation up, if {it:on} is specified.
Default is {cmd:usegtools(off)}.
{p_end}

{p 4 8} For more information about the package {cmd:gtools}, please see {browse "https://gtools.readthedocs.io/en/latest/index.html":https://gtools.readthedocs.io/en/latest/index.html}.
{p_end}

{p 4 8} {opt noplot} omits binscatter plotting.
{p_end}

{p 4 8} {opt savedata(filename)} specifies a filename for saving all data underlying the binscatter plot (and more).
{p_end}

{p 4 8} {opt replace} overwrites the existing file when saving the graph data.
{p_end}

{p 4 8} {opt plotxrange(min max)} specifies the range of the x-axis for plotting. Observations outside the range are dropped in the plot.
{p_end}

{p 4 8} {opt plotyrange(min max)} specifies the range of the y-axis for plotting. Observations outside the range are dropped in the plot.
{p_end}

{p 4 8} {it:{help twoway_options}} any unrecognized options are appended to the end of the twoway command generating the binned scatter plot.
{p_end}


{marker examples}{...}
{title:Examples}

{p 4 8} Setup{p_end}
{p 8 8} . {stata  sysuse auto}{p_end}

{p 4 8} Run a binscatter regression and report the plot{p_end}
{p 8 8} . {stata binsreg mpg weight foreign}{p_end}

{p 4 8} Add confidence intervals and confidence band{p_end}
{p 8 8} . {stata binsreg mpg weight foreign, ci(3 3) cb(3 3)}{p_end}

{p 4 8} Run binscatter regression by group{p_end}
{p 8 8} . {stata binsreg mpg weight, by(foreign)}{p_end}

{marker stored_results}{...}
{title:Stored results}

{synoptset 17 tabbed}{...}
{p2col 5 17 21 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(level)}}confidence level{p_end}
{synopt:{cmd:e(dots_p)}}degree of polynomial for dots{p_end}
{synopt:{cmd:e(dots_s)}}smoothness of polynomial for dots{p_end}
{synopt:{cmd:e(line_p)}}degree of polynomial for line{p_end}
{synopt:{cmd:e(line_s)}}smoothness of polynomial for line{p_end}
{synopt:{cmd:e(ci_p)}}degree of polynomial for confidence interval{p_end}
{synopt:{cmd:e(ci_s)}}smoothness of polynomial for confidence interval{p_end}
{synopt:{cmd:e(cb_p)}}degree of polynomial for confidence band{p_end}
{synopt:{cmd:e(cb_s)}}smoothness of polynomial for confidence band{p_end}
{p2col 5 17 21 2: Matrices}{p_end}
{synopt:{cmd:e(N_by)}}number of observations for each group{p_end}
{synopt:{cmd:e(Ndist_by)}}number of distinct values for each group{p_end}
{synopt:{cmd:e(Nclust_by)}}number of clusters for each group{p_end}
{synopt:{cmd:e(nbins_by)}}number of bins for each group{p_end}
{synopt:{cmd:e(cval_by)}}critical value for each group, used for confidence bands{p_end}

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

{p 4 8} Yingjie Feng, Tsinghua University, Beijing, China.
{browse "mailto:fengyingjiepku@gmail.com":fengyingjiepku@gmail.com}.
{p_end}


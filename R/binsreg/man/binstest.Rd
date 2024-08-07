% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binstest.R
\name{binstest}
\alias{binstest}
\title{Data-Driven Nonparametric Shape Restriction and Parametric Model Specification Testing using Binscatter}
\usage{
binstest(y, x, w = NULL, data = NULL, estmethod = "reg",
  family = gaussian(), quantile = NULL, deriv = 0, at = NULL,
  nolink = F, testmodel = NULL, testmodelparfit = NULL,
  testmodelpoly = NULL, testshape = NULL, testshapel = NULL,
  testshaper = NULL, testshape2 = NULL, lp = Inf, bins = NULL,
  nbins = NULL, pselect = NULL, sselect = NULL, binspos = "qs",
  binsmethod = "dpi", nbinsrot = NULL, randcut = NULL, nsims = 500,
  simsgrid = 20, simsseed = NULL, vce = NULL, cluster = NULL,
  asyvar = F, dfcheck = c(20, 30), masspoints = "on", weights = NULL,
  subset = NULL, numdist = NULL, numclust = NULL, estmethodopt = NULL,
  ...)
}
\arguments{
\item{y}{outcome variable. A vector.}

\item{x}{independent variable of interest. A vector.}

\item{w}{control variables. A matrix, a vector or a \code{\link{formula}}.}

\item{data}{an optional data frame containing variables used in the model.}

\item{estmethod}{estimation method. The default is \code{estmethod="reg"} for tests based on binscatter least squares regression. Other options are \code{"qreg"} for quantile regression and \code{"glm"} for generalized linear regression. If \code{estmethod="glm"}, the option \code{family} must be specified.}

\item{family}{a description of the error distribution and link function to be used in the generalized linear model when \code{estmethod="glm"}. (See \code{\link{family}} for details of family functions.)}

\item{quantile}{the quantile to be estimated. A number strictly between 0 and 1.}

\item{deriv}{derivative order of the regression function for estimation, testing and plotting.
The default is \code{deriv=0}, which corresponds to the function itself.}

\item{at}{value of \code{w} at which the estimated function is evaluated.  The default is \code{at="mean"}, which corresponds to
the mean of \code{w}. Other options are: \code{at="median"} for the median of \code{w}, \code{at="zero"} for a vector of zeros.
\code{at} can also be a vector of the same length as the number of columns of \code{w} (if \code{w} is a matrix) or a data frame containing the same variables as specified in \code{w} (when
\code{data} is specified). Note that when \code{at="mean"} or \code{at="median"}, all factor variables (if specified) are excluded from the evaluation (set as zero).}

\item{nolink}{if true, the function within the inverse link function is reported instead of the conditional mean function for the outcome.}

\item{testmodel}{a vector or a logical value. It sets the degree of polynomial and the number of smoothness constraints for parametric model specification
testing. If \code{testmodel=c(p,s)} is specified, a piecewise polynomial of degree \code{p} with \code{s} smoothness constraints is used.
If \code{testmodel=T} or \code{testmodel=NULL} (default) is specified, \code{testmodel=c(1,1)} is used unless the degree \code{p} or the smoothness \code{s}
selection is requested via the option \code{pselect} or \code{sselect} (see more details in the explanation of \code{pselect} and \code{sselect}).}

\item{testmodelparfit}{a data frame or matrix which contains the evaluation grid and fitted values of the model(s) to be
tested against.  The column contains a series of evaluation points
at which the binscatter model and the parametric model of interest are compared with
each other.  Each parametric model is represented by other columns, which must
contain the fitted values at the corresponding evaluation points.}

\item{testmodelpoly}{degree of a global polynomial model to be tested against.}

\item{testshape}{a vector or a logical value. It sets the degree of polynomial and the number of smoothness constraints for nonparametric shape restriction
testing. If \code{testshape=c(p,s)} is specified, a piecewise polynomial of degree \code{p} with \code{s} smoothness constraints is used.
If \code{testshape=T} or \code{testshape=NULL} (default) is specified, \code{testshape=c(1,1)} is used unless the degree \code{p} or smoothness \code{s} selection
is requested via the option \code{pselect} or \code{sselect} (see more details in the explanation of \code{pselect} and \code{sselect}).}

\item{testshapel}{a vector of null boundary values for hypothesis testing. Each number \code{a} in the vector
corresponds to one boundary of a one-sided hypothesis test to the left of the form
\code{H0: sup_x mu(x)<=a}.}

\item{testshaper}{a vector of null boundary values for hypothesis testing. Each number \code{a} in the vector
corresponds to one boundary of a one-sided hypothesis test to the right of the form
\code{H0: inf_x mu(x)>=a}.}

\item{testshape2}{a vector of null boundary values for hypothesis testing. Each number \code{a} in the vector
corresponds to one boundary of a two-sided hypothesis test of the form
\code{H0: sup_x |mu(x)-a|=0}.}

\item{lp}{an Lp metric used for parametric model specification testing and/or shape restriction testing. The default is \code{lp=Inf}, which
corresponds to the sup-norm of the t-statistic. Other options are \code{lp=q} for a positive number \code{q>=1}. Note that
\code{lp=Inf} ("sup-norm") has to be used for testing one-sided shape restrictions.}

\item{bins}{a vector. If \code{bins=c(p,s)}, it sets the piecewise polynomial of degree \code{p} with \code{s} smoothness constraints
for data-driven (IMSE-optimal) selection of the partitioning/binning scheme.
The default is \code{bins=c(0,0)}, which corresponds to the piecewise constant.}

\item{nbins}{number of bins for partitioning/binning of \code{x}.  If \code{nbins=T} or \code{nbins=NULL} (default) is specified, the number
of bins is selected via the companion command \code{\link{binsregselect}} in a data-driven, optimal way whenever possible.
If a vector with more than one number is specified, the number of bins is selected within this vector via the companion command \code{\link{binsregselect}}.}

\item{pselect}{vector of numbers within which the degree of polynomial \code{p} for point estimation is selected.
 If the selected optimal degree is \code{p}, then piecewise polynomials of degree
 \code{p+1} are used to conduct testing for nonparametric shape restrictions or parametric model specifications.
 \emph{Note:} To implement the degree or smoothness selection, in addition to \code{pselect} or \code{sselect},
\code{nbins=#} must be specified.}

\item{sselect}{vector of numbers within which the number of smoothness constraints \code{s} for point estimation is selected.
If the selected optimal smoothness is \code{s}, then piecewise polynomials of \code{s+1} smoothness constraints
are used to conduct testing for nonparametric shape restrictions or parametric model specifications.
If not specified, for each value \code{p} supplied in the option \code{pselect}, only the piecewise polynomial
with the maximum smoothness is considered, i.e., \code{s=p}.}

\item{binspos}{position of binning knots. The default is \code{binspos="qs"}, which corresponds to quantile-spaced
binning (canonical binscatter).  The other options are \code{"es"} for evenly-spaced binning, or
a vector for manual specification of the positions of inner knots (which must be within the range of
\code{x}).}

\item{binsmethod}{method for data-driven selection of the number of bins. The default is \code{binsmethod="dpi"},
which corresponds to the IMSE-optimal direct plug-in rule.  The other option is: \code{"rot"}
for rule of thumb implementation.}

\item{nbinsrot}{initial number of bins value used to construct the DPI number of bins selector.
If not specified, the data-driven ROT selector is used instead.}

\item{randcut}{upper bound on a uniformly distributed variable used to draw a subsample for bins/degree/smoothness selection.
Observations for which \code{runif()<=#} are used. # must be between 0 and 1. By default, \code{max(5000, 0.01n)} observations
are used if the samples size \code{n>5000}.}

\item{nsims}{number of random draws for hypothesis testing. The default is
\code{nsims=500}, which corresponds to 500 draws from a standard Gaussian random vector of size
\code{[(p+1)*J - (J-1)*s]}. Setting at least \code{nsims=2000} is recommended to obtain the final results.}

\item{simsgrid}{number of evaluation points of an evenly-spaced grid within each bin used for evaluation of
the supremum (infimum or Lp metric) operation needed to construct hypothesis testing
procedures. The default is \code{simsgrid=20}, which corresponds to 20 evenly-spaced
evaluation points within each bin for approximating the supremum (infimum or Lp metric) operator.
Setting at least \code{simsgrid=50} is recommended to obtain the final results.}

\item{simsseed}{seed for simulation.}

\item{vce}{procedure to compute the variance-covariance matrix estimator. For least squares regression and generalized linear regression, the allowed options are the same as that for \code{\link{binsreg}} or \code{\link{binsqreg}}.
For quantile regression, the allowed options are the same as that for \code{\link{binsqreg}}.}

\item{cluster}{cluster ID. Used for compute cluster-robust standard errors.}

\item{asyvar}{if true, the standard error of the nonparametric component is computed and the uncertainty related to control
variables is omitted. Default is \code{asyvar=FALSE}, that is, the uncertainty related to control variables is taken into account.}

\item{dfcheck}{adjustments for minimum effective sample size checks, which take into account number of unique
values of \code{x} (i.e., number of mass points), number of clusters, and degrees of freedom of
the different stat models considered. The default is \code{dfcheck=c(20, 30)}.
See \href{https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2024_Stata.pdf}{Cattaneo, Crump, Farrell and Feng (2024c)} for more details.}

\item{masspoints}{how mass points in \code{x} are handled. Available options:
\itemize{
\item \code{"on"} all mass point and degrees of freedom checks are implemented. Default.
\item \code{"noadjust"} mass point checks and the corresponding effective sample size adjustments are omitted.
\item \code{"nolocalcheck"} within-bin mass point and degrees of freedom checks are omitted.
\item \code{"off"} "noadjust" and "nolocalcheck" are set simultaneously.
\item \code{"veryfew"} forces the function to proceed as if \code{x} has only a few number of mass points (i.e., distinct values).
                       In other words, forces the function to proceed as if the mass point and degrees of freedom checks were failed.
}}

\item{weights}{an optional vector of weights to be used in the fitting process. Should be \code{NULL} or
a numeric vector. For more details, see \code{\link{lm}}.}

\item{subset}{optional rule specifying a subset of observations to be used.}

\item{numdist}{number of distinct values for selection. Used to speed up computation.}

\item{numclust}{number of clusters for selection. Used to speed up computation.}

\item{estmethodopt}{a list of optional arguments used by \code{\link[quantreg]{rq}} (for quantile regression) or \code{\link{glm}} (for fitting generalized linear models).}

\item{...}{optional arguments to control bootstrapping if \code{estmethod="qreg"} and \code{vce="boot"}. See \code{\link[quantreg]{boot.rq}}.}
}
\value{
\item{\code{testshapeL}}{Results for \code{testshapel}, including: \code{testvalL}, null boundary values;
                                \code{stat.shapeL}, test statistics; and \code{pval.shapeL}, p-value.}
       \item{\code{testshapeR}}{Results for \code{testshaper}, including: \code{testvalR}, null boundary values;
                                \code{stat.shapeR}, test statistics; and \code{pval.shapeR}, p-value.}
       \item{\code{testshape2}}{Results for \code{testshape2}, including: \code{testval2}, null boundary values;
                                \code{stat.shape2}, test statistics; and \code{pval.shape2}, p-value.}
       \item{\code{testpoly}}{Results for \code{testmodelpoly}, including: \code{testpoly}, the degree of global polynomial;
                              \code{stat.poly}, test statistic; \code{pval.poly}, p-value.}
       \item{\code{testmodel}}{Results for \code{testmodelparfit}, including: \code{stat.model}, test statistics;
                               \code{pval.model}, p-values.}
       \item{\code{imse.var.rot}}{Variance constant in IMSE, ROT selection.}
       \item{\code{imse.bsq.rot}}{Bias constant in IMSE, ROT selection.}
       \item{\code{imse.var.dpi}}{Variance constant in IMSE, DPI selection.}
       \item{\code{imse.bsq.dpi}}{Bias constant in IMSE, DPI selection.}
       \item{\code{opt}}{ A list containing options passed to the function, as well as total sample size \code{n},
                          number of distinct values \code{Ndist} in \code{x}, number of clusters \code{Nclust}, and
                          number of bins \code{nbins}.}
}
\description{
\code{binstest} implements binscatter-based hypothesis testing procedures for parametric functional
            forms of and nonparametric shape restrictions on the regression function of interest, following the results
            in \href{https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2024_AER.pdf}{Cattaneo, Crump, Farrell and Feng (2024a)} and
            \href{https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2024_NonlinearBinscatter.pdf}{Cattaneo, Crump, Farrell and Feng (2024b)}.
            If the binning scheme is not set by the user,
            the companion function \code{\link{binsregselect}} is used to implement binscatter in a
            data-driven way and inference procedures are based on robust bias correction.
            Binned scatter plots based on different methods can be constructed using the companion functions \code{\link{binsreg}},
            \code{\link{binsqreg}} or \code{\link{binsglm}}.
}
\examples{
 x <- runif(500); y <- sin(x)+rnorm(500)
 est <- binstest(y,x, testmodelpoly=1)
 summary(est)
}
\references{
Cattaneo, M. D., R. K. Crump, M. H. Farrell, and Y. Feng. 2024a: \href{https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2024_AER.pdf}{On Binscatter}. American Economic Review 114(5): 1488-1514.

Cattaneo, M. D., R. K. Crump, M. H. Farrell, and Y. Feng. 2024b: \href{https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2024_NonlinearBinscatter.pdf}{Nonlinear Binscatter Methods}. Working Paper.

Cattaneo, M. D., R. K. Crump, M. H. Farrell, and Y. Feng. 2024c: \href{https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2024_Stata.pdf}{Binscatter Regressions}. Working Paper.
}
\seealso{
\code{\link{binsreg}}, \code{\link{binsqreg}}, \code{\link{binsglm}}, \code{\link{binsregselect}}.
}
\author{
Matias D. Cattaneo, Princeton University, Princeton, NJ. \email{cattaneo@princeton.edu}.

Richard K. Crump, Federal Reserve Bank of New York, New York, NY. \email{richard.crump@ny.frb.org}.

Max H. Farrell, UC Santa Barbara, Santa Barbara, CA. \email{mhfarrell@gmail.com}.

Yingjie Feng (maintainer), Tsinghua University, Beijing, China. \email{fengyingjiepku@gmail.com}.
}

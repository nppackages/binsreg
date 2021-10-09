################################################################################################
#'@title  Data-Driven Binscatter Quantile Regression with Robust Inference Procedures and Plots
#'@description \code{binsqreg} implements binscatter quantile regression with robust inference procedures and plots, following the
#'             results in \href{https://arxiv.org/abs/1902.09608}{Cattaneo, Crump, Farrell and Feng (2021a)}.
#'             Binscatter provides a flexible way to describe the quantile relationship between two variables, after
#'             possibly adjusting for other covariates, based on partitioning/binning of the independent variable of interest.
#'             The main purpose of this function is to generate binned scatter plots with curve estimation with robust pointwise confidence intervals and
#'             uniform confidence band.  If the binning scheme is not set by the user, the companion function
#'             \code{\link{binsregselect}} is used to implement binscatter in a data-driven way. Hypothesis testing about the function of interest
#'             can be conducted via the companion function \code{\link{binstest}}.
#'@param  y outcome variable. A vector.
#'@param  x independent variable of interest. A vector.
#'@param  w control variables. A matrix, a vector or a \code{\link{formula}}.
#'@param  data an optional data frame containing variables in the model.
#'@param  at value of \code{w} at which the estimated function is evaluated.  The default is \code{at="mean"}, which corresponds to
#'           the mean of \code{w}. Other options are: \code{at="median"} for the median of \code{w}, \code{at="zero"} for a vector of zeros.
#'           \code{at} can also be a vector of the same length as the number of columns of \code{w} (if \code{w} is a matrix) or a data frame containing the same variables as specified in \code{w} (when
#'           \code{data} is specified). Note that when \code{at="mean"} or \code{at="median"}, all factor variables (if specified) are excluded from the evaluation (set as zero).
#'@param  quantile the quantile to be estimated. A number strictly between 0 and 1.
#'@param  deriv  derivative order of the regression function for estimation, testing and plotting.
#'               The default is \code{deriv=0}, which corresponds to the function itself.
#'@param  dots a vector. \code{dots=c(p,s)} sets a piecewise polynomial of degree \code{p} with \code{s} smoothness constraints for
#'             point estimation and plotting as "dots". The default is \code{dots=c(0,0)}, which corresponds to
#'             piecewise constant (canonical binscatter)
#'@param  dotsgrid number of dots within each bin to be plotted. Given the choice, these dots are point estimates
#'                 evaluated over an evenly-spaced grid within each bin. The default is \code{dotsgrid=0}, and only
#'                 the point estimates at the mean of \code{x} within each bin are presented.
#'@param  dotsgridmean If true, the dots corresponding to the point estimates evaluated at the mean of \code{x} within each bin
#'                     are presented. By default, they are presented, i.e., \code{dotsgridmean=T}.
#'@param  line a vector. \code{line=c(p,s)} sets a piecewise polynomial of degree \code{p} with \code{s} smoothness constraints
#'             for plotting as a "line".  By default, the line is not included in the plot unless explicitly
#'             specified.  Recommended specification is \code{line=c(3,3)}, which adds a cubic B-spline estimate
#'             of the regression function of interest to the binned scatter plot.
#'@param  linegrid number of evaluation points of an evenly-spaced grid within each bin used for evaluation of
#'                 the point estimate set by the \code{line=c(p,s)} option. The default is \code{linegrid=20},
#'                 which corresponds to 20 evenly-spaced evaluation points within each bin for fitting/plotting the line.
#'@param  ci a vector. \code{ci=c(p,s)} sets a piecewise polynomial of degree \code{p} with \code{s} smoothness constraints used for
#'             constructing confidence intervals. By default, the confidence intervals are not included in the plot
#'             unless explicitly specified.  Recommended specification is \code{ci=c(3,3)}, which adds confidence
#'             intervals based on cubic B-spline estimate of the regression function of interest to the binned scatter plot.
#'@param  cigrid number of evaluation points of an evenly-spaced grid within each bin used for evaluation of the point
#'               estimate set by the \code{ci=c(p,s)} option. The default is \code{cigrid=1}, which corresponds to 1
#'               evenly-spaced evaluation point within each bin for confidence interval construction.
#'@param  cigridmean If true, the confidence intervals corresponding to the point estimates evaluated at the mean of \code{x} within each bin
#'                   are presented. The default is \code{cigridmean=T}.
#'@param  cb a vector. \code{cb=c(p,s)} sets a the piecewise polynomial of degree \code{p} with \code{s} smoothness constraints used for
#'           constructing the confidence band. By default, the confidence band is not included in the plot unless
#'           explicitly specified.  Recommended specification is \code{cb=c(3,3)}, which adds a confidence band
#'           based on cubic B-spline estimate of the regression function of interest to the binned scatter plot.
#'@param  cbgrid number of evaluation points of an evenly-spaced grid within each bin used for evaluation of the point
#'               estimate set by the \code{cb=c(p,s)} option. The default is \code{cbgrid=20}, which corresponds
#'               to 20 evenly-spaced evaluation points within each bin for confidence interval construction.
#'@param  polyreg degree of a global polynomial regression model for plotting. By default, this fit is not included
#'                in the plot unless explicitly specified. Recommended specification is \code{polyreg=3}, which
#'                adds a cubic (global) polynomial fit of the regression function of interest to the binned scatter plot.
#'@param  polyreggrid number of evaluation points of an evenly-spaced grid within each bin used for evaluation of
#'                    the point estimate set by the \code{polyreg=p} option. The default is \code{polyreggrid=20},
#'                    which corresponds to 20 evenly-spaced evaluation points within each bin for confidence
#'                    interval construction.
#'@param  polyregcigrid number of evaluation points of an evenly-spaced grid within each bin used for constructing
#'                      confidence intervals based on polynomial regression set by the \code{polyreg=p} option.
#'                      The default is \code{polyregcigrid=0}, which corresponds to not plotting confidence
#'                      intervals for the global polynomial regression approximation.
#'@param  by a vector containing the group indicator for subgroup analysis; both numeric and string variables
#'           are supported. When \code{by} is specified, \code{binsreg} implements estimation and inference for each subgroup
#'           separately, but produces a common binned scatter plot. By default, the binning structure is selected for each
#'           subgroup separately, but see the option \code{samebinsby} below for imposing a common binning structure across subgroups.
#'@param  bycolors an ordered list of colors for plotting each subgroup series defined by the option \code{by}.
#'@param  bysymbols an ordered list of symbols for plotting each subgroup series defined by the option \code{by}.
#'@param  bylpatterns an ordered list of line patterns for plotting each subgroup series defined by the option \code{by}.
#'@param  legendTitle String, title of legend.
#'@param  legendoff If true, no legend is added.
#'@param  nbins number of bins for partitioning/binning of \code{x}.  If not specified, the number of bins is
#'              selected via the companion function \code{binsregselect} in a data-driven, optimal way whenever possible.
#'@param  binspos position of binning knots. The default is \code{binspos="qs"}, which corresponds to quantile-spaced
#'                binning (canonical binscatter).  The other options are \code{"es"} for evenly-spaced binning, or
#'                a vector for manual specification of the positions of inner knots (which must be within the range of
#'                \code{x}).
#'@param  binsmethod method for data-driven selection of the number of bins. The default is \code{binsmethod="dpi"},
#'                   which corresponds to the IMSE-optimal direct plug-in rule.  The other option is: \code{"rot"}
#'                   for rule of thumb implementation.
#'@param  nbinsrot initial number of bins value used to construct the DPI number of bins selector.
#'                 If not specified, the data-driven ROT selector is used instead.
#'@param  samebinsby if true, a common partitioning/binning structure across all subgroups specified by the option \code{by} is forced.
#'                   The knots positions are selected according to the option \code{binspos} and using the full sample. If \code{nbins}
#'                   is not specified, then the number of bins is selected via the companion command \code{\link{binsregselect}} and
#'                   using the full sample.
#'@param  randcut upper bound on a uniformly distributed variable used to draw a subsample for bins selection.
#'                Observations for which \code{runif()<=#} are used. # must be between 0 and 1.
#'@param  nsims number of random draws for constructing confidence bands. The default is
#'              \code{nsims=500}, which corresponds to 500 draws from a standard Gaussian random vector of size
#'              \code{[(p+1)*J - (J-1)*s]}.
#'@param  simsgrid number of evaluation points of an evenly-spaced grid within each bin used for evaluation of
#'                 the supremum operation needed to construct confidence bands.
#'                 The default is \code{simsgrid=20}, which corresponds to 20 evenly-spaced
#'                 evaluation points within each bin for approximating the supremum operator.
#'@param  simsseed  seed for simulation.
#'@param  vce Procedure to compute the variance-covariance matrix estimator (see \code{\link{summary.rq}} for more details). Options are
#'           \itemize{
#'           \item \code{"iid"} which presumes that the errors are iid and computes an estimate of the asymptotic covariance matrix as in KB(1978).
#'           \item \code{"nid"} which presumes local (in quantile) linearity of the the conditional quantile functions and computes a Huber sandwich estimate using a local estimate of the sparsity.
#'           \item \code{"ker"} which uses a kernel estimate of the sandwich as proposed by Powell (1991).
#'           \item \code{"boot"} which implements one of several possible bootstrapping alternatives for estimating standard errors including a variate of the wild bootstrap for clustered response. See \code{\link{boot.rq}} for further details.
#'           }
#'@param  cluster cluster ID. Used for compute cluster-robust standard errors.
#'@param  asyvar  If true, the standard error of the nonparametric component is computed and the uncertainty related to control
#'                variables is omitted. Default is \code{asyvar=FALSE}, that is, the uncertainty related to control variables is taken into account.
#'@param  level nominal confidence level for confidence interval and confidence band estimation. Default is \code{level=95}.
#'@param  noplot If true, no plot produced.
#'@param  dfcheck adjustments for minimum effective sample size checks, which take into account number of unique
#'                values of \code{x} (i.e., number of mass points), number of clusters, and degrees of freedom of
#'                the different statistical models considered. The default is \code{dfcheck=c(20, 30)}.
#'                See \href{https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2021_Stata.pdf}{Cattaneo, Crump, Farrell and Feng (2021b)} for more details.
#'@param  masspoints how mass points in \code{x} are handled. Available options:
#'                   \itemize{
#'                   \item \code{"on"} all mass point and degrees of freedom checks are implemented. Default.
#'                   \item \code{"noadjust"} mass point checks and the corresponding effective sample size adjustments are omitted.
#'                   \item \code{"nolocalcheck"} within-bin mass point and degrees of freedom checks are omitted.
#'                   \item \code{"off"} "noadjust" and "nolocalcheck" are set simultaneously.
#'                   \item \code{"veryfew"} forces the function to proceed as if \code{x} has only a few number of mass points (i.e., distinct values).
#'                                          In other words, forces the function to proceed as if the mass point and degrees of freedom checks were failed.
#'                   }
#'@param  weights an optional vector of weights to be used in the fitting process. Should be \code{NULL} or
#'                a numeric vector. For more details, see \code{\link{lm}}.
#'@param  subset  optional rule specifying a subset of observations to be used.
#'@param  plotxrange a vector. \code{plotxrange=c(min, max)} specifies a range of the x-axis for plotting. Observations outside the range are dropped in the plot.
#'@param  plotyrange a vector. \code{plotyrange=c(min, max)} specifies a range of the y-axis for plotting. Observations outside the range are dropped in the plot.
#'@param  ...     optional arguments to control bootstrapping. See \code{\link{boot.rq}}.
#'@return \item{\code{bins_plot}}{A \code{ggplot} object for binscatter plot.}
#'        \item{\code{data.plot}}{A list containing data for plotting. Each item is a sublist of data frames for each group. Each sublist may contain the following data frames:
#'        \itemize{
#'        \item \code{data.dots} Data for dots. It contains: \code{x}, evaluation points; \code{bin}, the indicator of bins;
#'                               \code{isknot}, indicator of inner knots; \code{mid}, midpoint of each bin; and \code{fit}, fitted values.
#'        \item \code{data.line} Data for line. It contains: \code{x}, evaluation points; \code{bin}, the indicator of bins;
#'                                \code{isknot}, indicator of inner knots; \code{mid}, midpoint of each bin; and \code{fit}, fitted values.
#'        \item \code{data.ci} Data for CI. It contains: \code{x}, evaluation points; \code{bin}, the indicator of bins;
#'                                \code{isknot}, indicator of inner knots; \code{mid}, midpoint of each bin;
#'                                \code{ci.l} and \code{ci.r}, left and right boundaries of each confidence intervals.
#'        \item \code{data.cb} Data for CB. It contains: \code{x}, evaluation points; \code{bin}, the indicator of bins;
#'                                \code{isknot}, indicator of inner knots; \code{mid}, midpoint of each bin;
#'                                \code{cb.l} and \code{cb.r}, left and right boundaries of the confidence band.
#'        \item \code{data.poly} Data for polynomial regression. It contains: \code{x}, evaluation points;
#'                                \code{bin}, the indicator of bins;
#'                                \code{isknot}, indicator of inner knots; \code{mid}, midpoint of each bin; and
#'                                \code{fit}, fitted values.
#'        \item \code{data.polyci} Data for confidence intervals based on polynomial regression. It contains: \code{x}, evaluation points;
#'                                \code{bin}, the indicator of bins;
#'                                \code{isknot}, indicator of inner knots; \code{mid}, midpoint of each bin;
#'                                \code{polyci.l} and \code{polyci.r}, left and right boundaries of each confidence intervals.}}
#'        \item{\code{cval.by}}{A vector of critical values for constructing confidence band for each group.}
#'        \item{\code{opt}}{ A list containing options passed to the function, as well as \code{N.by} (total sample size for each group),
#'                           \code{Ndist.by} (number of distinct values in \code{x} for each group), \code{Nclust.by} (number of clusters for each group),
#'                           and \code{nbins.by} (number of bins for each group), and \code{byvals} (number of distinct values in \code{by}).}
#'
#'@author
#' Matias D. Cattaneo, Princeton University, Princeton, NJ. \email{cattaneo@princeton.edu}.
#'
#' Richard K. Crump, Federal Reserve Bank of New York, New York, NY. \email{richard.crump@ny.frb.org}.
#'
#' Max H. Farrell, University of Chicago, Chicago, IL. \email{max.farrell@chicagobooth.edu}.
#'
#' Yingjie Feng (maintainer), Tsinghua University, Beijing, China. \email{fengyingjiepku@gmail.com}.
#'
#'@references
#' Cattaneo, M. D., R. K. Crump, M. H. Farrell, and Y. Feng. 2021a: \href{https://arxiv.org/abs/1902.09608}{On Binscatter}. Working Paper.
#'
#' Cattaneo, M. D., R. K. Crump, M. H. Farrell, and Y. Feng. 2021b: \href{https://arxiv.org/abs/1902.09615}{Binscatter Regressions}. Working Paper.
#'
#'@seealso \code{\link{binsregselect}}, \code{\link{binstest}}.
#'
#'@examples
#'  x <- runif(500); y <- sin(x)+rnorm(500)
#'  ## Binned scatterplot
#'  binsqreg(y,x)
#'@export

binsqreg <- function(y, x, w=NULL, data=NULL, at=NULL, quantile=0.5, deriv=0,
                     dots=c(0,0), dotsgrid=0, dotsgridmean=T, line=NULL, linegrid=20,
                     ci=NULL, cigrid=0, cigridmean=T, cb=NULL, cbgrid=20,
                     polyreg=NULL, polyreggrid=20, polyregcigrid=0,
                     by=NULL, bycolors=NULL, bysymbols=NULL, bylpatterns=NULL,
                     legendTitle=NULL, legendoff=F,
                     nbins=NULL, binspos="qs", binsmethod="dpi", nbinsrot=NULL, samebinsby=F, randcut=NULL,
                     nsims=500, simsgrid=20, simsseed=NULL,
                     vce="nid",cluster=NULL, asyvar=F, level=95,
                     noplot=F, dfcheck=c(20,30), masspoints="on", weights=NULL, subset=NULL, plotxrange=NULL, plotyrange=NULL, ...) {

  # param for internal use
  qrot <- 2

  ####################
  ### prepare data ###
  ####################

  # variable name
  xname <- deparse(substitute(x))
  yname <- deparse(substitute(y))
  byname <- deparse(substitute(by))
  weightsname <- deparse(substitute(weights))
  subsetname  <- deparse(substitute(subset))
  clustername <- deparse(substitute(cluster))

  # extract y, x, w, weights, subset, if needed (w as a matrix, others as vectors)
  # generate design matrix for covariates W
  w.factor <- NULL
  if (is.null(data)) {
    if (is.data.frame(y)) y <- y[,1]
    if (is.data.frame(x)) x <- x[,1]
    if (!is.null(w)) {
      if (is.vector(w)) {
        w <- as.matrix(w)
        if (is.character(w)) {
          w <- model.matrix(~w)[,-1,drop=F]
          w.factor <- rep(T, ncol(w))
        }
      } else if (is.formula(w)) {
        w.model <- binsreg.model.mat(w)
        w <- w.model$design
        w.factor <- w.model$factor.colnum
      }
    }
  } else {
    if (xname %in% names(data)) x <- data[,xname]
    if (yname %in% names(data)) y <- data[,yname]
    if (byname != "NULL") if (byname %in% names(data)) {
      by <- data[,byname]
    }
    if (weightsname != "NULL") if (weightsname %in% names(data)) {
      weights <- data[,weightsname]
    }
    if (subsetname != "NULL") if (subsetname %in% names(data))  {
      subset  <- data[,subsetname]
    }
    if (clustername != "NULL") if (clustername %in% names(data)) {
      cluster <- data[,clustername]
    }
    if (deparse(substitute(w))!="NULL") {
      if (is.formula(w)) {
        if (is.data.frame(at)) {
          if (ncol(data)!=ncol(at)) {
            data[setdiff(names(at), names(data))] <- NA
            at[setdiff(names(data), names(at))] <- NA
          }
          data <- rbind(data, at)
        }
        w.model <- binsreg.model.mat(w, data=data)
        w <- w.model$design
        w.factor <- w.model$factor.colnum
        if (is.data.frame(at)) {
          eval.w <- w[nrow(w),]
          w <- w[-nrow(w),, drop=F]
        }
      } else {
        if (deparse(substitute(w)) %in% names(data)) w <- data[,deparse(substitute(w))]
        w <- as.matrix(w)
        if (is.character(w)) {
          w <- model.matrix(~w)[,-1,drop=F]
          w.factor <- rep(T, ncol(w))
        }
      }
    }
  }

  if (!is.null(w)) nwvar <- ncol(w)
  else             nwvar <- 0

  # extract subset
  if (!is.null(subset)) {
    y <- y[subset]
    x <- x[subset]
    w <- w[subset, , drop = F]
    by <- by[subset]
    if (!is.null(cluster)) cluster <- cluster[subset]
    if (!is.null(weights)) weights <- weights[subset]
  }

  na.ok <- complete.cases(x) & complete.cases(y)
  if (!is.null(w))  na.ok <- na.ok & complete.cases(w)
  if (!is.null(by)) na.ok <- na.ok & complete.cases(by)
  if (!is.null(weights)) na.ok <- na.ok & complete.cases(weights)
  if (!is.null(cluster)) na.ok <- na.ok & complete.cases(cluster)

  y  <- y[na.ok]
  x  <- x[na.ok]
  w  <- w[na.ok, , drop = F]
  by <- by[na.ok]
  weights <- weights[na.ok]
  cluster <- cluster[na.ok]

  xmin <- min(x); xmax <- max(x)

  # evaluation point of w
  if (is.null(at))  at <- "mean"

  if (vce=="iid") vce.select <- "const"
  else            vce.select <- "HC1"

  ##################################################
  ######## Error Checking ##########################
  ##################################################
  exit <- 0
  if (deriv < 0) {
    print("derivative incorrectly specified.")
    exit <- 1
  }
  if (dotsgrid<0|linegrid<0|cigrid<0|cbgrid<0|polyreggrid<0|polyregcigrid<0) {
    print("# of evaluation points incorrectly specified.")
    exit <- 1
  }
  if (!is.null(nbins)) {
    if (nbins < 0) {
      print("# of bins incorrectly specified.")
      exit <- 1
    }
  }
  if (!is.character(binspos)) {
    if (min(binspos)<=xmin|max(binspos)>=xmax) {
      print("knots out of allowed range")
      exit <- 1
    }
  } else {
    if (binspos!="es" & binspos!="qs") {
      print("binspos incorrectly specified.")
      exit <- 1
    }
  }
  if (length(dots)==2) if (dots[1] < dots[2]) {
    print("p<s not allowed.")
    exit <- 1
  }
  if (!is.null(line)) if (length(line)==2) if (line[1]<line[2]) {
    print("p<s not allowed.")
    exit <- 1
  }
  if (!is.null(ci)) if (length(ci)==2) if (ci[1]<ci[2]) {
    print("p<s not allowed.")
    exit <- 1
  }
  if (!is.null(cb)) if (length(cb)==2) if (cb[1]<cb[2]) {
    print("p<s not allowed.")
    exit <- 1
  }
  if (dots[1] < deriv) {
    print("p<deriv not allowed.")
    exit <- 1
  }
  if (!is.null(line)) if (line[1] < deriv) {
    print("p<deriv not allowed.")
    exit <- 1
  }
  if (!is.null(ci)) if (ci[1] < deriv) {
    print("p<deriv not allowed.")
    exit <- 1
  }
  if (!is.null(cb)) if (cb[1] < deriv) {
    print("p<deriv not allowed.")
    exit <- 1
  }
  if (binsmethod!="dpi" & binsmethod!="rot") {
    print("bin selection method incorrectly specified.")
    exit <- 1
  }
  if (!is.null(w)) {
    if (is.character(at)) {
      if (at!="mean" & at!="median" & at!="zero") {
        print("at incorrectly specified.")
        exit <- 1
      }
    } else {
      if (is.vector(at)) {
        if (length(at)!=nwvar) {
          print("length of at not equal to # of w variables.")
          exit <- 1
        }
      } else {
        if (!is.data.frame(at)) {
          print("at incorrectly specified.")
          exit <- 1
        }
      }
    }
  }
  if (exit > 0) {
    stop()
  }

  ##################################################
  # Prepare options
  dots.p <- dots[1]; dots.s <- dots[2]
  dotsmean <- 0
  if (dotsgridmean) dotsmean <- 1

  if (is.null(line)) {
    linegrid <- 0
  } else {
    line.p <- line[1]
    if (length(line)==1) {
      line.s <- line.p
    } else {
      line.s <- line[2]
    }
  }

  cimean <- 0
  if (cigridmean) cimean <- 1
  if (is.null(ci)) {
    cigrid <- 0
    cimean <- 0
  } else {
    ci.p <- ci[1]
    if (length(ci)==1) {
      ci.s <- ci.p
    } else {
      ci.s <- ci[2]
    }
  }

  if (is.null(cb)) {
    cbgrid <- 0
  } else {
    cb.p <- cb[1]
    if (length(cb)==1) {
      cb.s <- cb.p
    } else {
      cb.s <- cb[2]
    }
  }

  if (is.null(polyreg)) {
    polyreggrid   <- 0
    polyregcigrid <- 0
  }


  localcheck <- massadj <- T; fewmasspoints <- F
  if (masspoints=="on") {
    localcheck <- T; massadj <- T
  } else if (masspoints=="off") {
    localcheck <- F; massadj <- F
  } else if (masspoints=="noadjust") {
    localcheck <- T; massadj <- F
  } else if (masspoints=="nolocalcheck") {
    localcheck <- F; massadj <- T
  } else if (masspoints=="veryfew") {
    fewmasspoints <- T
  }

  #############################
  # Extract byvals in by ######
  if (!is.null(by)) {
    byvals <- sort(unique(by))
    ngroup <- length(byvals)
  } else {
    byvals <- "Full Sample"
    ngroup <- 1
  }

  ########################################
  # Default plotting options
  if (length(bycolors)==0) {
    bycolors <- c("navy", "maroon", "forestgreen", "darkorange", "lavenderblush3",
                  "khaki", "sienna", "steelblue", "brown", "gold", "gray45")
    bycolors <- rep(bycolors, length.out=ngroup)
  } else {
    bycolors <- rep(bycolors, length.out=ngroup)
  }

  if (length(bysymbols)==0) {
    bysymbols <- c(19, 15:18, 0:14)
    bysymbols <- rep(bysymbols, length.out=ngroup)
  } else {
    bysymbols <- rep(bysymbols, length.out=ngroup)
  }

  if (length(bylpatterns)==0) {
    bylpatterns <- rep(1, ngroup)
  } else {
    bylpatterns <- rep(bylpatterns, length.out=ngroup)
  }

  # Legend
  if (ngroup==1) {
    legendoff <- T    # turn off legend when only one group is available
  } else {
    legendGroups <- paste(byvals)
    if (length(legendTitle) == 0) {
      legendTitle <- paste(byname, "=")
    } else {
      legendTitle <- legendTitle[1]
    }
  }

  #########################################
  if (binsmethod=="dpi") {
    selectmethod <- "IMSE direct plug-in"
  } else {
    selectmethod <- "IMSE rule-of-thumb"
  }
  nbins_all <- nbins         # "nbins" is reserved for use within loop
  if (!is.null(nbins)) {
    selectmethod <- "User-specified"
  }

  knot <- NULL
  knotlistON <- F
  if (!is.character(binspos)) {
    nbins <- length(binspos)+1
    knot <- c(xmin, sort(binspos), xmax)
    position <- "User-specified"
    es <- F
    selectmethod <- "User-specified"
    knotlistON <- T
  } else {
    if (binspos == "es") {
      es <- T
      position <- "Evenly-spaced"
    } else {
      es <- F
      position <- "Quantile-spaced"
    }
  }
  Ndist <- Nclust <- NA

  ################################################
  ### Bin selection using full sample if needed ##
  fullfewobs <- fewobs <- selectfullON <- F
  if (fewmasspoints) fullfewobs <- fewobs <- T
  if (!fullfewobs & is.null(nbins) & (is.null(by) | (!is.null(by) & samebinsby))) selectfullON <- T
  if (selectfullON) {
    # effective size
    eN <- N <- length(x)

    if (massadj) {
      Ndist <- length(unique(x))
      eN <- min(eN, Ndist)
    }

    if (!is.null(cluster)) {
      Nclust <- length(unique(cluster))
      eN <- min(eN, Nclust)
    }

    # check if rot can be implemented
    if (is.null(nbinsrot)) {
      if (eN <= dfcheck[1]+dots.p+1+qrot) {
        warning("too small effective sample size for bin selection. # of mass of points or clusters used and by option ignored.")
        fewobs <- fullfewobs <- T
        byvals <- "Full Sample"
        es <- F
        position <- "Quantile-spaced"
      }
    }
    if (!fewobs) {
      if (is.na(Ndist))  {
        Ndist.sel <- NULL
      } else {
        Ndist.sel <- Ndist
      }
      if (is.na(Nclust)) {
        Nclust.sel <- NULL
      } else {
        Nclust.sel <- Nclust
      }
      binselect <- binsregselect(y, x, w, deriv=deriv,
                                 bins=dots, binspos=binspos,
                                 binsmethod=binsmethod, nbinsrot=nbinsrot,
                                 vce=vce.select, cluster=cluster, randcut=randcut,
                                 dfcheck=dfcheck, masspoints=masspoints, weights=weights,
                                 numdist=Ndist.sel, numclust=Nclust.sel)
      if (is.na(binselect$nbinsrot.regul)) {
        print("bin selection fails.")
        stop()
      }
      if (binsmethod == "rot") {
        nbins <- binselect$nbinsrot.regul
      } else if (binsmethod == "dpi") {
        nbins <- binselect$nbinsdpi
        if (is.na(nbins)) {
          warning("DPI selection fails. ROT choice used.")
          nbins <- binselect$nbinsrot.regul
        }
      }
    }
  }

  # Generate knot using the full sample if needed

  if ((selectfullON | (!is.null(nbins) & samebinsby)) & !fullfewobs & is.null(knot)) {
    knotlistON <- T
    if (es) {
      knot <- genKnot.es(xmin, xmax, nbins)
    } else {
      knot <- genKnot.qs(x, nbins)
    }
  }

  knot_all <- NULL
  if (knotlistON) {
     knot_all <- knot    # universal knot sequence
  }

  ##########################################
  if (!is.null(simsseed)) set.seed(simsseed)
  alpha <- (100-(100 - level)/2)/100

  ##################################################################
  N.by <- Ndist.by <- Nclust.by <- nbins.by <- cval.by <- NULL   # save results
  data.plot <- list()     # list storing graph data

  ##################################################################
  ##################### ENTER the loop #############################
  ##################################################################
  for (i in 1:ngroup) {
    # Take subsample
    if (ngroup != 1) {
      sub <- (by == byvals[i])
      y.sub <- y[sub]; x.sub <- x[sub]; w.sub <- w[sub, , drop = F]
      cluster.sub <- cluster[sub]; weights.sub <- weights[sub]
    } else {
      y.sub <- y; x.sub <- x; w.sub <- w; cluster.sub <- cluster
      weights.sub <- weights
    }

    # Effective size
    xmin <- min(x.sub); xmax <- max(x.sub)
    eN <- N <- length(x.sub)
    N.by <- c(N.by, N)

    Ndist <- NA
    if (massadj) {
      Ndist <- length(unique(x.sub))
      eN <- min(eN, Ndist)
    }
    Ndist.by <- c(Ndist.by, Ndist)

    Nclust <- NA
    if (!is.null(cluster.sub)) {
      Nclust <- length(unique(cluster.sub))
      eN <- min(eN, Nclust)
    }
    Nclust.by <- c(Nclust.by, Nclust)

    #######################################################
    ############### Bin selection if needed ###############
    nbins <- NULL; knot <- NULL                       # initialize again
    if (is.null(nbins_all) & !knotlistON & !fullfewobs) {
      # check if rot can be implemented
      if (is.null(nbinsrot)) if (eN <= dfcheck[1]+dots.p+1+qrot) {
          warning("too small effective sample size for bin selection.
                  # of mass points or clusters used.")
          fewobs <- T
          nbins <- eN
          es <- F
      }
      if (!fewobs) {
        if (is.na(Ndist))  {
          Ndist.sel <- NULL
        } else {
          Ndist.sel <- Ndist
        }
        if (is.na(Nclust)) {
          Nclust.sel <- NULL
        } else {
          Nclust.sel <- Nclust
        }
        binselect <- binsregselect(y.sub, x.sub, w.sub, deriv=deriv,
                                   bins=dots, binspos=binspos,
                                   binsmethod=binsmethod, nbinsrot=nbinsrot,
                                   vce=vce.select, cluster=cluster.sub, randcut=randcut,
                                   dfcheck=dfcheck, masspoints=masspoints, weights=weights.sub,
                                   numdist=Ndist.sel, numclust=Nclust.sel)
        if (is.na(binselect$nbinsrot.regul)) {
          print("bin selection fails.")
          stop()
        }
        if (binsmethod == "rot") {
          nbins <- binselect$nbinsrot.regul
        } else if (binsmethod == "dpi") {
          nbins <- binselect$nbinsdpi
          if (is.na(nbins)) {
            warning("DPI selection fails. ROT choice used.")
            nbins <- binselect$nbinsrot.regul
          }
        }
      }
    }
    if (!is.null(nbins_all)) nbins <- nbins_all
    if (knotlistON) {
      nbins <- length(knot_all)-1
      knot  <- knot_all
    }
    if (fullfewobs) {
      fewobs <- T
      nbins <- eN
    }

    ###########################################
    # Checking for each case
    dots.fewobs <- line.fewobs <- ci.fewobs <- cb.fewobs <- polyreg.fewobs <- F
    if (!fewobs) {
      if ((nbins-1)*(dots.p-dots.s+1)+dots.p+1+dfcheck[2]>=eN) {
        fewobs <- T
        nbins <- eN
        es <- F
        warning("too small effective sample size for dots. # of mass points or clusters used.")
      }
      if (!is.null(line)) {
        if ((nbins-1)*(line.p-line.s+1)+line.p+1+dfcheck[2]>=eN) {
          line.fewobs <- T
          warning("too small effective sample size for line.")
        }
      }
      if (!is.null(ci)) {
        if ((nbins-1)*(ci.p-ci.s+1)+ci.p+1+dfcheck[2]>=eN) {
          ci.fewobs <- T
          warning("too small effective sample size for ci.")
        }
      }
      if (!is.null(cb)) {
        if ((nbins-1)*(cb.p-cb.s+1)+cb.p+1+dfcheck[2]>=eN) {
          cb.fewobs <- T
          warning("too small effective sample size for line.")
        }
      }
    }
    if (!is.null(polyreg)) if (polyreg+1>eN) {
       polyreg.fewobs <- T
       warning("too small effective sample size for polynomial fit.")
    }

    ####################################
    ########### Generate knot ##########
    ####################################
    fewmass <- F
    if (fewobs) if (!is.na(Ndist)) if (eN==Ndist) fewmass <- T
    if (fewmasspoints) fewmass <- T

    if (is.null(knot)) {
      if (!fewmass) {
          if (es) knot <- genKnot.es(xmin, xmax, nbins)
          else    knot <- genKnot.qs(x.sub, nbins)
      }
    } else {
      if (fewobs) if (!is.na(Ndist)) if (eN!=Ndist) {
         knot <- genKnot.qs(x.sub, nbins)
      }
    }

    # knot for few mass points
    if (fewmass) {
      knot <- sort(unique(x.sub))
      if (fewmasspoints) {
        eN <- nbins <- Ndist <- length(knot)
      }
    } else {
      knot <- c(knot[1], unique(knot[-1]))
      if (nbins!=length(knot)-1) {
        warning("repeated knots. Some bins dropped.")
        nbins <- length(knot)-1
      }
    }

    # check local mass points
    if (!fewobs & localcheck) {
      uniqmin <- binsreg.checklocalmass(x.sub, nbins, es, knot=knot) # mimic STATA
      if (!is.null(dots)) {
        if (uniqmin < dots.p+1) {
          dots.fewobs <- T
          warning("some bins have too few distinct values of x for dots.")
        }
      }
      if (!is.null(line)) {
        if (uniqmin < line.p+1) {
          line.fewobs <- T
          warning("some bins have too few distinct values of x for line.")
        }
      }
      if (!is.null(ci)) {
        if (uniqmin < ci.p+1) {
          ci.fewobs <- T
          warning("some bins have too few distinct values of x for CI.")
        }
      }
      if (!is.null(cb)) {
        if (uniqmin < cb.p+1) {
          cb.fewobs <- T
          warning("some bins have too few distinct values of x for CB.")
        }
      }
    }

    # NOW, save nbins
    nbins.by <- c(nbins.by, nbins)

    #######################################
    ###### Prepare Plots ##################
    #######################################
    data.by <- list()    # initialize data list
    # data.dots <- data.line <- data.ci<- data.cb <- data.poly <- data.polyci <- NULL

    # adjust w variables
    if (!is.null(w.sub)) {
      if (is.character(at)) {
        if (at=="mean") {
          eval.w <- colWeightedMeans(x=w.sub, w=weights.sub)
          if (!is.null(w.factor)) eval.w[w.factor] <- 0
        } else if (at=="median") {
          eval.w <- colWeightedMedians(x=w.sub, w=weights.sub)
          if (!is.null(w.factor)) eval.w[w.factor] <- 0
        } else if  (at=="zero") {
          eval.w <- rep(0, nwvar)
        }
      } else if (is.vector(at)) {
        eval.w <- at
      } else if (is.data.frame(at)) {
        eval.w <- eval.w
      }
    } else {
      eval.w <- NULL
    }

    ############################################
    # Dots and CI for Small eN case
    ############################################

    if (dotsmean+dotsgrid!=0 & !noplot & fewobs) {
      warning("dots=c(0,0) used.")
      k <- nbins
      if (!is.na(Ndist)) {
        if (eN == Ndist) {
          dots.x <- knot
          # RENEW knot, each value forms a category
          xcat.few  <- findInterval(x.sub, knot)
        } else {
          dots.x <- (knot[-1]+knot[-length(knot)])/2
          xcat.few <- findInterval(x.sub, knot, rightmost.closed = T, left.open = T)
        }
      } else {
        dots.x <- (knot[-1]+knot[-length(knot)])/2
        xcat.few  <- findInterval(x.sub, knot, rightmost.closed = T, left.open = T)
      }

      if (is.null(w.sub)) {
        model <- rq(y.sub~-1+factor(xcat.few), tau=quantile, weights=weights.sub)
      } else {
        model <- rq(y.sub~-1+factor(xcat.few)+w.sub, tau=quantile, weights=weights.sub)
      }
      beta  <- model$coeff[1:k]
      beta[is.na(beta)] <- 0
      vcv   <- binsreg.vcov(model, type=vce, cluster=cluster.sub, is.qreg=TRUE, ...)

      dots.fit <- beta
      if (!is.null(eval.w)) {
        coeff.w <- model$coeff[-(1:k)]
        coeff.w[is.na(coeff.w)] <- 0
        dots.fit <- dots.fit + sum(eval.w * coeff.w)
      }
      data.dots <- data.frame(group=paste(byvals[i]), x=dots.x, fit=dots.fit)
      data.by$data.dots <- data.dots

      if (cigrid+cimean!=0) {
        warning("ci=c(0,0) used.")
        basis.all <- cbind(diag(length(dots.x)), outer(rep(1, length(dots.x)), eval.w))
        dots.se  <- sqrt(rowSums((basis.all %*% vcv) * basis.all))
        ci.l <- dots.fit - dots.se * qnorm(alpha)
        ci.r <- dots.fit + dots.se * qnorm(alpha)
        data.ci <- data.frame(group=byvals[i], x=dots.x, ci.l=ci.l, ci.r=ci.r)
        data.by$data.ci <- data.ci
      }
    }

    ##########################################
    ########## Usual Case ####################
    ##########################################

    ###############################################
    dotsON <- lineON <- ciON <- cbON <- polyON <- F
    xmean <- NULL

    ################ Dots ####################
    if (dotsmean+dotsgrid !=0 & !noplot & !dots.fewobs & !fewobs) {
      dotsON <- T
    }
    if (dotsON) {
      dots.x <- dots.bin <- dots.isknot <- dots.mid <- c()
      if (dotsmean!=0) {
        xcat   <- findInterval(x.sub, knot, rightmost.closed = T, left.open = T)
        xmean  <- as.vector(tapply(x.sub, xcat, FUN=mean))
        dots.x <- c(dots.x, xmean)
        dots.bin <- c(dots.bin, 1:nbins)
        dots.isknot <- c(dots.isknot, rep(0, nbins))
        dots.mid <- c(dots.mid, rep(0, nbins))
      }
      if (dotsgrid!=0) {
        grid <- binsreg.grid(knot, dotsgrid, addmore=T)
        dots.x <- c(dots.x, grid$eval); dots.bin <- c(dots.bin, grid$bin)
        dots.isknot <- c(dots.isknot, grid$isknot); dots.mid <- c(dots.mid, grid$mid)
      }
      B    <- binsreg.spdes(eval=x.sub, p=dots.p, s=dots.s, deriv=0, knot=knot)
      P    <- cbind(B, w.sub)         # full design matrix
      model.dots <- rq(y.sub~-1+P, tau=quantile, weights=weights.sub)
      check.drop(model.dots$coeff, ncol(B))

      basis <- binsreg.spdes(eval=dots.x, p=dots.p, s=dots.s, knot=knot, deriv=deriv)
      pred.dots  <- binsreg.pred(basis, model.dots, type = "xb", vce=vce, cluster=cluster.sub, deriv=deriv, wvec=eval.w)
      dots.fit <- pred.dots$fit
      dots.fit[dots.isknot==1] <- NA
      data.dots <- data.frame(group=paste(byvals[i]), x=dots.x, bin=dots.bin, isknot=dots.isknot,
                              mid=dots.mid, fit=dots.fit)
      data.by$data.dots <- data.dots
    }

    ################ Line ####################
    if (linegrid !=0 & !noplot & !line.fewobs & !fewobs) {
      lineON <- T
    }
    if (lineON) {
      grid <- binsreg.grid(knot, linegrid, addmore=T)
      line.x <- grid$eval; line.bin <- grid$bin
      line.isknot <- grid$isknot; line.mid <- grid$mid

      line.reg.ON <- T
      if (dotsON) if(line.p==dots.p & line.s==dots.s) {
        model.line <- model.dots; line.reg.ON <- F
      }
      if (line.reg.ON) {
        B     <- binsreg.spdes(eval=x.sub, p=line.p, s=line.s, deriv=0, knot=knot)
        P     <- cbind(B, w.sub)         # full design matrix
        model.line <- rq(y.sub~-1+P, tau=quantile, weights=weights.sub)
        check.drop(model.line$coeff, ncol(B))
      }

      basis <- binsreg.spdes(eval=line.x, p=line.p, s=line.s, knot=knot, deriv=deriv)
      pred.line  <- binsreg.pred(basis, model.line, type = "xb", vce=vce, cluster=cluster.sub, deriv=deriv, wvec=eval.w)
      line.fit <- pred.line$fit
      if (line.s == 0 | line.s-deriv <= 0) {
        line.fit[line.isknot==1] <- NA
      }
      data.line <- data.frame(group=byvals[i], x=line.x, bin=line.bin, isknot=line.isknot,
                              mid=line.mid, fit=line.fit)
      data.by$data.line <- data.line
    }

    ############### Poly fit #########################
    if (polyreggrid!=0 & !noplot & !polyreg.fewobs) {
      polyON <- T
    }
    if (polyON) {
      grid <- binsreg.grid(knot, polyreggrid, addmore=T)
      poly.x <- grid$eval; poly.bin <- grid$bin
      poly.isknot <- grid$isknot; poly.mid <- grid$mid

      # Run a poly reg
      x.p <- matrix(NA, N, polyreg+1)
      for (j in 1:(polyreg+1))  x.p[,j] <- x.sub^(j-1)
      P.poly <- cbind(x.p, w.sub)
      model.poly <- rq(y.sub~-1+P.poly, tau=quantile, weights=weights.sub)
      beta.poly <- model.poly$coefficients
      beta.poly[is.na(beta.poly)] <- 0
      poly.fit  <- 0
      for (j in deriv:polyreg) {
        poly.fit <- poly.fit + poly.x^(j-deriv)*beta.poly[j+1]*factorial(j)/factorial(j-deriv)
      }
      if (!is.null(eval.w) & deriv==0) poly.fit <- poly.fit + sum(beta.poly[-(1:(polyreg+1))]*eval.w)

      data.poly <- data.frame(group=byvals[i], x=poly.x, bin=poly.bin, isknot=poly.isknot,
                              mid=poly.mid, fit=poly.fit)
      data.by$data.poly <- data.poly

      # add CI?
      if (polyregcigrid!=0) {
        grid <- binsreg.grid(knot, polyregcigrid, addmore=T)
        polyci.x <- grid$eval; polyci.bin <- grid$bin
        polyci.isknot <- grid$isknot; polyci.mid <- grid$mid

        npolyci.x <- length(polyci.x)
        basis.polyci <- matrix(NA, npolyci.x, polyreg+1)
        for (j in 1:(polyreg+1))  {
          if (j-1>=deriv) {
            basis.polyci[,j] <- polyci.x^(j-deriv-1)*factorial(j-1)/factorial(j-1-deriv)
          } else {
            basis.polyci[,j] <- rep(0, npolyci.x)
          }
        }
        if (!is.null(eval.w)) {
          if (deriv==0) basis.polyci <- cbind(basis.polyci, outer(rep(1, nrow(basis.polyci)), eval.w))
          else          basis.polyci <- cbind(basis.polyci, outer(rep(1, nrow(basis.polyci)), rep(0, nwvar)))
        }

        polyci.pred <- binsreg.pred(basis.polyci, model=model.poly, type="all",
                                    vce=vce, cluster=cluster.sub, is.qreg=TRUE, avar=T, ...)
        polyci.l <- polyci.pred$fit - qnorm(alpha) * polyci.pred$se
        polyci.r <- polyci.pred$fit + qnorm(alpha) * polyci.pred$se

        data.polyci <- data.frame(group=byvals[i], x=polyci.x, bin=polyci.bin, isknot=polyci.isknot,
                                  mid=polyci.mid, polyci.l=polyci.l, polyci.r=polyci.r)
        data.by$data.polyci <- data.polyci
      }
    }


    ################ CI ####################
    if (cimean+cigrid !=0 & !noplot & !ci.fewobs & !fewobs) {
      ciON <- T
    }
    if (ciON) {
      ci.x <- ci.bin <- ci.isknot <- ci.mid <- c()
      if (cimean!=0) {
        if (!is.null(xmean)) {
          ci.x <- c(ci.x, xmean)
        } else {
          xcat   <- findInterval(x.sub, knot, rightmost.closed = T, left.open = T)
          ci.x <- c(ci.x, as.vector(tapply(x.sub, xcat, FUN=mean)))
        }
        ci.bin <- c(ci.bin, 1:nbins)
        ci.isknot <- c(ci.isknot, rep(0, nbins))
        ci.mid <- c(ci.mid, rep(0, nbins))
      }
      if (cigrid!=0) {
        grid <- binsreg.grid(knot, cigrid, addmore=T)
        ci.x <- c(ci.x, grid$eval); ci.bin <- c(ci.bin, grid$bin)
        ci.isknot <- c(ci.isknot, grid$isknot); ci.mid <- c(ci.mid, grid$mid)
      }

      ci.reg.ON <- T
      if (lineON) if (ci.p==line.p & ci.s==line.s) {
        model.ci <- model.line; ci.reg.ON <- F
      }
      if (ci.reg.ON) if (dotsON) if (ci.p==dots.p & ci.s==dots.s) {
        model.ci <- model.dots; ci.reg.ON <- F
      }
      if (ci.reg.ON) {
        B    <- binsreg.spdes(eval=x.sub, p=ci.p, s=ci.s, deriv=0, knot=knot)
        P    <- cbind(B, w.sub)            # full design matrix
        model.ci <- rq(y.sub~P-1, tau=quantile, weights=weights.sub)
        check.drop(model.ci$coeff, ncol(B))
      }

      basis <- binsreg.spdes(eval=ci.x, p=ci.p, s=ci.s, knot=knot, deriv=deriv)
      ci.pred <- binsreg.pred(X=basis, model=model.ci, type="all", vce=vce,
                              cluster=cluster.sub, deriv=deriv, wvec=eval.w, is.qreg=TRUE, avar = asyvar, ...)
      ci.l <- ci.pred$fit - qnorm(alpha)*ci.pred$se
      ci.r <- ci.pred$fit + qnorm(alpha)*ci.pred$se
      ci.l[ci.isknot==1] <- NA
      ci.r[ci.isknot==1] <- NA
      data.ci <- data.frame(group=byvals[i], x=ci.x, bin=ci.bin, isknot=ci.isknot,
                            mid=ci.mid, ci.l=ci.l, ci.r=ci.r)
      data.by$data.ci <- data.ci
    }

    ################ CB ###############################
    cval <- NA
    if (cbgrid !=0 & !noplot & !cb.fewobs & !fewobs) {
      cbON <- T
    }
    if (cbON) {
      grid <- binsreg.grid(knot, cbgrid, addmore=T)
      cb.x <- grid$eval; cb.bin <- grid$bin
      cb.isknot <- grid$isknot; cb.mid <- grid$mid

      cb.reg.ON <- T
      if (ciON) if (cb.p==ci.p & cb.s==ci.s) {
        model.cb <- model.ci; cb.reg.ON <- F
      }
      if (cb.reg.ON) if (lineON) if (cb.p==line.p & cb.s==line.s) {
        model.cb <- model.line; cb.reg.ON <- F
      }
      if (cb.reg.ON) if (dotsON) if (cb.p==dots.p & cb.s==dots.s) {
        model.cb <- model.dots; cb.reg.ON <- F
      }

      if (cb.reg.ON) {
        B    <- binsreg.spdes(eval=x.sub, p=cb.p, s=cb.s, deriv=0, knot=knot)
        P    <- cbind(B, w.sub)            # full design matrix
        model.cb <- rq(y.sub~P-1, tau=quantile, weights=weights.sub)
        check.drop(model.cb$coeff, ncol(B))
      }

      basis <- binsreg.spdes(eval=cb.x, p=cb.p, s=cb.s, knot=knot, deriv=deriv)
      pos <- !is.na(model.cb$coeff[1:ncol(basis)])
      k.new <- sum(pos)
      cb.pred <- binsreg.pred(X=basis, model=model.cb, type="all", vce=vce,
                              cluster=cluster.sub, deriv=deriv, wvec=eval.w, is.qreg=TRUE, avar=asyvar, ...)

      ### Compute cval ####
      x.grid <- binsreg.grid(knot, simsgrid)$eval
      basis.sim <- binsreg.spdes(eval=x.grid, p=cb.p, s=cb.s, knot=knot, deriv=deriv)
      sim.pred <- binsreg.pred(X=basis.sim, model=model.cb, type="all",
                               vce=vce, cluster=cluster.sub, is.qreg=TRUE, avar=T, ...)
      vcv <- binsreg.vcov(model.cb, type=vce, cluster=cluster.sub, is.qreg=TRUE, ...)[1:k.new, 1:k.new]
      Sigma.root <- lssqrtm(vcv)
      num        <- basis.sim[,pos,drop=F] %*% Sigma.root
      cval <- binsreg.pval(num, sim.pred$se, rep=nsims, tstat=NULL, side="two", alpha=level)$cval

      cb.l <- cb.pred$fit - cval*cb.pred$se
      cb.r <- cb.pred$fit + cval*cb.pred$se
      if (cb.s == 0 | cb.s - deriv <=0) {
        cb.l[cb.isknot==1] <- NA
        cb.r[cb.isknot==1] <- NA
      }
      data.cb <- data.frame(group=byvals[i], x=cb.x, bin=cb.bin, isknot=cb.isknot,
                            mid=cb.mid, cb.l=cb.l, cb.r=cb.r)
      data.by$data.cb <- data.cb
    }
    cval.by <- c(cval.by, cval)

    # Save all data for each group
    data.plot[[i]] <- data.by
    names(data.plot)[i] <- paste("Group", byvals[i], sep=" ")

  }

  ########################################
  ############# Ploting ? ################
  binsplot <- NULL
  if (!noplot) {
    binsplot <- ggplot() + theme_bw()
    x <- fit <- ci.l <- ci.r <- cb.l <- cb.r <- polyci.l <- polyci.r <- group <- NULL
    xr.min <- yr.min <- -Inf; xr.max <- yr.max <- Inf
    if (!is.null(plotxrange)) {
      xr.min <- plotxrange[1]
      if (length(plotxrange)==2) xr.max <- plotxrange[2]
    }
    if (!is.null(plotyrange)) {
      yr.min <- plotyrange[1]
      if (length(plotyrange)==2) yr.max <- plotyrange[2]
    }

    for (i in 1:ngroup) {
      data.by <- data.plot[[i]]
      if (!is.null(data.by$data.dots)) {
        index <- complete.cases(data.by$data.dots[c("x", "fit")]) & (data.by$data.dots["x"]>=xr.min) & (data.by$data.dots["x"]<=xr.max) &
                   (data.by$data.dots["fit"]>=yr.min) & (data.by$data.dots["fit"]<=yr.max)
        if (!legendoff) {
          binsplot <- binsplot +
            geom_point(data=data.by$data.dots[index,],
                       aes(x=x, y=fit, colour=group), shape=bysymbols[i], size=2)
        } else {
          binsplot <- binsplot +
            geom_point(data=data.by$data.dots[index,],
                       aes(x=x, y=fit), col=bycolors[i], shape=bysymbols[i], size=2)
        }
      }
      if (!is.null(data.by$data.line)) {
        index <- (data.by$data.line["x"]>=xr.min) & (data.by$data.line["x"]<=xr.max) &
                 (((data.by$data.line["fit"]>=yr.min) & (data.by$data.line["fit"]<=yr.max)) | is.na(data.by$data.line["fit"]))
        if (!legendoff) {
          binsplot <- binsplot +
            geom_line(data=data.by$data.line[index,], aes(x=x, y=fit, colour=as.factor(group)),
                      linetype=bylpatterns[i])
        } else {
        binsplot <- binsplot +
          geom_line(data=data.by$data.line[index,], aes(x=x, y=fit),
                    col=bycolors[i], linetype=bylpatterns[i])
        }
      }
      if (!is.null(data.by$data.poly)) {
        index <- (data.by$data.poly["x"]>=xr.min) & (data.by$data.poly["x"]<=xr.max) & (data.by$data.poly["fit"]>=yr.min) & (data.by$data.poly["fit"]<=yr.max)
        binsplot <- binsplot +
          geom_line(data=data.by$data.poly[index,], aes(x=x, y=fit),
                    col=bycolors[i], linetype=bylpatterns[i])
      }
      if (!is.null(data.by$data.polyci)) {
        index <- (data.by$data.polyci["x"]>=xr.min) & (data.by$data.polyci["x"]<=xr.max) & (data.by$data.polyci["polyci.l"]>=yr.min) & (data.by$data.polyci["polyci.r"]<=yr.max)
        binsplot <- binsplot +
          geom_errorbar(data=data.by$data.polyci[index,],
                        aes(x=x, ymin=polyci.l, ymax=polyci.r),
                        alpha=1, col=bycolors[i], linetype=bylpatterns[i])
      }
      if (!is.null(data.by$data.ci)) {
        index <- complete.cases(data.by$data.ci[c("x", "ci.l", "ci.r")]) & (data.by$data.ci["x"]>=xr.min) & (data.by$data.ci["x"]<=xr.max) &
                   (data.by$data.ci["ci.l"]>=yr.min) & (data.by$data.ci["ci.r"]<=yr.max)
        binsplot <- binsplot +
          geom_errorbar(data=data.by$data.ci[index,],
                        aes(x=x, ymin=ci.l, ymax=ci.r),
                        alpha=1, col=bycolors[i], linetype=bylpatterns[i])
      }
      if (!is.null(data.by$data.cb)) {
        index <- (data.by$data.cb["x"]>=xr.min) & (data.by$data.cb["x"]<=xr.max) &
                 (((data.by$data.cb["cb.l"]>=yr.min) & (data.by$data.cb["cb.r"]<=yr.max)) | is.na(data.by$data.cb["cb.l"]))
        binsplot <- binsplot +
          geom_ribbon(data=data.by$data.cb[index,], aes(x=x, ymin=cb.l, ymax=cb.r),
                      alpha=0.2, fill=bycolors[i])
      }
    }

    # Add legend ?
    if (!legendoff) {
       binsplot <- binsplot +
                   scale_color_manual(name=legendTitle, values = bycolors[1:ngroup],
                              guide=guide_legend(override.aes = list(
                              linetype=bylpatterns[1:ngroup], shape=bysymbols[1:ngroup])))
    } else {
       binsplot <- binsplot + theme(legend.position="none")
    }
    binsplot <- binsplot + labs(x=paste(xname), y=paste(yname))
    print(binsplot)
  }

  ######################################
  ########### Output ###################
  ######################################
  out <- list(bins_plot=binsplot, data.plot=data.plot, cval.by=cval.by,
              opt=list(dots=dots, line=line, ci=ci, cb=cb,
                       polyreg=polyreg, deriv=deriv, quantile=quantile,
                       binspos=position, binsmethod=selectmethod,
                       N.by=N.by, Ndist.by=Ndist.by, Nclust.by=Nclust.by,
                       nbins.by=nbins.by, byvals=byvals))
  out$call <- match.call()
  class(out) <- "CCFFbinsqreg"
  return(out)
}


##########################################################################
#' Internal function.
#'
#' @param x Class \code{CCFFbinsqreg} objects.
#'
#' @keywords internal
#' @export
#'

print.CCFFbinsqreg <- function(x, ...) {
  cat("Call: binsqreg\n\n")

  cat("Binscatter Plot\n")
  cat(paste("Bin selection method (binsmethod)  =  ", x$opt$binsmethod,   "\n", sep=""))
  cat(paste("Placement (binspos)                =  ", x$opt$binspos,      "\n", sep=""))
  cat(paste("Derivative (deriv)                 =  ", x$opt$deriv,        "\n", sep=""))
  cat("\n")
  for (i in 1:length(x$opt$byvals)) {
  cat(paste("Group (by)                         =  ", x$opt$byvals[i],    "\n", sep=""))
  cat(paste("Sample size (n)                    =  ", x$opt$N.by[i],      "\n", sep=""))
  cat(paste("# of distinct values (Ndist)       =  ", x$opt$Ndist.by[i],  "\n", sep=""))
  cat(paste("# of clusters (Nclust)             =  ", x$opt$Nclust.by[i], "\n", sep=""))
  cat(paste("dots, degree (p)                   =  ", x$opt$dots[1],      "\n", sep=""))
  cat(paste("dots, smooth (s)                   =  ", x$opt$dots[2],      "\n", sep=""))
  cat(paste("# of bins (nbins)                  =  ", x$opt$nbins.by[i],  "\n", sep=""))
  cat("\n")
  }
}

################################################################################
#' Internal function.
#'
#' @param object Class \code{CCFFbinsqreg} objects.
#'
#' @keywords internal
#' @export
summary.CCFFbinsqreg <- function(object, ...) {
  x <- object
  args <- list(...)

  cat("Call: binsqreg\n\n")

  cat("Binscatter Plot\n")
  cat(paste("Bin selection method (binsmethod)  =  ", x$opt$binsmethod,   "\n", sep=""))
  cat(paste("Placement (binspos)                =  ", x$opt$binspos,      "\n", sep=""))
  cat(paste("Derivative (deriv)                 =  ", x$opt$deriv,        "\n", sep=""))
  cat("\n")
  for (i in 1:length(x$opt$byvals)) {
  cat(paste("Group (by)                         =  ", x$opt$byvals[i],    "\n", sep=""))
  cat(paste("Sample size (n)                    =  ", x$opt$N.by[i],      "\n", sep=""))
  cat(paste("# of distinct values (Ndist)       =  ", x$opt$Ndist.by[i],  "\n", sep=""))
  cat(paste("# of clusters (Nclust)             =  ", x$opt$Nclust.by[i], "\n", sep=""))
  cat(paste("dots, degree (p)                   =  ", x$opt$dots[1],      "\n", sep=""))
  cat(paste("dots, smooth (s)                   =  ", x$opt$dots[2],      "\n", sep=""))
  cat(paste("# of bins (nbins)                  =  ", x$opt$nbins.by[i],  "\n", sep=""))
  cat("\n")

  cat(paste(rep("=", 8 + 10 + 10 + 10), collapse="")); cat("\n")
  cat(format(" ",  width= 8, justify="right"))
  cat(format("p",  width=10, justify="right"))
  cat(format("s",  width=10, justify="right"))
  cat(format("df", width=10, justify="right"))
  cat("\n")

  cat(paste(rep("-", 8 + 10 + 10 + 10), collapse="")); cat("\n")
  cat(format("dots",  width= 8, justify="right"))
  cat(format(sprintf("%3.0f", x$opt$dots[1]), width=10 , justify="right"))
  cat(format(sprintf("%3.0f", x$opt$dots[2]), width=10 , justify="right"))
  dots.df <- x$opt$dots[1]+1+(x$opt$nbins.by[i]-1)*(x$opt$dots[1]-x$opt$dots[2]+1)
  cat(format(sprintf("%3.0f", dots.df), width=10 , justify="right"))

  if (!is.null(x$opt$line)) {
    cat("\n")
    cat(format("line",  width= 8, justify="right"))
    cat(format(sprintf("%3.0f", x$opt$line[1]), width=10 , justify="right"))
    cat(format(sprintf("%3.0f", x$opt$line[2]), width=10 , justify="right"))
    line.df <- x$opt$line[1]+1+(x$opt$nbins.by[i]-1)*(x$opt$line[1]-x$opt$line[2]+1)
    cat(format(sprintf("%3.0f", line.df), width=10 , justify="right"))
  }

  if (!is.null(x$opt$ci)) {
    cat("\n")
    cat(format("CI",  width= 8, justify="right"))
    cat(format(sprintf("%3.0f", x$opt$ci[1]), width=10 , justify="right"))
    cat(format(sprintf("%3.0f", x$opt$ci[2]), width=10 , justify="right"))
    ci.df <- x$opt$ci[1]+1+(x$opt$nbins.by[i]-1)*(x$opt$ci[1]-x$opt$ci[2]+1)
    cat(format(sprintf("%3.0f", ci.df), width=10 , justify="right"))
  }

  if (!is.null(x$opt$cb)) {
    cat("\n")
    cat(format("CB",  width= 8, justify="right"))
    cat(format(sprintf("%3.0f", x$opt$cb[1]), width=10 , justify="right"))
    cat(format(sprintf("%3.0f", x$opt$cb[2]), width=10 , justify="right"))
    cb.df <- x$opt$cb[1]+1+(x$opt$nbins.by[i]-1)*(x$opt$cb[1]-x$opt$cb[2]+1)
    cat(format(sprintf("%3.0f", cb.df), width=10 , justify="right"))
  }

  if (!is.null(x$opt$polyreg)) {
    cat("\n")
    cat(format("polyreg",  width= 8, justify="right"))
    cat(format(sprintf("%3.0f", x$opt$polyreg), width=10 , justify="right"))
    cat(format(NA, width=10 , justify="right"))
    cat(format(sprintf("%3.0f", x$opt$polyreg+1), width=10 , justify="right"))
  }

  cat("\n")
  cat(paste(rep("-", 8 + 10 + 10 + 10), collapse=""))
  cat("\n")
  cat("\n")
  }

  #cat("\n")
}

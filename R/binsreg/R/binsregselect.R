######################################################################################
#'@title Data-Driven IMSE-Optimal Partitioning/Binning Selection for Binscatter
#'@description \code{binsregselect} implements data-driven procedures for selecting the number of bins for binscatter
#'             estimation. The selected number is optimal in minimizing integrated mean squared error (IMSE).
#'@param  y outcome variable. A vector.
#'@param  x independent variable of interest. A vector.
#'@param  w control variables. A matrix, a vector or a \code{\link{formula}}.
#'@param  data an optional data frame containing variables used in the model.
#'@param  deriv  derivative order of the regression function for estimation, testing and plotting.
#'               The default is \code{deriv=0}, which corresponds to the function itself.
#'@param  bins a vector. \code{bins=c(p,s)} set a piecewise polynomial of degree \code{p} with \code{s} smoothness constraints
#'             for data-driven (IMSE-optimal) selection of the partitioning/binning scheme. By default, the function sets
#'             \code{bins=c(0,0)}, which corresponds to piecewise constant (canonical binscatter).
#'@param  pselect vector of numbers within which the degree of polynomial \code{p} for point estimation is selected.
#'                \emph{Note:} To implement the degree or smoothness selection, in addition to \code{pselect} or \code{sselect},
#'                \code{nbins=#} must be specified.
#'@param  sselect vector of numbers within which the number of smoothness constraints \code{s} for point estimation is selected.
#'                If not specified, for each value \code{p} supplied in the option \code{pselect}, only the
#'                piecewise polynomial with the maximum smoothness is considered, i.e., \code{s=p}.
#'@param  binspos position of binning knots.  The default is \code{binspos="qs"}, which corresponds to quantile-spaced
#'                binning (canonical binscatter).  The other option is \code{binspos="es"} for evenly-spaced binning.
#'@param  nbins number of bins for degree/smoothness selection. If \code{nbins=T} or \code{nbins=NULL} (default) is specified,
#'              the function selects the number of bins instead, given the specified degree and smoothness.
#'              If a vector with more than one number is specified, the command selects the number of bins within this vector.
#'@param  binsmethod method for data-driven selection of the number of bins. The default is \code{binsmethod="dpi"},
#'                   which corresponds to the IMSE-optimal direct plug-in rule.  The other option is: \code{"rot"}
#'                   for rule of thumb implementation.
#'@param  nbinsrot initial number of bins value used to construct the DPI number of bins selector.
#'                 If not specified, the data-driven ROT selector is used instead.
#'@param  simsgrid number of evaluation points of an evenly-spaced grid within each bin used for evaluation of
#'                 the supremum (infimum or Lp metric) operation needed to construct confidence bands and hypothesis testing
#'                 procedures. The default is \code{simsgrid=20}, which corresponds to 20 evenly-spaced
#'                 evaluation points within each bin for approximating the supremum (infimum or Lp metric) operator.
#'@param  savegrid if true, a data frame produced containing grid.
#'@param  vce procedure to compute the variance-covariance matrix estimator. Options are
#'           \itemize{
#'           \item \code{"const"} homoskedastic variance estimator.
#'           \item \code{"HC0"} heteroskedasticity-robust plug-in residuals variance estimator
#'                              without weights.
#'           \item \code{"HC1"} heteroskedasticity-robust plug-in residuals variance estimator
#'                              with hc1 weights. Default.
#'           \item \code{"HC2"} heteroskedasticity-robust plug-in residuals variance estimator
#'                              with hc2 weights.
#'           \item \code{"HC3"} heteroskedasticity-robust plug-in residuals variance estimator
#'                              with hc3 weights.
#'           }
#'@param  useeffn effective sample size to be used when computing the (IMSE-optimal) number of bins. This option
#'                is useful for extrapolating the optimal number of bins to larger (or smaller) datasets than
#'                the one used to compute it.
#'@param  randcut upper bound on a uniformly distributed variable used to draw a subsample for bins/degree/smoothness selection.
#'                Observations for which \code{runif()<=#} are used. # must be between 0 and 1.
#'@param  cluster cluster ID. Used for compute cluster-robust standard errors.
#'@param  dfcheck adjustments for minimum effective sample size checks, which take into account number of unique
#'                values of \code{x} (i.e., number of mass points), number of clusters, and degrees of freedom of
#'                the different statistical models considered. The default is \code{dfcheck=c(20, 30)}.
#'                See \href{https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2024_Stata.pdf}{Cattaneo, Crump, Farrell and Feng (2024c)} for more details.
#'@param  masspoints how mass points in \code{x} are handled. Available options:
#'                   \itemize{
#'                   \item \code{"on"} all mass point and degrees of freedom checks are implemented. Default.
#'                   \item \code{"noadjust"} mass point checks and the corresponding effective sample size adjustments are omitted.
#'                   \item \code{"nolocalcheck"} within-bin mass point and degrees of freedom checks are omitted.
#'                   \item \code{"off"} "noadjust" and "nolocalcheck" are set simultaneously.
#'                   \item \code{"veryfew"} forces the function to proceed as if \code{x} has only a few number of mass points (i.e., distinct values).
#'                                          In other words, forces the function to proceed as if the mass point and degrees of freedom checks were failed.
#'                   }
#'@param  norotnorm if true, a uniform density rather than normal density used for ROT selection.
#'@param  numdist  number of distinct values for selection. Used to speed up computation.
#'@param  numclust number of clusters for selection. Used to speed up computation.
#'@param  weights an optional vector of weights to be used in the fitting process. Should be \code{NULL} or
#'                a numeric vector. For more details, see \code{\link{lm}}.
#'@param  subset optional rule specifying a subset of observations to be used.
#'@return \item{\code{nbinsrot.poly}}{ROT number of bins, unregularized.}
#'        \item{\code{nbinsrot.regul}}{ROT number of bins, regularized.}
#'        \item{\code{nbinsrot.uknot}}{ROT number of bins, unique knots.}
#'        \item{\code{nbinsdpi}}{DPI number of bins.}
#'        \item{\code{nbinsdpi.uknot}}{DPI number of bins, unique knots.}
#'        \item{\code{prot.poly}}{ROT degree of polynomials, unregularized.}
#'        \item{\code{prot.regul}}{ROT degree of polynomials, regularized.}
#'        \item{\code{prot.uknot}}{ROT degree of polynomials, unique knots.}
#'        \item{\code{pdpi}}{DPI degree of polynomials.}
#'        \item{\code{pdpi.uknot}}{DPI degree of polynomials, unique knots.}
#'        \item{\code{srot.poly}}{ROT number of smoothness constraints, unregularized.}
#'        \item{\code{srot.regul}}{ROT number of smoothness constraints, regularized.}
#'        \item{\code{srot.uknot}}{ROT number of smoothness constraints, unique knots.}
#'        \item{\code{sdpi}}{DPI number of smoothness constraints.}
#'        \item{\code{sdpi.uknot}}{DPI number of smoothness constraints, unique knots.}
#'        \item{\code{imse.var.rot}}{Variance constant in IMSE expansion, ROT selection.}
#'        \item{\code{imse.bsq.rot}}{Bias constant in IMSE expansion, ROT selection.}
#'        \item{\code{imse.var.dpi}}{Variance constant in IMSE expansion, DPI selection.}
#'        \item{\code{imse.bsq.dpi}}{Bias constant in IMSE expansion, DPI selection.}
#'        \item{\code{int.result}}{Intermediate results, including a matrix of degree and smoothness (\code{deg_mat}),
#'                                 the selected numbers of bins (\code{vec.nbinsrot.poly},\code{vec.nbinsrot.regul},
#'                                 \code{vec.nbinsrot.uknot}, \code{vec.nbinsdpi}, \code{vec.nbinsdpi.uknot}),
#'                                 and the bias and variance constants in IMSE (\code{vec.imse.b.rot},
#'                                 \code{vec.imse.v.rot}, \code{vec.imse.b.dpi}, \code{vec.imse.v.dpi})
#'                                 under each rule (ROT or DPI), corresponding to each pair of degree and smoothness
#'                                 (each row in \code{deg_mat}).}
#'        \item{\code{opt}}{ A list containing options passed to the function, as well as total sample size \code{n},
#'                           number of distinct values \code{Ndist} in \code{x}, and number of clusters \code{Nclust}.}
#'        \item{\code{data.grid}}{A data frame containing grid.}
#'@author
#' Matias D. Cattaneo, Princeton University, Princeton, NJ. \email{cattaneo@princeton.edu}.
#'
#' Richard K. Crump, Federal Reserve Bank of New York, New York, NY. \email{richard.crump@ny.frb.org}.
#'
#' Max H. Farrell, UC Santa Barbara, Santa Barbara, CA. \email{mhfarrell@gmail.com}.
#'
#' Yingjie Feng (maintainer), Tsinghua University, Beijing, China. \email{fengyingjiepku@gmail.com}.
#'
#'@references
#' Cattaneo, M. D., R. K. Crump, M. H. Farrell, and Y. Feng. 2024a: \href{https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2024_AER.pdf}{On Binscatter}. American Economic Review 114(5): 1488-1514.
#'
#' Cattaneo, M. D., R. K. Crump, M. H. Farrell, and Y. Feng. 2024b: \href{https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2024_NonlinearBinscatter.pdf}{Nonlinear Binscatter Methods}. Working Paper.
#'
#' Cattaneo, M. D., R. K. Crump, M. H. Farrell, and Y. Feng. 2024c: \href{https://nppackages.github.io/references/Cattaneo-Crump-Farrell-Feng_2024_Stata.pdf}{Binscatter Regressions}. Working Paper.
#'
#'@seealso \code{\link{binsreg}}, \code{\link{binstest}}.
#'
#'@examples
#'  x <- runif(500); y <- sin(x)+rnorm(500)
#'  est <- binsregselect(y,x)
#'  summary(est)
#'@export

binsregselect <- function(y, x, w=NULL, data=NULL, deriv=0,
                          bins=NULL, pselect=NULL, sselect=NULL,
                          binspos="qs", nbins=NULL, binsmethod="dpi", nbinsrot=NULL,
                          simsgrid=20, savegrid=F,
                          vce="HC1", useeffn=NULL, randcut=NULL, cluster=NULL,
                          dfcheck=c(20,30), masspoints="on", weights=NULL, subset=NULL,
                          norotnorm=F, numdist=NULL, numclust=NULL) {

  # param for internal use
  rot.lb <- 1
  qrot <- 2

  ####################
  ### prepare data ###
  ####################
  xname <- deparse(substitute(x), width.cutoff = 500L)
  yname <- deparse(substitute(y), width.cutoff = 500L)
  weightsname <- deparse(substitute(weights), width.cutoff = 500L)
  subsetname  <- deparse(substitute(subset), width.cutoff = 500L)
  clustername <- deparse(substitute(cluster), width.cutoff = 500L)

  # capture w names for grid file
  wname <- NULL

  # extract y, x, w, weights, subset, if needed (w as a matrix, others as vectors)
  # generate design matrix for covariates W
  if (is.null(data)) {
    if (is.data.frame(y)) y <- y[,1]
    if (is.data.frame(x)) x <- x[,1]
    if (!is.null(w)) {
      if (is.matrix(w)) wname <- colnames(w)
      if (is.vector(w)) {
        wname <- deparse(substitute(w), width.cutoff = 500L)
        w <- as.matrix(w)
        if (is.character(w)) {
          w <- model.matrix(~w)[,-1,drop=F]
        }
      } else if (is.formula(w)) {
        wname <- all.vars(w)
        w.model <- binsreg.model.mat(w)
        w <- w.model$design
      } else if (is.data.frame(w)) {
        w <- as.matrix(w)
      }
    }
  } else {
    if (xname %in% names(data)) x <- data[,xname]
    if (yname %in% names(data)) y <- data[,yname]
    if (weightsname != "NULL") if (weightsname %in% names(data)) {
      weights <- data[,weightsname]
    }
    if (subsetname != "NULL") if (subsetname %in% names(data))  {
      subset  <- data[,subsetname]
    }
    if (clustername != "NULL") if (clustername %in% names(data)) {
      cluster <- data[,clustername]
    }
    if (deparse(substitute(w), width.cutoff = 500L)[1]!="NULL") {
      if (is.formula(w)) {
        wname <- all.vars(w)
        w.model <- binsreg.model.mat(w, data=data)
        w <- w.model$design
      } else {
        wname <- deparse(substitute(w), width.cutoff = 500L)
        if (wname %in% names(data)) w <- data[,wname]
        w <- as.matrix(w)
        if (is.character(w)) w <- model.matrix(~w)[,-1,drop=F]
      }
    }
  }

  if (!is.null(w)) nwvar <- ncol(w)
  else             nwvar <- 0

  # extract subset
  if (!is.null(subset)) {
    y <- y[subset]
    x <- x[subset]
    if (!is.null(w)) w <- w[subset, , drop = F]
    if (!is.null(cluster)) cluster <- cluster[subset]
    if (!is.null(weights)) weights <- weights[subset]
  }

  na.ok <- complete.cases(x) & complete.cases(y)
  if (!is.null(w)) na.ok <- na.ok & complete.cases(w)
  if (!is.null(weights)) na.ok <- na.ok & complete.cases(weights)
  if (!is.null(cluster)) na.ok <- na.ok & complete.cases(cluster)

  y <- y[na.ok]
  x <- x[na.ok]
  w <- w[na.ok, , drop = F]
  weights <- weights[na.ok]
  cluster <- cluster[na.ok]


  ## analyze bins- and degree-related options
  if (!is.null(bins)) {
    bins.p <- bins[1]
    if (length(bins)==1) {
      bins.s <- bins.p
    } else if (length(bins)==2) {
      bins.s <- bins[2]
    } else {
      print("Bins not correctly specified.")
      stop()
    }
    plist <- bins[1]; slist <- bins[2]
  } else {
    plist <- pselect; slist <- sselect
  }

  len_p <- length(plist); len_s <- length(slist)
  if (len_p==1 & len_s==0) {
    slist <- plist; len_s <- 1
  }
  if (len_p==0 & len_s==1) {
    plist <- slist; len_p <- 1
  }


  selectJ <- NULL
  if (is.logical(nbins)) if (nbins) selectJ <- T
  if (length(nbins)==1) if (nbins==0) selectJ <- T
  if (is.null(nbins)) selectJ <- T
  if (length(nbins)>1|!is.null(bins)|(len_p==1&len_s==1)) {
    selectJ <- T                      # turn on J selection
  }

  if (is.logical(selectJ)) if (selectJ) {
    if (len_p>1|len_s>1) {
      print("Only one p and one s are allowed for J selection.")
      stop()
    }
    if (is.null(plist)) {
      plist <- deriv; len_p <- 1
    }
    if (is.null(slist)) {
      slist <- plist; len_s <- 1
    }
  }

  if (is.null(selectJ) & (len_p>1|len_s>1)) {
    selectJ <- F
  }

  if (is.null(selectJ)) {
    print("Degree, smoothness, or # of bins are not correctly specified.")
    stop()
  }

  # find all compatible pairs
  if (selectJ) {
    deg_mat <- t(c(plist, slist))
  } else {
    if (len_p>0 & len_s==0) deg_mat <- cbind(plist, plist)
    if (len_p==0 & len_s>0) deg_mat <- cbind(slist, slist)
    if (len_p>0 & len_s>0)  {
      deg_mat <- as.matrix(expand.grid(plist, slist))
      comp.ind <- (deg_mat[,1]>=deg_mat[,2] & deg_mat[,1]>deriv)     # p>s and p>deriv
      if (sum(comp.ind)==0) {
        print("Degree and smoothness incompatible")
        stop()
      }
      deg_mat <- deg_mat[comp.ind,,drop=F]
    }
  }
  colnames(deg_mat) <- c("degree", "smoothness")


  #############################error checking
  exit <- 0
  if (length(bins)==2) if (bins[1]<bins[2]) {
    print("p<s not allowed.")
    exit <- 1
  }
  if (binsmethod!="dpi" & binsmethod!="rot") {
    print("bin selection method incorrectly specified.")
    exit <- 1
  }
  if (binspos!="es" & binspos!="qs") {
    print("binspos incorrectly specified.")
    exit <- 1
  }
  if (exit>0) stop()

  #####################################
  rot.fewobs <- dpi.fewobs <- F
  localcheck <- massadj <- T
  if (masspoints=="on") {
    localcheck <- T; massadj <- T
  } else if (masspoints=="off") {
    localcheck <- F; massadj <- F
  } else if (masspoints=="noadjust") {
    localcheck <- T; massadj <- F
  } else if (masspoints=="nolocalcheck") {
    localcheck <- F; massadj <- T
  } else if (masspoints=="veryfew") {
    rot.fewobs <- dpi.fewobs <- T
  }

  # effective size
  eN <- N <- length(x)
  Ndist <- NA
  if (massadj) {
    if (!is.null(numdist)) {
      Ndist <- numdist
    } else {
      Ndist <- length(unique(x))
    }
    eN <- min(eN, Ndist)
  }
  Nclust <- NA
  if (!is.null(cluster)) {
    if (!is.null(numclust)) {
      Nclust <- numclust
    } else {
      Nclust <- length(unique(cluster))
    }
    eN <- min(eN, Nclust)
  }

  # take a subsample?
  x.full <- x
  eN.sub <- eN; Ndist.sub <- Ndist; Nclust.sub <- Nclust
  if (!is.null(randcut)) {
    subsample <- (runif(N)<=randcut)
    x <- x[subsample]
    y <- y[subsample]
    w <- w[subsample,, drop=F]
    weights <- weights[subsample]
    cluster <- cluster[subsample]

    eN.sub <- length(x)
    if (massadj) {
      Ndist.sub <- length(unique(x))
      eN.sub    <- min(eN.sub, Ndist.sub)
    }
    if (!is.null(cluster)) {
      Nclust.sub <- length(unique(cluster))
      eN.sub     <- min(eN.sub, Nclust.sub)
    }
  }

  # redefine cluster if needed
  if (massadj) {
    if (is.null(cluster)) {
      cluster <- x
      # note: clustered at mass point level
    } else {
      if (Nclust.sub > Ndist.sub) {
        cluster <- x
        warning("# mass points < # clusters. Clustered at mass point level.")
      }
    }
  }

  # Prepare params
  #p <- bins[1]; s <- bins[2]

  if (binspos == "es") {
    es <- T
  } else {
    es <- F
  }

  # Store options
  if (es)  {
    position <- "Evenly-spaced"
  } else {
    position <- "Quantile-spaced"
  }

  if (binsmethod == "dpi") {
    selectmethod <- "IMSE direct plug-in"
  } else {
    selectmethod <- "IMSE rule-of-thumb"
  }
  if (selectJ) selectmethod <- paste(selectmethod, "(select # of bins)")
  else         selectmethod <- paste(selectmethod, "(select degree and smoothness)")

  vec.J.rot.poly <- vec.J.rot.regul <- vec.J.rot.uniq  <- vec.J.dpi <- vec.J.dpi.uniq <- c()
  vec.imse.v.rot <- vec.imse.b.rot <- vec.imse.v.dpi <- vec.imse.b.dpi <- c()
  #### START loop here #########
  for (j in 1:nrow(deg_mat)) {
    p <- deg_mat[j,1]; s <- deg_mat[j,2]

    # Run rot selection
    J.rot.regul <- J.rot.poly <- NA; imse.v.rot <- imse.b.rot <- NA
    if (!is.null(nbinsrot)) J.rot.regul <- nbinsrot
    if (is.na(J.rot.regul) & !rot.fewobs) {
      # checking
      if (eN.sub <= dfcheck[1]+p+1+qrot) {
        rot.fewobs <- T
        warning("Too small effective sample size for bin selection.")
      }

      if (!rot.fewobs) {
        J.rot.poly <- binsregselect.rot(y, x, w, p, s, deriv, es=es, eN=eN.sub, qrot=qrot, norotnorm=norotnorm, weights=weights)
        imse.v.rot <- J.rot.poly$imse.v
        imse.b.rot <- J.rot.poly$imse.b
        J.rot.poly <- J.rot.poly$J.rot
      }
      J.rot.regul <- max(J.rot.poly, ceiling((2*(p+1-deriv)/(1+2*deriv)*rot.lb*eN.sub)^(1/(2*p+3))))
    }
    # repeated knots?
    J.rot.uniq <- J.rot.regul
    if (!es & !is.na(J.rot.regul)) {
       J.rot.uniq <- length(unique(genKnot.qs(x, J.rot.regul)[-1]))
    }

    # Run dpi selection
    J.dpi <- NA; imse.v.dpi <- imse.b.dpi <- NA
    if (binsmethod == "dpi" & !dpi.fewobs) {
      # check if dpi can be implemented
      if (!is.na(J.rot.uniq)) {
        if ((p-s+1)*(J.rot.uniq-1)+p+2+dfcheck[2]>=eN.sub) {
           dpi.fewobs <- T
           warning("Too small effective sample size for DPI selection.")
        }

        # check empty bins
        if (localcheck) {
           uniqmin <- binsreg.checklocalmass(x, J.rot.regul, es, knot=NULL) # mimic STATA
           if (uniqmin < p+2) {
             dpi.fewobs <- T
             warning("Some bins have too few distinct values of x for DPI selection.")
           }
        }
      } else {
        dpi.fewobs <- T
      }

      if (!dpi.fewobs) {
        J.dpi <- binsregselect.dpi(y, x, w, p, s, deriv, es=es, vce=vce, cluster=cluster, nbinsrot=J.rot.uniq, weights=weights)
        imse.b.dpi <- J.dpi$imse.b
        imse.v.dpi <- J.dpi$imse.v * eN.sub
        J.dpi      <- J.dpi$J.dpi
      }
    }
    J.dpi.uniq <- J.dpi

    if (!is.null(useeffn)|!is.null(randcut)) {
      if (!is.null(useeffn)) scaling <- (useeffn/eN)^(1/(2*p+2+1))
      if (!is.null(randcut)) scaling <- (eN/eN.sub)^(1/(2*p+2+1))

      if (!is.na(J.rot.poly))  J.rot.poly  <- ceiling(J.rot.poly * scaling)
      if (!is.na(J.rot.regul)) J.rot.regul <- ceiling(J.rot.regul * scaling)
      if (!is.na(J.rot.uniq))  J.rot.uniq  <- ceiling(J.rot.uniq * scaling)
      if (!is.na(J.dpi))       J.dpi       <- ceiling(J.dpi * scaling)
      if (!is.na(J.dpi.uniq))  J.dpi.uniq  <- ceiling(J.dpi.uniq * scaling)
    }

    # save results
    vec.J.rot.poly[j]  <- J.rot.poly
    vec.J.rot.regul[j] <- J.rot.regul
    vec.J.rot.uniq[j]  <- J.rot.uniq
    vec.J.dpi[j]       <- J.dpi
    vec.J.dpi.uniq[j] <- J.dpi.uniq
    vec.imse.b.rot[j] <- imse.b.rot
    vec.imse.v.rot[j] <- imse.v.rot
    vec.imse.b.dpi[j] <- imse.b.dpi
    vec.imse.v.dpi[j] <- imse.v.dpi
  }
  ##### END loop here  ########

  # extract results
  imse.b.dpi.upd <- imse.v.dpi.upd <- NA
  if (selectJ) {
    if (length(nbins)>1) {
      J.rot.poly  <- nbins[which.min(abs(vec.J.rot.poly-nbins))]
      J.rot.regul <- nbins[which.min(abs(vec.J.rot.regul-nbins))]
      J.rot.uniq  <- nbins[which.min(abs(vec.J.rot.uniq-nbins))]
      J.dpi       <- nbins[which.min(abs(vec.J.dpi-nbins))]
      J.dpi.uniq  <- nbins[which.min(abs(vec.J.dpi.uniq-nbins))]
    }
    ord.rot.poly <- ord.rot.regul <- ord.rot.uniq <- ord.dpi <- ord.dpi.uniq <- deg_mat
  } else {
    ind.rot.poly <- which.min(abs(vec.J.rot.poly-nbins))
    ind.dpi      <- which.min(abs(vec.J.dpi-nbins))

    ord.rot.poly  <- deg_mat[ind.rot.poly,]
    ord.rot.regul <- deg_mat[which.min(abs(vec.J.rot.regul-nbins)),]
    ord.rot.uniq  <- deg_mat[which.min(abs(vec.J.rot.uniq-nbins)),]
    ord.dpi       <- deg_mat[ind.dpi,]
    ord.dpi.uniq  <- deg_mat[which.min(abs(vec.J.dpi.uniq-nbins)),]
    J.rot.poly <- J.rot.regul <- J.rot.uniq <- J.dpi <- J.dpi.uniq <- nbins

    imse.v.rot <- vec.imse.v.rot[ind.rot.poly]
    imse.b.rot <- vec.imse.b.rot[ind.rot.poly]
    imse.v.dpi <- vec.imse.v.dpi[ind.dpi]
    imse.b.dpi <- vec.imse.b.dpi[ind.dpi]

    if (nbins!=vec.J.dpi[ind.dpi]) {
      fit <- binsregselect.dpi(y, x, w, p=ord.dpi[1], s=ord.dpi[2], deriv, es=es, vce=vce, cluster=cluster, nbinsrot=nbins, weights=weights)
      imse.b.dpi.upd <- fit$imse.b
      imse.v.dpi.upd <- fit$imse.v * eN.sub
      imse.b.dpi <- imse.b.dpi.upd
      imse.v.dpi <- imse.v.dpi.upd
    }
  }

  # Generate a knot vector
  if (binsmethod == "rot") {
    Jselect <- J.rot.uniq
  } else {
    Jselect <- J.dpi
  }
  knot <- NA; data.grid <- NA
  if (!is.na(Jselect) & is.null(useeffn)) {
    if (es) {
      knot <- genKnot.es(min(x.full), max(x.full), Jselect)
    } else {
      knot <- genKnot.qs(x.full, Jselect)
    }
    knot <- c(knot[1], unique(knot[-1]))
    Jselect <- length(knot)-1
    if (binsmethod=="dpi") {
      J.dpi.uniq <- Jselect
    }

    # a grid dataset
    if (savegrid) {
      grid <- binsreg.grid(knot=knot, ngrid=simsgrid, addmore=T)
      data.grid <- cbind(grid$eval, grid$bin, grid$isknot)
      data.grid <- cbind(data.grid, matrix(0, nrow(data.grid), length(wname)))
      data.grid <- data.frame(data.grid)
      colnames(data.grid) <- c(xname, "binreg_bin", "binsreg_isknot", wname)
    }
  }

  ######################
  #######output#########
  ######################
  out <- list(nbinsrot.poly=J.rot.poly, nbinsrot.regul=J.rot.regul, nbinsrot.uknot=J.rot.uniq,
              nbinsdpi=J.dpi, nbinsdpi.uknot=J.dpi.uniq,
              imse.bsq.rot=imse.b.rot, imse.var.rot=imse.v.rot,
              imse.bsq.dpi=imse.b.dpi, imse.var.dpi=imse.v.dpi,
              int.result=list(vec.nbinsrot.poly=vec.J.rot.poly, vec.nbinsrot.regul=vec.J.rot.regul,
                              vec.nbinsrot.uknot=vec.J.rot.uniq, vec.nbinsdpi=vec.J.dpi, vec.nbinsdpi.uknot=vec.J.dpi.uniq,
                              vec.imse.b.rot=vec.imse.b.rot, vec.imse.v.rot=vec.imse.v.rot,
                              vec.imse.b.dpi=vec.imse.b.dpi, vec.imse.v.dpi=vec.imse.v.dpi,
                              deg_mat=deg_mat),
              prot.poly=ord.rot.poly[1], srot.poly=ord.rot.poly[2],
              prot.regul=ord.rot.regul[1], srot.regul=ord.rot.regul[2],
              prot.uknot=ord.rot.uniq[1], srot.uknot=ord.rot.uniq[2],
              pdpi=ord.dpi[1], sdpi=ord.dpi[2],
              pdpi.uknot=ord.dpi.uniq[1], sdpi.uknot=ord.dpi.uniq[2],
              opt = list(deriv=deriv, selectJ=selectJ,
                         binspos=position, binsmethod=selectmethod,
                         n=N, Ndist=Ndist, Nclust=Nclust),
              knot=knot, data.grid=data.grid)
  out$call   <- match.call()
  class(out) <- "CCFFbinsregselect"
  return(out)
}

##########################################################################
#' Internal function.
#'
#' @param x Class \code{CCFFbinsregselect} objects.
#'
#' @keywords internal
#' @export
#'
print.CCFFbinsregselect <- function(x, ...) {
  cat("Call: binsregselect\n\n")

  cat(paste("Sample size (n)                    =  ", x$opt$n,          "\n", sep=""))
  cat(paste("# of distinct values (Ndist)       =  ", x$opt$Ndist,      "\n", sep=""))
  cat(paste("# of clusters (Nclust)             =  ", x$opt$Nclust,     "\n", sep=""))
  cat(paste("Derivative (deriv)                 =  ", x$opt$deriv,      "\n", sep=""))
  cat(paste("Bin/Degree selection:", "\n"))
  cat(paste("  Method (binsmethod)              =  ", x$opt$binsmethod, "\n", sep=""))
  cat(paste("  Placement (binspos)              =  ", x$opt$binspos,    "\n", sep=""))
  if (x$opt$selectJ) {
    cat(paste("  degree (p)                       =  ", x$prot.poly,     "\n", sep=""))
    cat(paste("  smooth (s)                       =  ", x$srot.poly,     "\n", sep=""))
    cat(paste("  # of bins (ROT-POLY)             =  ", sprintf("%1.0f",x$nbinsrot.poly),  "\n", sep=""))
    cat(paste("  # of bins (ROT-REGUL)            =  ", sprintf("%1.0f",x$nbinsrot.regul), "\n", sep=""))
    cat(paste("  # of bins (ROT-UKNOT)            =  ", sprintf("%1.0f",x$nbinsrot.uknot), "\n", sep=""))
    if (x$opt$binsmethod=="IMSE direct plug-in (select # of bins)") {
    cat(paste("  # of bins (DPI)                  =  ", sprintf("%1.0f",x$nbinsdpi),       "\n", sep=""))
    cat(paste("  # of bins (DPI-UKNOT)            =  ", sprintf("%1.0f",x$nbinsdpi.uknot), "\n", sep=""))
    }
  } else {
    cat(paste("  degree (ROT-POLY)                =  ", sprintf("%1.0f",x$prot.poly),  "\n", sep=""))
    cat(paste("  degree (ROT-REGUL)               =  ", sprintf("%1.0f",x$prot.regul), "\n", sep=""))
    cat(paste("  degree (ROT-UKNOT)               =  ", sprintf("%1.0f",x$prot.uknot), "\n", sep=""))
    if (x$opt$binsmethod=="IMSE direct plug-in (select degree and smoothness)") {
      cat(paste("  degree (DPI)                     =  ", sprintf("%1.0f",x$pdpi),       "\n", sep=""))
      cat(paste("  degree (DPI-UKNOT)               =  ", sprintf("%1.0f",x$pdpi.uknot), "\n", sep=""))
    }
    cat(paste("  smoothness (ROT-POLY)            =  ", sprintf("%1.0f",x$srot.poly),  "\n", sep=""))
    cat(paste("  smoothness (ROT-REGUL)           =  ", sprintf("%1.0f",x$srot.regul), "\n", sep=""))
    cat(paste("  smoothness (ROT-UKNOT)           =  ", sprintf("%1.0f",x$srot.uknot), "\n", sep=""))
    if (x$opt$binsmethod=="IMSE direct plug-in (select degree and smoothness)") {
      cat(paste("  smoothness (DPI)                 =  ", sprintf("%1.0f",x$sdpi),       "\n", sep=""))
      cat(paste("  smoothness (DPI-UKNOT)           =  ", sprintf("%1.0f",x$sdpi.uknot), "\n", sep=""))
    }
  }
  cat("\n")
}

################################################################################
#' Internal function.
#'
#' @param object Class \code{CCFFbinsregselect} objects.
#'
#' @keywords internal
#' @export
summary.CCFFbinsregselect <- function(object, ...) {
  x <- object
  args <- list(...)

  cat("Call: binsregselect\n\n")

  cat(paste("Sample size (n)                    =  ", x$opt$n,          "\n", sep=""))
  cat(paste("# of distinct values (Ndist)       =  ", x$opt$Ndist,      "\n", sep=""))
  cat(paste("# of clusters (Nclust)             =  ", x$opt$Nclust,     "\n", sep=""))
  cat(paste("Derivative (deriv)                 =  ", x$opt$deriv,      "\n", sep=""))
  cat(paste("Bin/Degree selection:", "\n"))
  cat(paste("  Method (binsmethod)              =  ", x$opt$binsmethod, "\n", sep=""))
  cat(paste("  Placement (binspos)              =  ", x$opt$binspos,    "\n", sep=""))
  cat(paste("  degree (p)                       =  ", x$prot.poly,      "\n", sep=""))
  cat(paste("  smooth (s)                       =  ", x$srot.poly,      "\n", sep=""))
  cat("\n")

  if (x$opt$selectJ) {
    cat(paste(rep("=", 13 + 13 + 13 + 15 + 15), collapse="")); cat("\n")
    cat(format("method",     width=13, justify="right"))
    cat(format("# of bins",  width=13, justify="right"))
    cat(format("df",         width=13, justify="right"))
    cat(format("imse, bias^2",  width=15, justify="right"))
    cat(format("imse, var.",  width=15, justify="right"))
    cat("\n")
    cat(paste(rep("-", 13 + 13 + 13 + 15 + 15), collapse="")); cat("\n")
    cat(format("ROT-POLY",  width= 13, justify="right"))
    cat(format(sprintf("%3.0f", x$nbinsrot.poly), width=13 , justify="right"))
    ROT.poly.df <- NA
    if (!is.na(x$nbinsrot.poly)) ROT.poly.df <- x$prot.poly+1+(x$nbinsrot.poly-1)*(x$prot.poly-x$srot.poly+1)
    cat(format(sprintf("%3.0f", ROT.poly.df), width=13, justify="right"))
    cat(format(sprintf("%3.3f", x$imse.bsq.rot), width=15 , justify="right"))
    cat(format(sprintf("%3.3f", x$imse.var.rot), width=15 , justify="right"))
    cat("\n")

    cat(format("ROT-REGUL",  width= 13, justify="right"))
    cat(format(sprintf("%3.0f", x$nbinsrot.regul), width=13 , justify="right"))
    ROT.regul.df <- NA
    if (!is.na(x$nbinsrot.regul)) ROT.regul.df <- x$prot.regul+1+(x$nbinsrot.regul-1)*(x$prot.regul-x$srot.regul+1)
    cat(format(sprintf("%3.0f", ROT.regul.df), width=13, justify="right"))
    cat(format(sprintf("%3.0f", NA), width=15 , justify="right"))
    cat(format(sprintf("%3.0f", NA), width=15 , justify="right"))
    cat("\n")

    cat(format("ROT-UKNOT",  width= 13, justify="right"))
    cat(format(sprintf("%3.0f", x$nbinsrot.uknot), width=13 , justify="right"))
    ROT.uknot.df <- NA
    if (!is.na(x$nbinsrot.uknot)) ROT.uknot.df <- x$prot.uknot+1+(x$nbinsrot.uknot-1)*(x$prot.uknot-x$srot.uknot+1)
    cat(format(sprintf("%3.0f", ROT.uknot.df), width=13, justify="right"))
    cat(format(sprintf("%3.0f", NA), width=15 , justify="right"))
    cat(format(sprintf("%3.0f", NA), width=15 , justify="right"))
    cat("\n")

    if (x$opt$binsmethod=="IMSE direct plug-in (select # of bins)") {
      cat(format("DPI",  width=13, justify="right"))
      cat(format(sprintf("%3.0f", x$nbinsdpi), width=13 , justify="right"))
      DPI.df <- NA
      if (!is.na(x$nbinsdpi)) DPI.df <- x$pdpi+1+(x$nbinsdpi-1)*(x$pdpi-x$sdpi+1)
      cat(format(sprintf("%3.0f", DPI.df), width=13, justify="right"))
      cat(format(sprintf("%3.3f", x$imse.bsq.dpi), width=15 , justify="right"))
      cat(format(sprintf("%3.3f", x$imse.var.dpi), width=15 , justify="right"))
      cat("\n")

      cat(format("DPI-UKNOT",  width=13, justify="right"))
      cat(format(sprintf("%3.0f", x$nbinsdpi.uknot), width=13 , justify="right"))
      DPI.uknot.df <- NA
      if (!is.na(x$nbinsdpi.uknot)) DPI.uknot.df <- x$pdpi.uknot+1+(x$nbinsdpi.uknot-1)*(x$pdpi.uknot-x$sdpi.uknot+1)
      cat(format(sprintf("%3.0f", DPI.uknot.df), width=13, justify="right"))
      cat(format(sprintf("%3.0f", NA), width=15 , justify="right"))
      cat(format(sprintf("%3.0f", NA), width=15 , justify="right"))
      cat("\n")
    }
    cat(paste(rep("-", 13 + 13 + 13 + 15 + 15), collapse=""))
  } else {
    cat(paste(rep("=", 13 + 13 + 13 + 13 + 15 + 15), collapse="")); cat("\n")
    cat(format("method",      width=13, justify="right"))
    cat(format("degree",      width=13, justify="right"))
    cat(format("smoothness",  width=13, justify="right"))
    cat(format("df",          width=13, justify="right"))
    cat(format("imse, bias^2",  width=15, justify="right"))
    cat(format("imse, var.",  width=15, justify="right"))
    cat("\n")
    cat(paste(rep("-", 13 + 13 + 13 + 13 + 15 + 15), collapse="")); cat("\n")
    cat(format("ROT-POLY",  width= 13, justify="right"))
    cat(format(sprintf("%3.0f", x$prot.poly), width=13 , justify="right"))
    cat(format(sprintf("%3.0f", x$srot.poly), width=13 , justify="right"))
    ROT.poly.df <- NA
    if (!is.na(x$nbinsrot.poly)) ROT.poly.df <- x$prot.poly+1+(x$nbinsrot.poly-1)*(x$prot.poly-x$srot.poly+1)
    cat(format(sprintf("%3.0f", ROT.poly.df), width=13, justify="right"))
    cat(format(sprintf("%3.3f", x$imse.bsq.rot), width=15 , justify="right"))
    cat(format(sprintf("%3.3f", x$imse.var.rot), width=15 , justify="right"))
    cat("\n")

    cat(format("ROT-REGUL",  width= 13, justify="right"))
    cat(format(sprintf("%3.0f", x$prot.regul), width=13 , justify="right"))
    cat(format(sprintf("%3.0f", x$srot.regul), width=13 , justify="right"))
    ROT.regul.df <- NA
    if (!is.na(x$nbinsrot.regul)) ROT.regul.df <- x$prot.regul+1+(x$nbinsrot.regul-1)*(x$prot.regul-x$srot.regul+1)
    cat(format(sprintf("%3.0f", ROT.regul.df), width=13, justify="right"))
    cat(format(sprintf("%3.0f", NA), width=15 , justify="right"))
    cat(format(sprintf("%3.0f", NA), width=15 , justify="right"))
    cat("\n")

    cat(format("ROT-UKNOT",  width= 13, justify="right"))
    cat(format(sprintf("%3.0f", x$prot.uknot), width=13 , justify="right"))
    cat(format(sprintf("%3.0f", x$srot.uknot), width=13 , justify="right"))
    ROT.uknot.df <- NA
    if (!is.na(x$nbinsrot.uknot)) ROT.uknot.df <- x$prot.uknot+1+(x$nbinsrot.uknot-1)*(x$prot.uknot-x$srot.uknot+1)
    cat(format(sprintf("%3.0f", ROT.uknot.df), width=13, justify="right"))
    cat(format(sprintf("%3.0f", NA), width=15 , justify="right"))
    cat(format(sprintf("%3.0f", NA), width=15 , justify="right"))
    cat("\n")

    if (x$opt$binsmethod=="IMSE direct plug-in (select degree and smoothness)") {
      cat(format("DPI",  width=13, justify="right"))
      cat(format(sprintf("%3.0f", x$pdpi), width=13 , justify="right"))
      cat(format(sprintf("%3.0f", x$sdpi), width=13 , justify="right"))
      DPI.df <- NA
      if (!is.na(x$nbinsdpi)) DPI.df <- x$pdpi+1+(x$nbinsdpi-1)*(x$pdpi-x$sdpi+1)
      cat(format(sprintf("%3.0f", DPI.df), width=13, justify="right"))
      cat(format(sprintf("%3.3f", x$imse.bsq.dpi), width=15 , justify="right"))
      cat(format(sprintf("%3.3f", x$imse.var.dpi), width=15 , justify="right"))
      cat("\n")

      cat(format("DPI-UKNOT",  width=13, justify="right"))
      cat(format(sprintf("%3.0f", x$pdpi.uknot), width=13 , justify="right"))
      cat(format(sprintf("%3.0f", x$sdpi.uknot), width=13 , justify="right"))
      DPI.uknot.df <- NA
      if (!is.na(x$nbinsdpi.uknot)) DPI.uknot.df <- x$pdpi.uknot+1+(x$nbinsdpi.uknot-1)*(x$pdpi.uknot-x$sdpi.uknot+1)
      cat(format(sprintf("%3.0f", DPI.uknot.df), width=13, justify="right"))
      cat(format(sprintf("%3.0f", NA), width=15 , justify="right"))
      cat(format(sprintf("%3.0f", NA), width=15 , justify="right"))
      cat("\n")
    }
    cat(paste(rep("-", 13 + 13 + 13 + 13 + 15 + 15), collapse=""))
  }
  cat("\n")
}

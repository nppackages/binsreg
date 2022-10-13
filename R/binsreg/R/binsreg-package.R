######################################################################################################
#'@title Binsreg Package Document
#'@description  Binscatter provides a flexible, yet parsimonious way of visualizing and summarizing large data sets
#'              and has been a popular methodology in applied microeconomics and other social sciences. The binsreg package provides tools for
#'              statistical analysis using the binscatter methods developed in
#'              \href{https://arxiv.org/abs/1902.09608}{Cattaneo, Crump, Farrell and Feng (2022a)}.
#'              \code{\link{binsreg}} implements binscatter least squares regression with robust inference and plots, including
#'              curve estimation, pointwise confidence intervals and uniform confidence band.
#'              \code{\link{binsqreg}} implements binscatter quantile regression with robust inference and plots, including
#'              curve estimation, pointwise confidence intervals and uniform confidence band.
#'              \code{\link{binsglm}} implements binscatter generalized linear regression with robust inference and plots, including
#'              curve estimation, pointwise confidence intervals and uniform confidence band.
#'              \code{\link{binstest}} implements binscatter-based hypothesis testing procedures for parametric specifications
#'              of and shape restrictions on the unknown function of interest.
#'              \code{\link{binspwc}} implements hypothesis testing procedures for pairwise group comparison of binscatter estimators.
#'              \code{\link{binsregselect}} implements data-driven number of bins selectors for binscatter
#'              implementation using either quantile-spaced or evenly-spaced binning/partitioning.
#'              All the commands allow for covariate adjustment, smoothness restrictions, and clustering,
#'              among other features.
#'
#'              The companion software article,
#'              \href{https://arxiv.org/abs/1902.09615}{Cattaneo, Crump, Farrell and Feng (2022b)},
#'              provides further implementation details and empirical illustration. For related Stata, R and Python packages
#'              useful for nonparametric data analysis and statistical inference, visit
#'              \href{https://nppackages.github.io/}{https://nppackages.github.io/}.
#'@importFrom stats complete.cases quantile weighted.mean sd rnorm qnorm dnorm lm glm vcov gaussian model.matrix terms.formula runif
#'@importFrom splines splineDesign
#'@importFrom sandwich vcovCL
#'@importFrom matrixStats colWeightedMeans colWeightedMedians
#'@importFrom quantreg rq summary.rq
#'@import ggplot2
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
#' Cattaneo, M. D., R. K. Crump, M. H. Farrell, and Y. Feng. 2022a: \href{https://arxiv.org/abs/1902.09608}{On Binscatter}. Working Paper.
#'
#' Cattaneo, M. D., R. K. Crump, M. H. Farrell, and Y. Feng. 2022b: \href{https://arxiv.org/abs/1902.09615}{Binscatter Regressions}. Working Paper.
#'
#'@aliases binsreg-package
"_PACKAGE"

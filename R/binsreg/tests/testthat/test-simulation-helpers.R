legacy_binsreg_pval <- function(num, denom, reps, tstat=NULL, side=NULL, alpha, lp=Inf) {
  tvec <- c()
  pval <- NA
  if (!is.null(tstat)) pval <- numeric(nrow(tstat))
  cval <- NA
  k <- ncol(num)

  for (i in seq_len(reps)) {
    eps <- matrix(rnorm(k, 0, 1), ncol = 1)
    tx <- (num %*% eps) / denom

    if (!is.null(side)) {
      if (side == "two") {
        if (is.infinite(lp)) tvec[i] <- max(abs(tx))
        else                 tvec[i] <- mean(abs(tx)^lp)^(1/lp)
      } else if (side == "left") {
        tvec[i] <- max(tx)
      } else if (side == "right") {
        tvec[i] <- min(tx)
      }
    }

    if (!is.null(tstat)) {
      for (j in seq_len(nrow(tstat))) {
        if (tstat[j, 2] == 1) {
          pval[j] <- pval[j] + (max(tx) >= tstat[j, 1])
        } else if (tstat[j, 2] == 2) {
          pval[j] <- pval[j] + (min(tx) <= tstat[j, 1])
        } else if (tstat[j, 2] == 3) {
          if (is.infinite(lp)) pval[j] <- pval[j] + (max(abs(tx)) >= tstat[j, 1])
          else                 pval[j] <- pval[j] + (mean(abs(tx)^lp)^(1/lp) >= tstat[j, 1])
        }
      }
    }
  }

  if (!is.null(tstat)) pval <- pval / reps
  if (!is.null(side)) cval <- quantile(tvec, alpha/100, na.rm=TRUE, names=FALSE, type=2)
  list(pval=pval, cval=cval)
}

legacy_binspwc_pval <- function(nummat1, nummat2, denom1, denom2, reps, tstat=NULL, testtype=NULL, lp=Inf, alpha=95) {
  pval <- 0
  tvec <- c()
  k1 <- ncol(nummat1)
  k2 <- ncol(nummat2)

  for (i in seq_len(reps)) {
    eps1 <- matrix(rnorm(k1, 0, 1), ncol = 1)
    eps2 <- matrix(rnorm(k2, 0, 1), ncol = 1)
    tx <- (nummat1 %*% eps1 - nummat2 %*% eps2) / sqrt(denom1^2+denom2^2)

    if (testtype == "left") {
      pval <- pval + (max(tx) >= tstat)
    } else if (testtype == "right") {
      pval <- pval + (min(tx) <= tstat)
    } else {
      if (is.infinite(lp)) pval <- pval + (max(abs(tx)) >= tstat)
      else                 pval <- pval + (mean(abs(tx)^lp)^(1/lp) >= tstat)
    }

    tvec[i] <- max(abs(tx))
  }

  list(
    pval=pval/reps,
    cval.cb=quantile(tvec, alpha/100, na.rm=TRUE, names=FALSE, type=2)
  )
}

test_that("binsreg.pval preserves the legacy simulation stream", {
  num <- matrix(c(
    0.4, -0.2, 0.7,
    1.1, 0.5, -0.3,
    -0.6, 0.9, 0.2,
    0.3, -0.8, 0.6
  ), nrow=4, byrow=TRUE)
  denom <- c(1.2, 0.9, 1.5, 1.1)
  tstat <- matrix(c(0.4, 1, -0.3, 2, 0.8, 3), ncol=2, byrow=TRUE)

  set.seed(60201)
  expected <- legacy_binsreg_pval(num, denom, 37, tstat=tstat, side="two", alpha=95, lp=Inf)
  set.seed(60201)
  observed <- binsreg:::binsreg.pval(num, denom, 37, tstat=tstat, side="two", alpha=95, lp=Inf)

  expect_equal(observed, expected, tolerance=1e-12)
})

test_that("binspwc.pval preserves the legacy simulation stream", {
  nummat1 <- matrix(c(
    0.1, 0.3, -0.5,
    0.7, -0.2, 0.4,
    -0.6, 0.8, 0.2,
    0.4, -0.1, 0.9
  ), nrow=4, byrow=TRUE)
  nummat2 <- matrix(c(
    -0.2, 0.6,
    0.5, -0.4,
    0.8, 0.1,
    -0.3, 0.7
  ), nrow=4, byrow=TRUE)
  denom1 <- c(1.0, 1.3, 0.8, 1.1)
  denom2 <- c(0.9, 0.7, 1.2, 1.4)

  set.seed(60202)
  expected <- legacy_binspwc_pval(nummat1, nummat2, denom1, denom2, 41, tstat=0.75, testtype="two-sided", alpha=95)
  set.seed(60202)
  observed <- binsreg:::binspwc.pval(nummat1, nummat2, denom1, denom2, 41, tstat=0.75, testtype="two-sided", alpha=95)

  expect_equal(observed, expected, tolerance=1e-12)
})

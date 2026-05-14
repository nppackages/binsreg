test_that("binsregselect cache reuses deterministic selections only", {
  dat <- make_binsreg_regression_data()
  w <- as.matrix(dat[, c("w1", "w2")])
  binsregselect.cache.env$entries <- list()

  fit1 <- quiet_binsreg(
    binsregselect.cached(
      dat$y, dat$x, w,
      bins = c(2, 2), nbins = 8,
      vce = "HC1", masspoints = "off", weights = dat$weights
    )
  )
  expect_length(binsregselect.cache.env$entries, 1)

  fit2 <- quiet_binsreg(
    binsregselect.cached(
      dat$y, dat$x, w,
      bins = c(2, 2), nbins = 8,
      vce = "HC1", masspoints = "off", weights = dat$weights
    )
  )
  expect_length(binsregselect.cache.env$entries, 1)
  expect_equal(
    unname(c(fit2$nbinsrot.regul, fit2$nbinsdpi,
             fit2$imse.bsq.rot, fit2$imse.var.rot,
             fit2$imse.bsq.dpi, fit2$imse.var.dpi)),
    unname(c(fit1$nbinsrot.regul, fit1$nbinsdpi,
             fit1$imse.bsq.rot, fit1$imse.var.rot,
             fit1$imse.bsq.dpi, fit1$imse.var.dpi)),
    tolerance = 1e-12
  )

  binsregselect.cache.env$entries <- list()
  set.seed(123)
  quiet_binsreg(
    binsregselect.cached(
      dat$y, dat$x, w,
      bins = c(2, 2), nbins = 8, randcut = 1,
      vce = "HC1", masspoints = "off", weights = dat$weights
    )
  )
  expect_length(binsregselect.cache.env$entries, 0)
})

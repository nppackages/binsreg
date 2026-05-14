test_that("binsregselect numerical output is stable", {
  dat <- make_binsreg_regression_data()
  w <- as.matrix(dat[, c("w1", "w2")])

  fit <- quiet_binsreg(
    binsregselect(
      dat$y, dat$x, w,
      bins = c(2, 2), nbins = 8,
      vce = "HC1", masspoints = "off", weights = dat$weights
    )
  )

  observed <- unname(c(
    fit$nbinsrot.regul, fit$nbinsdpi,
    fit$imse.bsq.rot, fit$imse.var.rot,
    fit$imse.bsq.dpi, fit$imse.var.dpi
  ))
  expected <- c(
    6, 6, 9.30877220225128, 0.0865275353823115,
    9.73806156422583, 0.0467717841532798
  )

  expect_equal(observed, expected, tolerance = 1e-8)
})

test_that("binsreg numerical output is stable", {
  dat <- make_binsreg_regression_data()
  w <- as.matrix(dat[, c("w1", "w2")])

  fit <- quiet_binsreg(
    binsreg(
      dat$y, dat$x, w,
      nbins = 8, dots = c(2, 2), line = c(2, 2),
      ci = c(3, 3), cb = c(3, 3),
      nsims = 99, simsgrid = 12, simsseed = 12345,
      vce = "HC1", masspoints = "off", weights = dat$weights
    )
  )

  expect_equal(unname(c(fit$opt$nbins.by, fit$cval.by)),
               c(8, 3.02049419884025), tolerance = 1e-8)
  expect_equal(first_numeric(fit$data.plot[[1]]$data.dots$fit),
               c(0.187508894588, -0.218167236007, -0.607516609981,
                 0.239242004003, 1.02611170345),
               tolerance = 1e-8)
  expect_equal(first_numeric(fit$data.plot[[1]]$data.line$fit),
               c(0.442402001168, 0.410211807802, 0.379074686922,
                 0.348990638528, 0.319959662619),
               tolerance = 1e-8)
  expect_equal(first_numeric(c(fit$data.plot[[1]]$data.ci$ci.l,
                               fit$data.plot[[1]]$data.ci$ci.r)),
               c(0.00187366873839, -0.31300172053, -0.651158005877,
                 0.181366404978, 0.933063006236),
               tolerance = 1e-8)
  expect_equal(first_numeric(c(fit$data.plot[[1]]$data.cb$cb.l,
                               fit$data.plot[[1]]$data.cb$cb.r)),
               c(0.307996996503, 0.275249053045, 0.221095144575,
                 0.159508431889, 0.100453972728),
               tolerance = 1e-8)
})

test_that("plot data report raw bin counts", {
  x <- seq(0, 1, length.out = 90)
  y <- 1 + 2 * x + 0.1 * sin(seq_along(x))
  expected <- c(30, 30, 30)

  check_counts <- function(fit) {
    dat <- fit$data.plot[[1]]
    expect_equal(dat$data.bin$n, expected)
    expect_equal(dat$data.dots$n, expected[dat$data.dots$bin])
  }

  check_counts(quiet_binsreg(
    binsreg(y, x, nbins = 3, binspos = "es", dots = c(0, 0), line = c(0, 0),
            ci = NULL, cb = NULL, masspoints = "off", printplot = FALSE)
  ))
  check_counts(quiet_binsreg(
    binsqreg(y, x, nbins = 3, binspos = "es", dots = c(0, 0), line = c(0, 0),
             ci = NULL, cb = NULL, masspoints = "off", printplot = FALSE)
  ))
  check_counts(quiet_binsreg(
    binsglm(y, x, nbins = 3, binspos = "es", dots = c(0, 0), line = c(0, 0),
            ci = NULL, cb = NULL, masspoints = "off", printplot = FALSE)
  ))
})

test_that("noplot suppresses rendering but keeps requested polynomial data", {
  x <- seq(0, 1, length.out = 90)
  y <- 1 + 2 * x + 0.1 * sin(seq_along(x))

  check_poly <- function(fit) {
    expect_null(fit$bins_plot)
    expect_s3_class(fit$data.plot[[1]]$data.poly, "data.frame")
    expect_gt(nrow(fit$data.plot[[1]]$data.poly), 0)
    expect_true(all(c("x", "fit", "n") %in% names(fit$data.plot[[1]]$data.poly)))
  }

  check_poly(quiet_binsreg(
    binsreg(y, x, nbins = 3, binspos = "es", dots = c(0, 0), line = NULL,
            ci = NULL, cb = NULL, polyreg = 2, polyreggrid = 4,
            masspoints = "off", noplot = TRUE, printplot = FALSE)
  ))
  check_poly(quiet_binsreg(
    binsqreg(y, x, nbins = 3, binspos = "es", dots = c(0, 0), line = NULL,
             ci = NULL, cb = NULL, polyreg = 2, polyreggrid = 4,
             masspoints = "off", noplot = TRUE, printplot = FALSE)
  ))
  check_poly(quiet_binsreg(
    binsglm(y, x, nbins = 3, binspos = "es", dots = c(0, 0), line = NULL,
            ci = NULL, cb = NULL, polyreg = 2, polyreggrid = 4,
            masspoints = "off", noplot = TRUE, printplot = FALSE)
  ))
})

test_that("binsqreg numerical output is stable", {
  dat <- make_binsreg_regression_data()
  w <- as.matrix(dat[, c("w1", "w2")])

  fit <- quiet_binsreg(
    binsqreg(
      dat$y, dat$x, w,
      nbins = 8, quantile = 0.5,
      dots = c(1, 1), line = c(1, 1), ci = c(2, 2), cb = c(2, 2),
      nsims = 49, simsgrid = 10, simsseed = 222,
      vce = "nid", masspoints = "off", weights = dat$weights
    )
  )

  expect_equal(unname(c(fit$opt$nbins.by, fit$cval.by)),
               c(8, 3.15318486770566), tolerance = 1e-7)
  expect_equal(first_numeric(fit$data.plot[[1]]$data.dots$fit),
               c(0.274021864802, -0.219841695995, -0.471057999428,
                 0.269198628547, 0.993859900129),
               tolerance = 1e-7)
  expect_equal(first_numeric(fit$data.plot[[1]]$data.line$fit),
               c(0.427124764911, 0.410361895565, 0.393599026219,
                 0.376836156872, 0.360073287526),
               tolerance = 1e-7)
  expect_equal(first_numeric(c(fit$data.plot[[1]]$data.ci$ci.l,
                               fit$data.plot[[1]]$data.ci$ci.r)),
               c(0.0994554818278, -0.289821156662, -0.736107499285,
                 0.143469258264, 0.886694108767),
               tolerance = 1e-7)
  expect_equal(first_numeric(c(fit$data.plot[[1]]$data.cb$cb.l,
                               fit$data.plot[[1]]$data.cb$cb.r)),
               c(0.289296378626, 0.273190329119, 0.250627625314,
                 0.222800459969, 0.191620038397),
               tolerance = 1e-7)
})

test_that("polynomial confidence interval covariance output is stable", {
  dat <- make_binsreg_regression_data()
  w <- as.matrix(dat[, c("w1", "w2")])

  fit <- quiet_binsreg(
    binsreg(
      dat$y, dat$x, w,
      nbins = 8, dots = NULL, line = NULL, ci = NULL, cb = NULL,
      polyreg = 2, polyreggrid = 4, polyregcigrid = 5,
      vce = "HC1", masspoints = "off", weights = dat$weights,
      cluster = dat$cluster
    )
  )
  expect_equal(first_numeric(fit$data.plot[[1]]$data.polyci$polyci.l, 10),
               c(-0.703613663739, -0.657489257774, -0.612530809347,
                 -0.568776798099, -0.526268733415, -0.485050573845,
                 -0.430741435731, -0.379055061242, -0.330093600756,
                 -0.283929805983),
               tolerance = 1e-8)

  fitq <- quiet_binsreg(
    binsqreg(
      dat$y, dat$x, w,
      nbins = 8, quantile = 0.5,
      dots = NULL, line = NULL, ci = NULL, cb = NULL,
      polyreg = 2, polyreggrid = 4, polyregcigrid = 5,
      vce = "nid", masspoints = "off", weights = dat$weights
    )
  )
  expect_equal(first_numeric(fitq$data.plot[[1]]$data.polyci$polyci.l, 10),
               c(-0.828295955029, -0.781617087682, -0.735702432793,
                 -0.690563408445, -0.646212162404, -0.602661505349,
                 -0.544302644923, -0.487507882146, -0.432312500983,
                 -0.378749709173),
               tolerance = 1e-8)
})

test_that("binsglm numerical output is stable", {
  dat <- make_binsreg_regression_data()
  w <- as.matrix(dat[, c("w1", "w2")])

  fit <- quiet_binsreg(
    binsglm(
      dat$d, dat$x, w,
      nbins = 8, family = binomial(),
      dots = c(1, 1), line = c(1, 1), ci = c(2, 2), cb = c(2, 2),
      nsims = 49, simsgrid = 10, simsseed = 333,
      vce = "HC1", masspoints = "off", weights = dat$weights
    )
  )

  expect_equal(unname(c(fit$opt$nbins.by, fit$cval.by)),
               c(8, 2.94050655661042), tolerance = 1e-7)
  expect_equal(first_numeric(fit$data.plot[[1]]$data.dots$fit),
               c(0.243763950545, 0.242384831832, 0.330227329688,
                 0.330660281857, 0.569754463839),
               tolerance = 1e-7)
  expect_equal(first_numeric(fit$data.plot[[1]]$data.line$fit),
               c(0.437922356371, 0.414297617092, 0.391060652878,
                 0.368307436235, 0.346125587107),
               tolerance = 1e-7)
  expect_equal(first_numeric(c(fit$data.plot[[1]]$data.ci$ci.l,
                               fit$data.plot[[1]]$data.ci$ci.r)),
               c(-0.0666702551551, 0.0551875640847, 0.00978450159187,
                 0.0966651352068, 0.243401806624),
               tolerance = 1e-7)
  expect_equal(first_numeric(c(fit$data.plot[[1]]$data.cb$cb.l,
                               fit$data.plot[[1]]$data.cb$cb.r)),
               c(-0.1489606982, -0.123317365498, -0.110441862323,
                 -0.110735637909, -0.120523418081),
               tolerance = 1e-7)
})

test_that("binstest numerical output is stable", {
  dat <- make_binsreg_regression_data()
  w <- as.matrix(dat[, c("w1", "w2")])

  fit <- quiet_binsreg(
    binstest(
      dat$y, dat$x, w,
      nbins = 8, testmodelpoly = 1, testmodel = c(2, 2),
      testshape = c(2, 2), testshapel = -1, testshaper = 1,
      testshape2 = 0, nsims = 99, simsgrid = 12, simsseed = 444,
      vce = "HC1", masspoints = "off", weights = dat$weights
    )
  )

  observed <- c(
    fit$testshapeL$stat.shapeL, fit$testshapeL$pval.shapeL,
    fit$testshapeR$stat.shapeR, fit$testshapeR$pval.shapeR,
    fit$testshape2$stat.shape2, fit$testshape2$pval.shape2,
    fit$testpoly$stat.poly, fit$testpoly$pval.poly,
    fit$testmodel$stat.model, fit$testmodel$pval.model
  )
  expected <- c(
    75.6289600143633, 0, -58.0583069369764, 0,
    40.7767733786001, 0, 20.2356974263774, 0, NA, NA
  )

  expect_equal(observed, expected, tolerance = 1e-8)
})

test_that("binspwc numerical output is stable", {
  dat <- make_binsreg_regression_data()
  w <- as.matrix(dat[, c("w1", "w2")])

  fit <- quiet_binsreg(
    binspwc(
      dat$y, dat$x, w,
      by = dat$by, bynbins = 8, pselect = 1, sselect = 1,
      pwc = c(2, 2), testtype = "two-sided",
      nsims = 99, simsgrid = 12, simsseed = 555,
      vce = "HC1", masspoints = "off", weights = dat$weights
    )
  )

  observed <- c(fit$tstat[, 1], fit$pval[, 1], fit$cval.cb[, 1])
  expected <- c(3.08417386477364, 0.0101010101010101, 2.65048187418525)

  expect_equal(observed, expected, tolerance = 1e-8)
})

test_that("plot objects can be returned without printing", {
  dat <- make_binsreg_regression_data()
  w <- as.matrix(dat[, c("w1", "w2")])

  reg <- quiet_binsreg(
    binsreg(
      dat$y, dat$x, w,
      nbins = 8, masspoints = "off", weights = dat$weights,
      printplot = FALSE
    )
  )
  expect_true(inherits(reg$bins_plot, "ggplot"))
  expect_true(!is.null(reg$data.plot[[1]]$data.dots))

  qreg <- quiet_binsreg(
    binsqreg(
      dat$y, dat$x, w,
      nbins = 8, quantile = 0.5, masspoints = "off",
      weights = dat$weights, printplot = FALSE
    )
  )
  expect_true(inherits(qreg$bins_plot, "ggplot"))
  expect_true(!is.null(qreg$data.plot[[1]]$data.dots))

  glm <- quiet_binsreg(
    binsglm(
      dat$d, dat$x, w,
      nbins = 8, family = binomial(), masspoints = "off",
      weights = dat$weights, printplot = FALSE
    )
  )
  expect_true(inherits(glm$bins_plot, "ggplot"))
  expect_true(!is.null(glm$data.plot[[1]]$data.dots))

  pwc <- quiet_binsreg(
    binspwc(
      dat$y, dat$x, w,
      by = dat$by, bynbins = 8, pwc = c(1, 1), plot = TRUE,
      nsims = 19, simsgrid = 5, simsseed = 123,
      vce = "HC1", masspoints = "off", weights = dat$weights,
      printplot = FALSE
    )
  )
  expect_true(inherits(pwc$bins_plot, "ggplot"))
  expect_true(!is.null(pwc$data.plot[[1]]$data.cb))
})

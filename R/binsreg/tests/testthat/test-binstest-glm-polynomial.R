test_that("GLM polynomial tests compare fits on the requested scale", {
  set.seed(20260820)
  n <- 2000
  x <- runif(n, -1, 1)
  w <- cbind(rnorm(n), rnorm(n))
  index <- -0.4 + 1.1*x + 0.35*w[, 1] - 0.2*w[, 2]
  y <- rbinom(n, 1, plogis(index))

  capture.output(
    response.level <- suppressWarnings(binstest(
      y, x, w=w, family=binomial(), testmodelpoly=1,
      nbins=16, testmodel=c(1, 1), nsims=99, simsgrid=20,
      simsseed=17, masspoints="off"
    ))
  )
  capture.output(
    index.level <- suppressWarnings(binstest(
      y, x, w=w, family=binomial(), nolink=TRUE, testmodelpoly=1,
      nbins=16, testmodel=c(1, 1), nsims=99, simsgrid=20,
      simsseed=17, masspoints="off"
    ))
  )
  capture.output(
    response.derivative <- suppressWarnings(binstest(
      y, x, w=w, family=binomial(), deriv=1, testmodelpoly=1,
      nbins=16, testmodel=c(1, 1), nsims=99, simsgrid=20,
      simsseed=17, masspoints="off"
    ))
  )

  expect_true(is.finite(response.level$testpoly$stat.poly))
  expect_true(is.finite(response.derivative$testpoly$stat.poly))
  expect_lt(response.level$testpoly$stat.poly, 10)
  expect_lt(response.derivative$testpoly$stat.poly, 10)
  expect_lt(abs(response.level$testpoly$stat.poly - index.level$testpoly$stat.poly), 2)
})

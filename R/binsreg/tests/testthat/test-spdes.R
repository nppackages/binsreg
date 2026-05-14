legacy_spdes_s0 <- function(eval, p, knot, deriv) {
  pos <- findInterval(eval, knot, rightmost.closed = TRUE, left.open = TRUE)
  polyx <- matrix(0, length(eval), p + 1)
  h <- diff(knot)
  eval.cen <- (eval - knot[-length(knot)][pos]) / h[pos]
  for (j in (deriv + 1):(p + 1)) {
    polyx[, j] <- eval.cen^(j - 1 - deriv) *
      factorial(j - 1) / factorial(j - 1 - deriv) / h[pos]^deriv
  }
  matrix(sapply(seq_len(length(knot) - 1), function(i) polyx * (pos == i)),
         nrow = length(eval))
}

test_that("s=0 spline design preserves legacy column layout", {
  knot <- c(-1, -0.25, 0.1, 0.7, 1.5)
  eval <- c(-0.95, -0.4, -0.1, 0.2, 0.65, 0.9, 1.4)

  for (p in 0:4) {
    for (deriv in 0:p) {
      expected <- legacy_spdes_s0(eval, p, knot, deriv)
      observed <- binsreg:::binsreg.spdes(eval, p=p, s=0, knot=knot, deriv=deriv)
      expect_equal(observed, expected, tolerance=0)
    }
  }
})

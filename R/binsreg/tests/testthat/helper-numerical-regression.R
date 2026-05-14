make_binsreg_regression_data <- function() {
  set.seed(20240513)
  n <- 160
  x <- sort(runif(n, -1, 1))
  w1 <- 0.5 * cos(seq_len(n) / 9) + rnorm(n, sd = 0.2)
  w2 <- as.numeric(seq_len(n) %% 3 == 0)
  pr <- plogis(-0.25 + 0.8 * x + 0.4 * w1 - 0.2 * w2)
  d <- rbinom(n, 1, pr)
  y <- 0.5 + sin(pi * x) + 0.3 * w1 - 0.15 * w2 + rnorm(n, sd = 0.15)
  by <- rep(c("a", "b"), length.out = n)
  weights <- 1 + (seq_len(n) %% 5) / 10
  cluster <- rep(seq_len(40), each = 4)
  data.frame(
    y = y, x = x, w1 = w1, w2 = w2, d = d, by = by,
    weights = weights, cluster = cluster
  )
}

quiet_binsreg <- function(expr) {
  result <- NULL
  invisible(capture.output({
    result <- suppressWarnings(suppressMessages(force(expr)))
  }))
  result
}

first_numeric <- function(x, n = 5) {
  unname(as.numeric(head(x, n)))
}

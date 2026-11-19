test_that("RemoveSeasonal() returns seasonal means and anomalies", {
  y <- c(1, 3, 10, 12, 100)
  t_mmdd <- c(101, 101, 102, 102, 103)
  out <- RemoveSeasonal(y, t_mmdd)

  expect_equal(out$mu, c("101" = 2, "102" = 11, "103" = 100))
  expect_equal(out$anomaly, c(-1, 1, -1, 1, 0))
})

test_that("RemoveSeasonal() removes the climatology exactly", {
  set.seed(42)
  t_mmdd <- rep(1:365, times = 3)
  y <- rep(rnorm(365), times = 3) + rnorm(365 * 3, sd = 0.01)
  out <- RemoveSeasonal(y, t_mmdd)

  expect_lt(max(abs(tapply(out$anomaly, t_mmdd, mean))), 1e-12)
})

test_that("RemoveSeasonal() keeps missing values and checks lengths", {
  y <- c(1, NA, 3, 5, 7, 9)
  t_mmdd <- c(101, 101, 101, 102, 102, 102)
  out <- RemoveSeasonal(y, t_mmdd)

  expect_equal(unname(out$mu), c(2, 7))
  expect_true(is.na(out$anomaly[2]))
  expect_error(RemoveSeasonal(1:3, 1:2), "same length")
})

test_that("LS_OneBreak() recovers a known piecewise model", {
  x <- 1:20
  y <- 1 + 2 * (x - mean(x)) + 5 * (x > 10)
  fit <- LS_OneBreak(y, x, 10)

  expect_equal(fit$coef, c(1, 2, 5))
  expect_equal(fit$fitted, y)
  expect_equal(fit$SSE, 0)
})

test_that("LS_OneBreak() without a break is ordinary linear regression", {
  x <- 1:20
  y <- 1 + 2 * (x - mean(x))
  fit <- LS_OneBreak(y, x, NA)
  ref <- stats::lm(y ~ I(x - mean(x)))

  expect_equal(fit$coef, unname(stats::coef(ref)))
})

test_that("LS_OneBreak() equals LS_MultiBreaks() with one break", {
  set.seed(1)
  x <- 1:120
  y <- 3 + 0.1 * x + ifelse(x > 60, 2, 0) + rnorm(120, sd = 0.1)
  one <- LS_OneBreak(y, x, 60)
  multi <- LS_MultiBreaks(y, x, c(60, 120))

  expect_equal(one$coef, multi$coef)
  expect_equal(one$fitted, multi$fitted)
})

test_that("LS_MultiBreaks() fits three segments", {
  x <- 1:90
  y <- c(rep(0, 30), rep(5, 30), rep(-2, 30)) + 0.1 * (x - mean(x))
  fit <- LS_MultiBreaks(y, x, c(30, 60, 90))

  expect_equal(fit$coef, c(0, 0.1, 5, -2))
  expect_equal(fit$SSE, 0)
})

test_that("LS_MultiBreaks() reduces SSE when adding a break", {
  set.seed(2)
  x <- 1:100
  y <- 0.5 * x + ifelse(x > 50, 3, 0) + rnorm(100, sd = 0.2)
  sse0 <- LS_MultiBreaks(y, x, 100, only.SSE = TRUE)$SSE
  sse1 <- LS_MultiBreaks(y, x, c(50, 100), only.SSE = TRUE)$SSE

  expect_lt(sse1, sse0)
})

test_that("piecewise regression validates its input", {
  expect_error(LS_OneBreak(1:10, 1:9, 5), "same length")
  expect_error(LS_OneBreak(1:10, 1:10, 0), "Ic")
  expect_error(LS_OneBreak(1:10, 1:10, 10), "Ic")
  expect_error(LS_MultiBreaks(1:10, 1:10, c(5, 9)), "Ips")
  expect_error(LS_MultiBreaks(1:10, 1:10, c(5, 5, 10)), "Ips")
  expect_error(LS_OneBreak(c(1, NA, 3), 1:3, 1), "missing")
})

test_that("lm_solve() returns the SSE only when asked", {
  fit <- lm_solve(cbind(1, 1:5), c(2, 4, 6, 8, 10), only.SSE = TRUE)

  expect_named(fit, "SSE")
  expect_equal(fit$SSE, 0)
})

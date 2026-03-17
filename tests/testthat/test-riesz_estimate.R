test_that("riesz_estimate errors for invalid rc input", {
  df <- data.frame(Y = 1:3, A = c(1, 0, 1))
  expect_error(riesz_estimate(df, rc = 1), "must be a RieszCurve or ComposedRieszCurve")
})

test_that("riesz_estimate matches manual calculations from fit_ic", {
  df <- data.frame(
    Y = c(1, 2, 4),
    A = c(1, 0, 1)
  )

  rc <- RieszCurve$new(
    nuis = list(
      g = c(0.5, 0.25, 0.8),
      m = c(1.2, 1.7, 3.5),
      ma = c(2.0, 2.0, 2.0)
    ),
    alpha = ~ A / g,
    f = ~ m,
    h = ~ ma,
    ic_expr = ~ h + alpha * (Y - f)
  )

  out <- riesz_estimate(df, rc)
  rc$fit(df)
  phi <- rc$fit_ic

  expect_equal(out$estimate, mean(phi))
  expect_equal(out$var, stats::var(phi) / nrow(df))
  expect_equal(out$se, sqrt(stats::var(phi) / nrow(df)))

  expected_ci_low <- mean(phi) + sqrt(stats::var(phi) / nrow(df)) * stats::qnorm(0.05 / 2)
  expected_ci_high <- mean(phi) + sqrt(stats::var(phi) / nrow(df)) * stats::qnorm(1 - 0.05 / 2)

  expect_equal(out$ci_low, expected_ci_low)
  expect_equal(out$ci_high, expected_ci_high)
})

test_that("catalog cfmean_a returns a RieszCurve", {
  rc <- riesz_representer_catalog$cfmean_a(
    L = "L",
    A = "A",
    a = 1,
    Y = "Y",
    nuis = function(data) {
      list(
        m = rep(0, nrow(data)),
        ma = rep(1, nrow(data)),
        g = rep(0.5, nrow(data))
      )
    }
  )

  expect_true(inherits(rc, "RieszCurve"))
})

test_that("catalog cfmean_a constructs correct alpha", {
  df <- data.frame(
    L = c(0, 1, 2),
    A = c(1, 0, 1),
    Y = c(3, 4, 5)
  )

  rc <- riesz_representer_catalog$cfmean_a(
    L = "L",
    A = "A",
    a = 1,
    Y = "Y",
    nuis = function(data) {
      list(
        m = c(0, 0, 0),
        ma = c(10, 10, 10),
        g = c(0.2, 0.4, 0.5)
      )
    }
  )

  rc$fit(df)
  expect_equal(rc$fit_alpha, c(1 / 0.2, 0, 1 / 0.5))
})

test_that("riesz_estimate for cfmean_a is close to truth in a large sample for a = 1", {
  set.seed(1)
  n <- 8000

  df <- data.frame(
    L = rnorm(n)
  )
  df$A <- rbinom(n, size = 1, prob = plogis(df$L))
  df$Y <- 2 + df$L + 3 * df$A + rnorm(n, sd = 0.1)

  rc <- riesz_representer_catalog$cfmean_a(
    L = "L",
    A = "A",
    a = 1,
    Y = "Y",
    nuis = function(data) {
      mfit <- stats::lm(Y ~ L + A, data = data)
      data1 <- data
      data1$A <- 1
      gfit <- stats::glm(I(A == 1) ~ L, family = stats::binomial(), data = data)

      list(
        m = stats::predict(mfit, newdata = data),
        ma = stats::predict(mfit, newdata = data1),
        g = stats::predict(gfit, newdata = data, type = "response")
      )
    }
  )

  out <- riesz_estimate(df, rc)
  expect_equal(out$estimate, 5, tolerance = 0.1)
})

test_that("riesz_estimate for cfmean_a is close to truth in a large sample for a = 0", {
  set.seed(1)
  n <- 8000

  df <- data.frame(
    L = rnorm(n)
  )
  df$A <- rbinom(n, size = 1, prob = plogis(df$L))
  df$Y <- 2 + df$L + 3 * df$A + rnorm(n, sd = 0.1)

  rc <- riesz_representer_catalog$cfmean_a(
    L = "L",
    A = "A",
    a = 0,
    Y = "Y",
    nuis = function(data) {
      mfit <- stats::lm(Y ~ L + A, data = data)
      data0 <- data
      data0$A <- 0
      gfit <- stats::glm(I(A == 0) ~ L, family = stats::binomial(), data = data)

      list(
        m = stats::predict(mfit, newdata = data),
        ma = stats::predict(mfit, newdata = data0),
        g = stats::predict(gfit, newdata = data, type = "response")
      )
    }
  )

  out <- riesz_estimate(df, rc)
  expect_equal(out$estimate, 2, tolerance = 0.1)
})

test_that("riesz_estimate for manually constructed ATE is close to truth", {
  set.seed(2)
  n <- 8000

  df <- data.frame(
    L = rnorm(n)
  )
  df$A <- rbinom(n, size = 1, prob = plogis(0.5 * df$L))
  df$Y <- 4 + 2 * df$L + 1 * df$A + rnorm(n, sd = 0.1)

  rc_ate <- RieszCurve$new(
    nuis = function(data) {
      mfit <- stats::lm(Y ~ L + A, data = data)
      gfit <- stats::glm(A ~ L, family = stats::binomial(), data = data)

      d1 <- data
      d0 <- data
      d1$A <- 1
      d0$A <- 0

      list(
        m = stats::predict(mfit, newdata = data),
        m1 = stats::predict(mfit, newdata = d1),
        m0 = stats::predict(mfit, newdata = d0),
        g = stats::predict(gfit, newdata = data, type = "response")
      )
    },
    alpha = ~ A / g + (1 - A) / (1 - g),
    f = ~ m,
    h = ~ m1 - m0,
    ic_expr = ~ h + alpha * (Y - f)
  )

  out <- riesz_estimate(df, rc_ate)
  expect_equal(out$estimate, 1, tolerance = 0.1)
})

test_that("riesz_estimate works on a simple ComposedRieszCurve", {
  df <- data.frame(
    L1 = c(0, 1, -1, 0.5),
    A1 = c(1, 0, 1, 1),
    L2 = c(0.2, 1.2, -0.5, 0.7),
    A2 = c(1, 1, 0, 1),
    Y = c(2.1, 3.0, 0.5, 2.4)
  )

  rc1 <- RieszCurve$new(
    nuis = list(
      m = c(2.0, 2.8, 0.6, 2.2),
      m1 = c(2.2, 3.1, 0.9, 2.5),
      g = c(0.8, 0.7, 0.4, 0.9)
    ),
    h = ~ m1,
    f = ~ m,
    alpha = ~ A2 / g,
    ic_expr = ~ h + alpha * (Y - f)
  )

  rc2 <- RieszCurve$new(
    nuis = function(data) {
      list(
        m = c(2.1, 2.7, 0.7, 2.3),
        m1 = c(2.3, 2.9, 1.0, 2.6),
        g = c(0.9, 0.5, 0.8, 0.7)
      )
    },
    h = ~ m1,
    f = ~ m,
    alpha = ~ A1 / g,
    ic_expr = ~ h + alpha * (`h[j=1]` - f)
  )

  rc_composed <- ComposedRieszCurve$new(rc_list = list(rc1, rc2))
  out <- riesz_estimate(df, rc_composed)

  expect_true(is.list(out))
  expect_true(all(c("estimate", "var", "se", "ci_low", "ci_high") %in% names(out)))
  expect_true(is.numeric(out$estimate))
  expect_true(length(out$estimate) == 1)
})

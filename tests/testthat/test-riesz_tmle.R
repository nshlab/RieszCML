test_that("riesz_tmle errors for invalid rc input", {
  df <- data.frame(Y = 1:3, A = c(1, 0, 1))
  expect_error(
    riesz_tmle(
      data = df,
      rc = 1,
      fluctuation_type = "identity",
      outcome_col = "Y"
    ),
    "`rc`"
  )
})

test_that("riesz_tmle errors for invalid outcome_col", {
  df <- data.frame(Y = 1:3, A = c(1, 0, 1))

  rc <- RieszCurve$new(
    nuis = list(m = c(1, 2, 3), ma = c(2, 2, 2), g = c(0.5, 0.5, 0.5)),
    alpha = ~ A / g,
    f = ~ m,
    h = ~ ma,
    ic_expr = ~ h + alpha * (Y - f)
  )

  expect_error(
    riesz_tmle(df, rc, fluctuation_type = "identity", outcome_col = ""),
    "outcome_col"
  )

  expect_error(
    riesz_tmle(df, rc, fluctuation_type = "identity", outcome_col = NA_character_),
    "outcome_col"
  )

  expect_error(
    riesz_tmle(df, rc, fluctuation_type = "identity", outcome_col = "Z"),
    "not a column"
  )
})

test_that("riesz_tmle logistic branch errors for missing or bad bounds", {
  df <- data.frame(Y = c(0.2, 0.4, 0.8), A = c(1, 0, 1))

  rc <- RieszCurve$new(
    nuis = list(m = c(0.3, 0.5, 0.7), ma = c(0.6, 0.6, 0.6), g = c(0.5, 0.5, 0.5)),
    alpha = ~ A / g,
    f = ~ m,
    h = ~ ma,
    ic_expr = ~ h + alpha * (Y - f),
    targeting_steps = list()
  )

  expect_error(
    riesz_tmle(df, rc, fluctuation_type = "logistic", outcome_col = "Y"),
  )

  expect_error(
    riesz_tmle(df, rc, fluctuation_type = "logistic", bounds = 1, outcome_col = "Y")
  )

  expect_error(
    riesz_tmle(df, rc, fluctuation_type = "logistic", bounds = c(0, 0), outcome_col = "Y")
  )

  expect_error(
    riesz_tmle(df, rc, fluctuation_type = "logistic", bounds = c(1, 0), outcome_col = "Y")
  )
})

test_that("identity riesz_tmle agrees approximately with riesz_estimate", {
  set.seed(1)
  n <- 3000

  df <- data.frame(
    L = rnorm(n)
  )
  df$A <- rbinom(n, 1, plogis(df$L))
  df$Y <- 1 + df$L + 2 * df$A + rnorm(n, sd = 0.1)

  rc <- riesz_curve_catalog$cfmean_a(
    a = 1,
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

  out_est <- riesz_estimate(df, rc)
  out_tmle <- riesz_tmle(
    data = df,
    rc = rc,
    fluctuation_type = "identity",
    outcome_col = "Y"
  )

  expect_lt(abs(out_tmle$estimate - out_est$estimate), 1e-3)
  expect_lt(abs(out_tmle$var - out_est$var), 1e-3)
  expect_lt(abs(out_tmle$se - out_est$se), 1e-3)
  expect_lt(abs(out_tmle$ci_low - out_est$ci_low), 1e-3)
  expect_lt(abs(out_tmle$ci_high - out_est$ci_high), 1e-3)
})

test_that("single-stage logistic riesz_tmle returns expected structure", {
  df <- data.frame(
    Y = c(0.15, 0.25, 0.5, 0.8),
    A = c(1, 0, 1, 1),
    L = c(-1, 0, 1, 2)
  )

  rc <- RieszCurve$new(
    nuis = list(
      m = c(0.2, 0.3, 0.45, 0.7),
      ma = c(0.55, 0.55, 0.55, 0.55),
      g = c(0.7, 0.4, 0.8, 0.9)
    ),
    alpha = ~ I(A == 1) / g,
    f = ~ m,
    h = ~ ma,
    ic_expr = ~ h + alpha * (Y - f),
    targeting_steps = list(m = list(), ma = list(A = 1))
  )

  out <- riesz_tmle(
    data = df,
    rc = rc,
    fluctuation_type = "logistic",
    bounds = c(0, 1),
    outcome_col = "Y"
  )

  expect_true(inherits(out, 'RieszFit'))
  expect_true(all(c("estimate", "var", "se", "ci_low", "ci_high", "eps", "ic", "f_star") %in% names(out)))
  expect_equal(length(out$ic), nrow(df))
  expect_equal(length(out$f_star), nrow(df))
  expect_true(is.numeric(out$eps))
  expect_equal(length(out$eps), 1)
})

test_that("single-stage logistic riesz_tmle keeps f_star within bounds", {
  df <- data.frame(
    Y = c(-0.8, -0.2, 0.1, 0.9),
    A = c(1, 0, 1, 1)
  )

  rc <- RieszCurve$new(
    nuis = list(
      m = c(-0.7, -0.1, 0.2, 0.8),
      ma = c(0.1, 0.1, 0.1, 0.1),
      g = c(0.6, 0.4, 0.7, 0.8)
    ),
    alpha = ~ I(A == 1) / g,
    f = ~ m,
    h = ~ ma,
    ic_expr = ~ h + alpha * (Y - f),
    targeting_steps = list(m = list(), ma = list(A = 1))
  )

  out <- riesz_tmle(
    data = df,
    rc = rc,
    fluctuation_type = "logistic",
    bounds = c(-1, 1),
    outcome_col = "Y"
  )

  expect_true(all(out$f_star >= -1))
  expect_true(all(out$f_star <= 1))
})

test_that("single-stage logistic riesz_tmle has near-zero epsilon when Y equals f exactly", {
  df <- data.frame(
    Y = c(0.2, 0.4, 0.7, 0.8),
    A = c(1, 0, 1, 1)
  )

  rc <- RieszCurve$new(
    nuis = list(
      m = c(0.2, 0.4, 0.7, 0.8),
      ma = c(0.5, 0.5, 0.5, 0.5),
      g = c(0.6, 0.5, 0.7, 0.9)
    ),
    alpha = ~ I(A == 1) / g,
    f = ~ m,
    h = ~ ma,
    ic_expr = ~ h + alpha * (Y - f),
    targeting_steps = list(m = list(), ma = list(A = 1))
  )

  out <- riesz_tmle(
    data = df,
    rc = rc,
    fluctuation_type = "logistic",
    bounds = c(0, 1),
    outcome_col = "Y"
  )

  expect_equal(out$eps, 0, tolerance = 1e-6)
  expect_equal(out$f_star, rc$nuis$m, tolerance = 1e-6)
})

test_that("single-stage logistic riesz_tmle approximately solves empirical score equation", {
  set.seed(3)
  n <- 500

  df <- data.frame(
    L = rnorm(n)
  )
  df$A <- rbinom(n, 1, plogis(0.4 * df$L))
  linpred <- -0.2 + 0.3 * df$L + 0.5 * df$A
  df$Y <- plogis(linpred)

  rc <- riesz_curve_catalog$cfmean_a(
    a = 1,
    nuis = function(data) {
      mfit <- stats::glm(Y ~ L + A, family = stats::gaussian(), data = data)
      data1 <- data
      data1$A <- 1
      gfit <- stats::glm(I(A == 1) ~ L, family = stats::binomial(), data = data)

      list(
        m = stats::predict(mfit, newdata = data, type = "response"),
        ma = stats::predict(mfit, newdata = data1, type = "response"),
        g = stats::predict(gfit, newdata = data, type = "response")
      )
    }
  )

  rc$fit(df)

  out <- riesz_tmle(
    data = df,
    rc = rc,
    fluctuation_type = "logistic",
    bounds = c(0, 1),
    outcome_col = "Y"
  )

  score_mean <- mean(rc$fit_alpha * (df$Y - out$f_star))
  expect_equal(score_mean, 0, tolerance = 1e-4)
})

test_that("single-stage logistic riesz_tmle is close to truth for bounded counterfactual mean", {
  set.seed(4)
  n <- 8000

  df <- data.frame(
    L = rnorm(n)
  )
  df$A <- rbinom(n, 1, plogis(df$L))

  # bounded outcome in [0, 1]
  df$Y <- plogis(-0.3 + 0.4 * df$L + 0.8 * df$A)

  # truth for E[Y^1] = E[ plogis(-0.3 + 0.4L + 0.8) ]
  truth <- mean(plogis(-0.3 + 0.4 * df$L + 0.8))

  rc <- riesz_curve_catalog$cfmean_a(
    a = 1,
    nuis = function(data) {
      mfit <- stats::glm(Y ~ L + A, family = stats::gaussian(), data = data)
      data1 <- data
      data1$A <- 1
      gfit <- stats::glm(I(A == 1) ~ L, family = stats::binomial(), data = data)

      list(
        m = stats::predict(mfit, newdata = data, type = "response"),
        ma = stats::predict(mfit, newdata = data1, type = "response"),
        g = stats::predict(gfit, newdata = data, type = "response")
      )
    }
  )

  out <- riesz_tmle(
    data = df,
    rc = rc,
    fluctuation_type = "logistic",
    bounds = c(0, 1),
    outcome_col = "Y"
  )

  expect_equal(out$estimate, truth, tolerance = 0.05)
})

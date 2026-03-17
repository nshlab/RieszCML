test_that("RieszCurve initialize errors for unnamed nuisance list", {
  expect_error(
    RieszCurve$new(
      nuis = list(c(1, 2, 3)),
      alpha = ~ A,
      f = ~ Y,
      h = ~ Y
    ),
    "if nuis is a list it must be named."
  )
})

test_that("RieszCurve initialize errors if nuis is neither list nor function", {
  expect_error(
    RieszCurve$new(
      nuis = 1,
      alpha = ~ A,
      f = ~ Y,
      h = ~ Y
    ),
    "must be either a named list or function"
  )
})

test_that("RieszCurve fit computes alpha f h and ic correctly with nuisance list", {
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

  rc$fit(df)

  expect_equal(rc$fit_alpha, c(2, 0, 1.25))
  expect_equal(rc$fit_f, c(1.2, 1.7, 3.5))
  expect_equal(rc$fit_h, c(2, 2, 2))

  expected_ic <- c(2 + 2 * (1 - 1.2),
                   2 + 0 * (2 - 1.7),
                   2 + 1.25 * (4 - 3.5))
  expect_equal(rc$fit_ic, expected_ic)
})

test_that("RieszCurve fit computes same results when nuis is a function", {
  df <- data.frame(
    Y = c(1, 2, 4),
    A = c(1, 0, 1)
  )

  nuis_fun <- function(data) {
    list(
      g = c(0.5, 0.25, 0.8),
      m = c(1.2, 1.7, 3.5),
      ma = c(2.0, 2.0, 2.0)
    )
  }

  rc <- RieszCurve$new(
    nuis = nuis_fun,
    alpha = ~ A / g,
    f = ~ m,
    h = ~ ma,
    ic_expr = ~ h + alpha * (Y - f)
  )

  rc$fit(df)

  expect_equal(rc$fit_alpha, c(2, 0, 1.25))
  expect_equal(rc$fit_f, c(1.2, 1.7, 3.5))
  expect_equal(rc$fit_h, c(2, 2, 2))

  expected_ic <- c(2 + 2 * (1 - 1.2),
                   2 + 0 * (2 - 1.7),
                   2 + 1.25 * (4 - 3.5))
  expect_equal(rc$fit_ic, expected_ic)
})

test_that("RieszCurve accepts numeric alpha f and h", {
  df <- data.frame(
    Y = c(1, 2, 3),
    A = c(1, 0, 1)
  )

  rc <- RieszCurve$new(
    nuis = list(),
    alpha = c(1, 2, 3),
    f = c(0.5, 1.5, 2.5),
    h = c(10, 20, 30),
    ic_expr = ~ h + alpha * (Y - f)
  )

  rc$fit(df)

  expect_equal(rc$fit_alpha, c(1, 2, 3))
  expect_equal(rc$fit_f, c(0.5, 1.5, 2.5))
  expect_equal(rc$fit_h, c(10, 20, 30))
  expect_equal(rc$fit_ic, c(10.5, 21, 31.5))
})

test_that("RieszCurve formulas can see both data columns and nuisance objects", {
  df <- data.frame(
    Y = c(1, 2, 3),
    A = c(1, 0, 1),
    L = c(2, 4, 6)
  )

  rc <- RieszCurve$new(
    nuis = list(g = c(0.5, 0.25, 0.75), m = c(1, 1, 1)),
    alpha = ~ A / g,
    f = ~ m + L,
    h = ~ 2 * m + L,
    ic_expr = ~ h + alpha * (Y - f)
  )

  rc$fit(df)

  expect_equal(rc$fit_alpha, c(2, 0, 4/3))
  expect_equal(rc$fit_f, c(3, 5, 7))
  expect_equal(rc$fit_h, c(4, 6, 8))
})

test_that("RieszCurve h formula can reference f", {
  df <- data.frame(
    Y = c(1, 2, 3),
    A = c(1, 0, 1)
  )

  rc <- RieszCurve$new(
    nuis = list(m = c(1, 2, 3), delta = c(10, 20, 30)),
    alpha = ~ c(1, 1, 1),
    f = ~ m,
    h = ~ f + delta,
    ic_expr = ~ h + alpha * (Y - f)
  )

  rc$fit(df)

  expect_equal(rc$fit_f, c(1, 2, 3))
  expect_equal(rc$fit_h, c(11, 22, 33))
})

test_that("RieszCurve fit_ic is unnamed", {
  df <- data.frame(
    Y = c(1, 2, 3),
    A = c(1, 0, 1)
  )

  rc <- RieszCurve$new(
    nuis = list(
      g = stats::setNames(c(0.5, 0.5, 0.5), c("a", "b", "c")),
      m = stats::setNames(c(1, 2, 3), c("a", "b", "c")),
      ma = stats::setNames(c(2, 2, 2), c("a", "b", "c"))
    ),
    alpha = ~ A / g,
    f = ~ m,
    h = ~ ma,
    ic_expr = ~ h + alpha * (Y - f)
  )

  rc$fit(df)
  expect_null(names(rc$fit_ic))
})

test_that("expit and logit are approximate inverses on the interior of (0, 1)", {
  x <- seq(0.01, 0.99, length.out = 25)
  expect_equal(expit(logit(x)), x, tolerance = 1e-8)
})

test_that("to01 and from01 are approximate inverses", {
  bounds <- c(-2, 3)
  x <- seq(-2, 3, length.out = 11)

  expect_equal(from01(to01(x, bounds), bounds), x, tolerance = 1e-10)
})

test_that("to01 maps endpoints correctly", {
  bounds <- c(-5, 5)
  expect_equal(to01(-5, bounds), 0)
  expect_equal(to01(5, bounds), 1)
  expect_equal(to01(0, bounds), 0.5)
})

test_that("from01 maps endpoints correctly", {
  bounds <- c(-5, 5)
  expect_equal(from01(0, bounds), -5)
  expect_equal(from01(1, bounds), 5)
  expect_equal(from01(0.5, bounds), 0)
})

test_that("clip01 clips outside values and leaves interior values alone", {
  x <- c(-0.1, 0, 0.001, 0.5, 0.999, 1, 1.2)
  out <- clip01(x, clip = 0.01)

  expect_equal(out, c(0.01, 0.01, 0.01, 0.5, 0.99, 0.99, 0.99))
})

test_that("to01 errors for bad bounds", {
  expect_error(to01(1, NULL), "`bounds`")
  expect_error(to01(1, 1), "`bounds`")
  expect_error(to01(1, c(0, 1, 2)), "`bounds`")
  expect_error(to01(1, c(1, 1)), "must be >")
  expect_error(to01(1, c(2, 1)), "must be >")
})

test_that("from01 errors for bad bounds", {
  expect_error(from01(0.5, NULL), "`bounds`")
  expect_error(from01(0.5, 1), "`bounds`")
  expect_error(from01(0.5, c(0, 1, 2)), "`bounds`")
  expect_error(from01(0.5, c(1, 1)), "must be >")
  expect_error(from01(0.5, c(2, 1)), "must be >")
})

test_that("clip01 errors for invalid clip values", {
  expect_error(clip01(c(0.1, 0.2), clip = 0), "`clip`")
  expect_error(clip01(c(0.1, 0.2), clip = -0.1), "`clip`")
  expect_error(clip01(c(0.1, 0.2), clip = 0.5), "`clip`")
  expect_error(clip01(c(0.1, 0.2), clip = 0.6), "`clip`")
  expect_error(clip01(c(0.1, 0.2), clip = NA), "`clip`")
})

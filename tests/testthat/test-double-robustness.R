# Tests for the correctness of the composed influence curve, the sequential
# TMLE, and double robustness. These address rOpenSci standards G5.4-G5.6
# (correctness and parameter-recovery tests, fixed seeds) and G5.2 (error
# behaviour).
#
# The DR stress tests are the decisive detectors for two historical bugs:
#   (1) prefix (rather than suffix) cumulative products of Riesz
#       representers in ComposedRieszCurve / riesz_tmle, which preserved
#       consistency under correct specification but destroyed double
#       robustness and invalidated the variance; and
#   (2) updating counterfactual predictions with the observed-data clever
#       covariate during TMLE, so off-support / unobserved units never
#       received the fluctuation.

# ---- helpers ----------------------------------------------------------------

sim_two_timepoint <- function(n) {
  # E[Y^{1,1}] = 1.0 + 1.5 * (0.5 + 2.0) + 3.0 = 7.75
  L1 <- rnorm(n)
  A1 <- rbinom(n, 1, plogis(-0.2 + 1.0 * L1))
  L2 <- 0.5 + 0.7 * L1 + 2.0 * A1 + rnorm(n)
  A2 <- rbinom(n, 1, plogis(-0.3 + 0.6 * A1 + 1.2 * L2))
  Y <- 1.0 + 1.5 * L2 + 3.0 * A2 + rnorm(n)
  data.frame(L1 = L1, A1 = A1, L2 = L2, A2 = A2, Y = Y)
}
TRUTH_MU11 <- 7.75

make_stage1 <- function(misspec_q = FALSE) {
  RieszCurve$new(
    nuis = function(data) {
      m_fit <- if (misspec_q) lm(Y ~ 1, data = data) else
        lm(Y ~ L1 + A1 + L2 + A2, data = data)
      d1 <- data; d1$A2 <- 1
      g_fit <- glm(A2 ~ L1 + A1 + L2, family = binomial(), data = data)
      list(
        m  = predict(m_fit, newdata = data),
        m1 = predict(m_fit, newdata = d1),
        g  = predict(g_fit, newdata = data, type = "response")
      )
    },
    h = ~ m1, f = ~ m,
    alpha = ~ A2 / g, alpha_star = ~ 1 / g,
    ic_expr = ~ h + alpha * (Y - f)
  )
}

make_stage2 <- function(misspec_q = FALSE) {
  RieszCurve$new(
    nuis = function(data) {
      m_fit <- if (misspec_q) lm(`h[j=1]` ~ 1, data = data) else
        lm(`h[j=1]` ~ L1 + A1, data = data)
      d1 <- data; d1$A1 <- 1
      g_fit <- glm(A1 ~ L1, family = binomial(), data = data)
      list(
        m  = predict(m_fit, newdata = data),
        m1 = predict(m_fit, newdata = d1),
        g  = predict(g_fit, newdata = data, type = "response")
      )
    },
    h = ~ m1, f = ~ m,
    alpha = ~ A1 / g, alpha_star = ~ 1 / g,
    ic_expr = ~ h + alpha * (`h[j=1]` - f)
  )
}

make_composed <- function(misspec_q = FALSE) {
  ComposedRieszCurve$new(rc_list = list(
    make_stage1(misspec_q), make_stage2(misspec_q)
  ))
}

# ---- 1. exact structure of the composed influence curve ---------------------

test_that("composed IC uses suffix cumulative products of representers", {
  # static nuisances so the IC is checkable by hand
  n <- 6
  set.seed(42)
  d <- data.frame(Y = rnorm(n))
  a1 <- c(2, 1, 0.5, 1, 3, 0.25)   # innermost stage representer
  f1 <- rnorm(n); h1 <- rnorm(n)
  a2 <- c(3, 2, 1, 1, 0.5, 2)      # outermost stage representer
  f2 <- rnorm(n); h2 <- rnorm(n)

  rc1 <- RieszCurve$new(nuis = list(a = a1, fv = f1, hv = h1),
                        alpha = ~ a, f = ~ fv, h = ~ hv,
                        ic_expr = ~ h + alpha * (Y - f))
  rc2 <- RieszCurve$new(nuis = list(a = a2, fv = f2, hv = h2),
                        alpha = ~ a, f = ~ fv, h = ~ hv,
                        ic_expr = ~ h + alpha * (`h[j=1]` - f))
  rc <- ComposedRieszCurve$new(rc_list = list(rc1, rc2))
  rc$fit(d)

  # Theorem 2 with stage 1 = innermost: the residual against the innermost
  # regression carries the FULL product a1 * a2; the outermost carries a2.
  ic_by_hand <- h2 + a2 * a1 * (d$Y - f1) + a2 * (h1 - f2)
  expect_equal(rc$fit_ic, ic_by_hand, tolerance = 1e-12)
  expect_equal(rc$fit_omega[[1]], a1 * a2, tolerance = 1e-12)
  expect_equal(rc$fit_omega[[2]], a2, tolerance = 1e-12)
})

# ---- 2. one-step: consistency and double robustness -------------------------

test_that("composed one-step recovers truth with correct nuisances", {
  set.seed(101)
  df <- sim_two_timepoint(20000)
  fit <- riesz_estimate(df, make_composed(misspec_q = FALSE))
  expect_lt(abs(fit$estimate - TRUTH_MU11), 4 * fit$se)
})

test_that("composed one-step is doubly robust (Q misspecified, alpha correct)", {
  # This is the test that detects reversed representer products: with prefix
  # products the estimate is biased by roughly -1.6 in this DGP.
  set.seed(102)
  df <- sim_two_timepoint(20000)
  fit <- riesz_estimate(df, make_composed(misspec_q = TRUE))
  expect_lt(abs(fit$estimate - TRUTH_MU11), 4 * fit$se)
})

test_that("composed one-step DR holds on average over replications", {
  skip_on_cran()
  set.seed(103)
  ests <- replicate(60, {
    df <- sim_two_timepoint(1500)
    riesz_estimate(df, make_composed(misspec_q = TRUE))$estimate
  })
  expect_lt(abs(mean(ests) - TRUTH_MU11), 0.08)
})

# ---- 3. sequential TMLE ------------------------------------------------------

test_that("composed TMLE recovers truth and solves its score equation", {
  set.seed(104)
  df <- sim_two_timepoint(20000)
  fit <- expect_no_warning(riesz_tmle(
    data = df, rc = make_composed(FALSE),
    fluctuation_type = "logistic",
    bounds = range(df$Y) + c(-1, 1),
    outcome_col = "Y"
  ))
  expect_lt(abs(fit$estimate - TRUTH_MU11), 4 * fit$se)
  # targeted EIF estimating equation solved
  expect_lt(abs(mean(fit$ic_star) - fit$estimate), 0.01 * fit$se)
})

test_that("composed TMLE is doubly robust (Q misspecified, alpha correct)", {
  # Detects both the reversed products AND the observed-covariate h* update:
  # under either bug the fluctuation cannot repair intercept-only regressions.
  set.seed(105)
  df <- sim_two_timepoint(20000)
  fit <- riesz_tmle(
    data = df, rc = make_composed(TRUE),
    fluctuation_type = "logistic",
    bounds = range(df$Y) + c(-1, 1),
    outcome_col = "Y"
  )
  expect_lt(abs(fit$estimate - TRUTH_MU11), 5 * fit$se)
  expect_lt(abs(mean(fit$ic_star) - fit$estimate), 0.01 * fit$se)
})

test_that("composed TMLE warns when a stage lacks alpha_star", {
  set.seed(106)
  df <- sim_two_timepoint(2000)
  s1 <- make_stage1(FALSE); s1$alpha_star <- NULL
  rc <- ComposedRieszCurve$new(rc_list = list(s1, make_stage2(FALSE)))
  expect_warning(
    riesz_tmle(df, rc, fluctuation_type = "logistic",
               bounds = range(df$Y) + c(-1, 1), outcome_col = "Y"),
    regexp = "alpha_star"
  )
})

# ---- 4. single-curve TMLE: counterfactual mean, binary outcome --------------

test_that("cfmean TMLE recovers truth; alpha_star and set targeting agree", {
  set.seed(107)
  n <- 20000
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.5 * L))
  Y <- rbinom(n, 1, plogis(-0.3 + 0.7 * A + 0.4 * L))
  df <- data.frame(L = L, A = A, Y = Y)
  truth <- integrate(function(l) plogis(-0.3 + 0.7 + 0.4 * l) * dnorm(l),
                     -Inf, Inf)$value

  nuis_fun <- function(data) {
    m_fit <- glm(Y ~ L + A, family = binomial(), data = data)
    g_fit <- glm(I(A == 1) ~ L, family = binomial(), data = data)
    d1 <- data; d1$A <- 1
    list(
      m  = predict(m_fit, newdata = data, type = "response"),
      ma = predict(m_fit, newdata = d1, type = "response"),
      g  = predict(g_fit, newdata = data, type = "response")
    )
  }

  fit_astar <- riesz_tmle(df, riesz_curve_catalog$cfmean_a(a = 1, nuis = nuis_fun),
                          fluctuation_type = "logistic", bounds = c(0, 1),
                          outcome_col = "Y")

  rc_set <- riesz_curve_catalog$cfmean_a(a = 1, nuis = nuis_fun)
  rc_set$targeting_steps <- list(m = list(), ma = list(set = list(A = 1)))
  fit_set <- riesz_tmle(df, rc_set, fluctuation_type = "logistic",
                        bounds = c(0, 1), outcome_col = "Y")

  expect_lt(abs(fit_astar$estimate - truth), 4 * fit_astar$se)
  expect_equal(fit_astar$estimate, fit_set$estimate, tolerance = 1e-10)
  expect_lt(abs(mean(fit_astar$ic_star) - fit_astar$estimate),
            0.01 * fit_astar$se)
})

# ---- 5. missingness: plug-in must be updated with 1/pi ----------------------

test_that("missing_mean TMLE de-biases a misspecified outcome regression", {
  # With m intercept-only and pi correct, the TMLE plug-in is consistent only
  # if unobserved units receive the fluctuation through 1/pi (not R/pi).
  set.seed(108)
  n <- 30000
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(L))
  Y_full <- 2 + L + 1.5 * A + rnorm(n, sd = 0.25)
  R <- rbinom(n, 1, plogis(-0.5 + 0.5 * L + 0.3 * A))
  Y <- ifelse(R == 1, Y_full, 0)
  df <- data.frame(L = L, A = A, R = R, Y = Y)
  truth <- 2 + 1.5 * 0.5  # E[L] = 0, E[A] = 0.5 by symmetry

  rc <- riesz_curve_catalog$missing_mean(nuis = function(data) {
    obs <- data[data$R == 1, , drop = FALSE]
    pi_fit <- glm(R ~ L + A, family = binomial(), data = data)
    list(
      m  = rep(mean(obs$Y), nrow(data)),  # deliberately misspecified
      pi = predict(pi_fit, newdata = data, type = "response")
    )
  })

  # the naive plug-in is visibly biased in this DGP ...
  rc$fit(df)
  expect_gt(abs(mean(rc$fit_h) - truth), 0.05)

  # ... but both the one-step and the TMLE are consistent
  fit_os <- riesz_estimate(df, rc)
  expect_lt(abs(fit_os$estimate - truth), 4 * fit_os$se)

  fit_tmle <- riesz_tmle(df, rc, fluctuation_type = "logistic",
                         bounds = range(df$Y) + c(-1, 1), outcome_col = "Y")
  expect_lt(abs(fit_tmle$estimate - truth), 5 * fit_tmle$se)
  expect_lt(abs(mean(fit_tmle$ic_star) - fit_tmle$estimate),
            0.01 * fit_tmle$se)
})

# ---- 6. ATT via ratio delta method ------------------------------------------

test_that("riesz_att recovers a constant additive effect", {
  set.seed(109)
  n <- 20000
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(L))
  Y <- 2 + L + 1.5 * A + rnorm(n, sd = 0.25)
  df <- data.frame(L = L, A = A, Y = Y)

  fit <- riesz_att(df, nuis = function(data) {
    g_fit  <- glm(A ~ L, family = binomial(), data = data)
    m0_fit <- lm(Y ~ L, data = subset(data, A == 0))
    list(
      m0 = predict(m0_fit, newdata = data),
      g  = predict(g_fit, newdata = data, type = "response")
    )
  })
  expect_lt(abs(fit$estimate - 1.5), 4 * fit$se)
})

test_that("the removed catalog att entry errors informatively", {
  expect_error(riesz_curve_catalog$att(nuis = identity), regexp = "riesz_att")
})

# ---- 7. input validation ------------------------------------------------------

test_that("malformed targeting_steps are rejected at construction", {
  nuis <- list(m = 1, m1 = 1, g = 0.5)
  # the historical silent failure mode: list(A2 = 1) without `set = `
  expect_error(
    RieszCurve$new(nuis = nuis, alpha = ~ A2/g, f = ~ m, h = ~ m1,
                   ic_expr = ~ h + alpha * (Y - f),
                   targeting_steps = list(m = list(), m1 = list(A2 = 1))),
    regexp = "set = list"
  )
  expect_error(
    RieszCurve$new(nuis = nuis, alpha = ~ A2/g, f = ~ m, h = ~ m1,
                   targeting_steps = list(m1 = "not a list")),
    regexp = "must be a list"
  )
})

test_that("nuisance functions that drop rows are rejected", {
  rc <- RieszCurve$new(
    nuis = function(data) {
      d <- data[data$A == 1, , drop = FALSE]  # silently drops rows
      list(m = d$Y, ma = d$Y, g = rep(0.5, nrow(d)))
    },
    alpha = ~ I(A == 1)/g, f = ~ m, h = ~ ma,
    ic_expr = ~ h + alpha * (Y - f)
  )
  df <- data.frame(A = c(0, 1, 1, 0), Y = 1:4)
  expect_error(rc$fit(df), regexp = "drop")
})

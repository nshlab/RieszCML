# The outcome (and treatment) column need not be named `Y` / `A`.
# Every entry point should give the same answer after renaming columns,
# and anything that cannot work should fail loudly rather than silently.

sim <- function(n = 2000, seed = 1) {
  set.seed(seed)
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(L))
  Y <- rbinom(n, 1, plogis(-0.5 + A + 0.8 * L))
  data.frame(L = L, A = A, Y = Y)
}

# nuisance factory that reads outcome/treatment names from arguments
make_nuis <- function(y, a) {
  function(data) {
    fm <- stats::reformulate(c("L", a), response = as.name(y))
    fg <- stats::reformulate("L", response = a)
    m_fit <- glm(fm, family = binomial(), data = data)
    g_fit <- glm(fg, family = binomial(), data = data)
    d1 <- data; d1[[a]] <- 1
    d0 <- data; d0[[a]] <- 0
    list(
      m  = predict(m_fit, newdata = data, type = "response"),
      ma = predict(m_fit, newdata = d1,   type = "response"),
      m1 = predict(m_fit, newdata = d1,   type = "response"),
      m0 = predict(m_fit, newdata = d0,   type = "response"),
      g  = predict(g_fit, newdata = data, type = "response")
    )
  }
}

rename <- function(df, y, a) {
  names(df)[names(df) == "Y"] <- y
  names(df)[names(df) == "A"] <- a
  df
}

# includes a non-syntactic name to exercise backtick quoting
outcome_names <- c("Z", "hiv_infected", "my outcome")

test_that("catalog cfmean_a is invariant to outcome/treatment names", {
  df <- sim()
  ref_os   <- riesz_estimate(df, riesz_curve_catalog$cfmean_a(1, make_nuis("Y", "A")))
  ref_tmle <- riesz_tmle(df, riesz_curve_catalog$cfmean_a(1, make_nuis("Y", "A")),
                         bounds = c(0, 1), outcome_col = "Y")
  for (y in outcome_names) {
    d <- rename(df, y, "trt")
    rc <- riesz_curve_catalog$cfmean_a(1, make_nuis(y, "trt"),
                                       outcome = y, treatment = "trt")
    expect_equal(riesz_estimate(d, rc)$estimate, ref_os$estimate)
    expect_equal(riesz_estimate(d, rc)$se, ref_os$se)
    fit <- riesz_tmle(d, rc, bounds = c(0, 1), outcome_col = y)
    expect_equal(fit$estimate, ref_tmle$estimate)
  }
})

test_that("catalog ate is invariant to outcome/treatment names", {
  df <- sim()
  ref <- riesz_tmle(df, riesz_curve_catalog$ate(make_nuis("Y", "A")),
                    bounds = c(0, 1), outcome_col = "Y")$estimate
  d <- rename(df, "Z", "trt")
  rc <- riesz_curve_catalog$ate(make_nuis("Z", "trt"), outcome = "Z", treatment = "trt")
  expect_equal(riesz_tmle(d, rc, bounds = c(0, 1), outcome_col = "Z")$estimate, ref)
})

test_that("riesz_att is invariant to outcome/treatment names", {
  df <- sim()
  nuis_att <- function(y, a) function(data) {
    g_fit  <- glm(stats::reformulate("L", a), family = binomial(), data = data)
    m0_fit <- glm(stats::reformulate("L", y), family = binomial(),
                  data = data[data[[a]] == 0, ])
    list(m0 = predict(m0_fit, newdata = data, type = "response"),
         g  = predict(g_fit,  newdata = data, type = "response"))
  }
  ref <- riesz_att(df, nuis_att("Y", "A"))
  d <- rename(df, "Z", "trt")
  got <- riesz_att(d, nuis_att("Z", "trt"), outcome = "Z", treatment = "trt")
  expect_equal(got$estimate, ref$estimate)
  expect_equal(got$se, ref$se)
})

test_that("missing_mean and subgroup_mean accept other column names", {
  df <- sim(); set.seed(2)
  df$R <- rbinom(nrow(df), 1, plogis(0.5 + 0.5 * df$L))
  df$Y[df$R == 0] <- 0
  nuis_mm <- function(y, r) function(data) {
    obs <- data[data[[r]] == 1, ]
    m_fit <- glm(stats::reformulate(c("L", "A"), y), binomial(), data = obs)
    p_fit <- glm(stats::reformulate(c("L", "A"), r), binomial(), data = data)
    list(m = predict(m_fit, newdata = data, type = "response"),
         pi = predict(p_fit, newdata = data, type = "response"))
  }
  ref <- riesz_estimate(df, riesz_curve_catalog$missing_mean(nuis_mm("Y", "R")))
  d <- df; names(d)[names(d) == "Y"] <- "Z"; names(d)[names(d) == "R"] <- "seen"
  got <- riesz_estimate(d, riesz_curve_catalog$missing_mean(
    nuis_mm("Z", "seen"), outcome = "Z", observed = "seen"))
  expect_equal(got$estimate, ref$estimate)

  df$V <- as.integer(df$L > 0)
  nuis_sg <- function(y, v) function(data) {
    m_fit <- lm(stats::reformulate(v, y), data = data)
    d1 <- data; d1[[v]] <- 1
    list(m = predict(m_fit, data), mv = predict(m_fit, d1), pv = mean(data[[v]] == 1))
  }
  ref <- riesz_estimate(df, riesz_curve_catalog$subgroup_mean(1, nuis_sg("Y", "V")))
  d <- df; names(d)[names(d) == "Y"] <- "Z"; names(d)[names(d) == "V"] <- "grp"
  got <- riesz_estimate(d, riesz_curve_catalog$subgroup_mean(
    1, nuis_sg("Z", "grp"), outcome = "Z", subgroup = "grp"))
  expect_equal(got$estimate, ref$estimate)
})

test_that("ComposedRieszCurve works with a non-Y outcome when told its name", {
  df <- rename(sim(), "Z", "A")
  stages <- function() list(
    riesz_curve_catalog$cfmean_a(1, make_nuis("Z", "A"), outcome = "Z"),
    RieszCurve$new(nuis = list(), alpha = 1, alpha_star = 1,
                   f = ~ `h[j=1]`, h = ~ `h[j=1]`, ic_expr = ~ h)
  )
  single <- riesz_estimate(df, riesz_curve_catalog$cfmean_a(1, make_nuis("Z", "A"), outcome = "Z"))
  rc <- ComposedRieszCurve$new(stages(), outcome_col = "Z")
  expect_equal(riesz_estimate(df, rc)$estimate, single$estimate)
  # riesz_tmle picks up the composed curve's outcome_col automatically ...
  expect_silent(fit <- riesz_tmle(df, ComposedRieszCurve$new(stages(), outcome_col = "Z"),
                                  bounds = c(0, 1)))
  expect_true(is.finite(fit$estimate))
  # ... and refuses a conflicting one
  expect_error(riesz_tmle(df, ComposedRieszCurve$new(stages(), outcome_col = "Z"),
                          bounds = c(0, 1), outcome_col = "L"), "does not match")
})

test_that("ComposedRieszCurve with default outcome_col errors instead of using a stray `Y`", {
  df <- rename(sim(), "Z", "A")
  df$Y <- rnorm(nrow(df))  # unrelated column that happens to be called Y
  stages <- list(
    riesz_curve_catalog$cfmean_a(1, make_nuis("Z", "A"), outcome = "Z"),
    RieszCurve$new(nuis = list(), alpha = 1, alpha_star = 1,
                   f = ~ `h[j=1]`, h = ~ `h[j=1]`, ic_expr = ~ h)
  )
  expect_error(ComposedRieszCurve$new(stages), "does not reference")
})

test_that("outcome columns colliding with nuisance/reserved names fail loudly", {
  for (bad in c("h", "f", "alpha", "m", "g")) {
    d <- rename(sim(), bad, "A")
    rc <- riesz_curve_catalog$cfmean_a(1, make_nuis(bad, "A"), outcome = bad)
    expect_error(riesz_estimate(d, rc), "share a name")
  }
})

test_that("catalog defaults still assume Y and A (backwards compatibility)", {
  df <- sim()
  rc <- riesz_curve_catalog$cfmean_a(1, make_nuis("Y", "A"))
  expect_true("Y" %in% all.vars(rc$ic_expr))
  expect_true("A" %in% all.vars(rc$alpha))
})

#' A Catalog of Pre-Specified Riesz Curves for Common Estimands
#'
#' `riesz_curve_catalog` is a list of factory functions, each returning a
#' [RieszCurve] for a commonly used statistical estimand. Every factory takes
#' a `nuis` argument: a function of `data` returning the named list of
#' nuisance vectors documented for that entry.
#'
#' @export
#' @examples
#' # Counterfactual means E[Y^a], a = 0, 1
#' rc_cfmean_1 <- riesz_curve_catalog$cfmean_a(
#'   a = 1,
#'   nuis = function(data) {
#'     m_fit <- lm(Y ~ L + A, data = data)
#'     g_fit <- glm(I(A == 1) ~ L, family = binomial(), data = data)
#'
#'     dataAeq_a <- data
#'     dataAeq_a$A <- 1
#'
#'     list(
#'       m = predict(m_fit, newdata = data),
#'       ma = predict(m_fit, newdata = dataAeq_a),
#'       g = predict(g_fit, newdata = data, type = "response"))
#'   })
#'
#' rc_cfmean_0 <- riesz_curve_catalog$cfmean_a(
#'   a = 0,
#'   nuis = function(data) {
#'     m_fit <- lm(Y ~ L + A, data = data)
#'     g_fit <- glm(I(A == 0) ~ L, family = binomial(), data = data)
#'
#'     dataAeq_a <- data
#'     dataAeq_a$A <- 0
#'
#'     list(
#'       m = predict(m_fit, newdata = data),
#'       ma = predict(m_fit, newdata = dataAeq_a),
#'       g = predict(g_fit, newdata = data, type = "response"))
#'   })
#'
#' set.seed(1)
#' n <- 5000
#' df <- data.frame(L = rnorm(n))
#' df$A <- rbinom(n, size = 1, prob = plogis(df$L))
#' df$Y <- df$L + df$A + rnorm(n, sd = .1)
#'
#' riesz_estimate(data = df, rc = rc_cfmean_0)
#' riesz_estimate(data = df, rc = rc_cfmean_1)
riesz_curve_catalog <- list()

#' @section Counterfactual mean `cfmean_a(a, nuis)`:
#' Target: \eqn{E[E(Y | A = a, L)]}. The `nuis` function must return
#' `m` (\eqn{E[Y | A, L]} at observed data), `ma` (\eqn{E[Y | A = a, L]}),
#' and `g` (\eqn{P(A = a | L)}). Data columns `A` and `Y` are assumed to be
#' named literally `A` and `Y`.
riesz_curve_catalog$cfmean_a <- function(
    a,
    nuis) {

  alpha_formula <- as.formula(paste0('~ I(A == ', a, ')/g'))

  return(RieszCurve$new(
    nuis = nuis,
    alpha = alpha_formula,
    alpha_star = ~ 1/g,
    f = ~ m,
    h = ~ ma,
    ic_expr = ~ h + alpha * (Y - f),
    targeting_steps = list(
      m = list(),
      ma = list(alpha_star = ~ 1/g)
    )
  ))
}


riesz_curve_catalog$cfmean_a1 <- function(nuis) {
  riesz_curve_catalog$cfmean_a(a = 1, nuis = nuis)
}

riesz_curve_catalog$cfmean_a0 <- function(nuis) {
  riesz_curve_catalog$cfmean_a(a = 0, nuis = nuis)
}

#' @section Average treatment effect `ate(nuis)`:
#' Target: \eqn{E[E(Y | A = 1, L) - E(Y | A = 0, L)]}. The `nuis` function
#' must return `m`, `m1`, `m0`, and `g` (\eqn{P(A = 1 | L)}). Because `h`
#' is a contrast of two counterfactual evaluations, the curve carries no
#' single `alpha_star` and is not intended for use as a stage of a
#' [ComposedRieszCurve]; compose the two `cfmean_a` curves instead and
#' contrast with [riesz_delta_difference()].
riesz_curve_catalog$ate <- function(nuis) {
  RieszCurve$new(
    nuis = nuis,
    alpha = ~ A / g - (1 - A) / (1 - g),
    f = ~ m,
    h = ~ m1 - m0,
    ic_expr = ~ h + alpha * (Y - f),
    targeting_steps = list(
      m  = list(),
      m1 = list(alpha_star = ~ 1 / g),
      m0 = list(alpha_star = ~ -1 / (1 - g))
    )
  )
}

#' @section Average treatment effect on the treated:
#' The ATT \eqn{\psi = E[Y^1 - Y^0 | A = 1] = E[A(Y - m_0(L))] / P(A = 1)}
#' is a *ratio* of two linear functionals, so it is estimated by combining
#' two [RieszCurve] fits with [riesz_delta_ratio()]; see [riesz_att()]. The
#' numerator curve is `att_theta(nuis)`, targeting
#' \eqn{\theta = E[A(Y - m_0(L))]} with the doubly robust uncentered
#' influence curve
#' \deqn{A(Y - m_0) - (1 - A)\frac{g}{1 - g}(Y - m_0).}
#' The `nuis` function must return `m0` (\eqn{E[Y | A = 0, L]} predicted for
#' all units) and `g` (\eqn{P(A = 1 | L)}). TMLE is not yet supported for
#' the ATT; use the one-step estimator via [riesz_att()].
riesz_curve_catalog$att_theta <- function(nuis) {
  RieszCurve$new(
    nuis = nuis,
    alpha = ~ -(1 - A) * g / (1 - g),
    f = ~ m0,
    h = ~ A * (Y - m0),
    ic_expr = ~ h + alpha * (Y - f)
  )
}

riesz_curve_catalog$treated_proportion <- function() {
  RieszCurve$new(
    nuis = list(),
    alpha = 0,
    f = 0,
    h = ~ A,
    ic_expr = ~ A
  )
}

# The previous `att` entry was statistically invalid (its influence-curve
# expression reduced to the non-doubly-robust A * (Y - m0) / mean(A), its
# `alpha` was not the ATT Riesz representer, and its documented `nuis` did
# not supply the components referenced by `f` and `h`). It is retained only
# as an informative error.
riesz_curve_catalog$att <- function(...) {
  stop("`riesz_curve_catalog$att` has been removed: the previous entry was ",
       "not doubly robust and its variance was incorrect. Use `riesz_att()`",
       ", which estimates the ATT as a ratio of two linear functionals via ",
       "the delta method.")
}

#' Estimate the Average Treatment Effect on the Treated
#'
#' One-step, doubly robust estimation of
#' \eqn{\psi = E[Y^1 - Y^0 | A = 1] = E[A(Y - m_0(L))] / P(A = 1)}
#' by combining the numerator functional \eqn{E[A(Y - m_0)]} (see
#' `riesz_curve_catalog$att_theta`) and \eqn{P(A = 1)} via
#' [riesz_delta_ratio()], which accounts for the estimation of the
#' denominator in the variance.
#'
#' @param data A `data.frame`-like object with columns `A` and `Y`.
#' @param nuis A function of `data` returning a named list with elements
#'   `m0` (predictions of \eqn{E[Y | A = 0, L]} for all units) and `g`
#'   (predictions of \eqn{P(A = 1 | L)}).
#' @param significance_alpha Significance level for Wald intervals.
#'
#' @export
#' @examples
#' set.seed(1)
#' n <- 5000
#' L <- rnorm(n)
#' A <- rbinom(n, 1, plogis(L))
#' Y <- 2 + L + 1.5 * A + rnorm(n, sd = 0.25)
#' df <- data.frame(L = L, A = A, Y = Y)
#'
#' riesz_att(df, nuis = function(data) {
#'   g_fit  <- glm(A ~ L, family = binomial(), data = data)
#'   m0_fit <- lm(Y ~ L, data = subset(data, A == 0))
#'   list(
#'     m0 = predict(m0_fit, newdata = data),
#'     g  = predict(g_fit, newdata = data, type = "response")
#'   )
#' })
riesz_att <- function(data, nuis, significance_alpha = 0.05) {
  fit_theta <- riesz_estimate(
    data, riesz_curve_catalog$att_theta(nuis),
    significance_alpha = significance_alpha
  )
  fit_p <- riesz_estimate(
    data, riesz_curve_catalog$treated_proportion(),
    significance_alpha = significance_alpha
  )
  riesz_delta_ratio(
    fit_theta, fit_p,
    significance_alpha = significance_alpha,
    parameter = "ATT = E[A(Y - m0)] / P(A = 1)"
  )
}

#' @section Subgroup mean `subgroup_mean(v, nuis)`:
#' Target: \eqn{E[Y | V = v]}. The `nuis` function must return `m`
#' (\eqn{E[Y | V]} at observed data), `mv` (\eqn{E[Y | V = v]}), and `pv`
#' (a scalar, \eqn{P(V = v)}).
riesz_curve_catalog$subgroup_mean <- function(v, nuis) {

  alpha_formula <- as.formula(
    paste0("~ I(V == ", v, ") / pv")
  )

  RieszCurve$new(
    nuis = nuis,
    alpha = alpha_formula,
    alpha_star = ~ 1 / pv,
    f = ~ m,
    h = ~ mv,
    ic_expr = ~ h + alpha * (Y - f),
    targeting_steps = list(
      m = list(),
      mv = list(alpha_star = ~ 1 / pv)
    )
  )
}

#' @section Mean with outcome missing at random `missing_mean(nuis)`:
#' Target: \eqn{E[Y] = E[E(Y | R = 1, V)]} under MAR given `V`. The `nuis`
#' function must return `m` (predictions of \eqn{E[Y | R = 1, V]} for
#' *every* unit, fit among `R == 1`) and `pi` (\eqn{P(R = 1 | V)}).
#' Note the TMLE plug-in updates `m` with the clever covariate
#' \eqn{1/\pi(V)} (all units), while the fluctuation score uses
#' \eqn{R/\pi(V)}; unobserved units must receive the update too.
riesz_curve_catalog$missing_mean <- function(nuis) {

  RieszCurve$new(
    nuis = nuis,
    alpha = ~ R / pi,
    alpha_star = ~ 1 / pi,
    f = ~ m,
    h = ~ m,
    ic_expr = ~ h + alpha * (Y - f),
    targeting_steps = list(
      m = list(alpha_star = ~ 1 / pi)
    )
  )
}

#' @section Modified treatment policy `mtp(nuis, alpha_star)`:
#' Target: \eqn{E[E(Y | A^d, L)]} for a modified treatment policy
#' \eqn{A^d = d(A, L)}. The `nuis` function must return `m` (natural-value
#' predictions), `m_shift` (predictions at \eqn{A^d}), and the density-ratio
#' components `g_shift` and `g` such that `g_shift / g` is the Riesz
#' representer at observed data. For TMLE (and for use as a stage in a
#' composed TMLE), `alpha_star` must additionally give the representer
#' *evaluated at the shifted exposure* (for an additive shift with density
#' ratio estimated via the classifier trick of Diaz et al. (2021), this is
#' the ratio evaluated at \eqn{A + \delta}); with `alpha_star = NULL`
#' (default) only the one-step estimator [riesz_estimate()] is supported and
#' [riesz_tmle()] will error.
riesz_curve_catalog$mtp <- function(nuis, alpha_star = NULL) {

  targeting <- NULL
  if (!is.null(alpha_star)) {
    targeting <- list(
      m = list(),
      m_shift = list(alpha_star = alpha_star)
    )
  }

  RieszCurve$new(
    nuis = nuis,
    alpha = ~ g_shift / g,
    alpha_star = alpha_star,
    f = ~ m,
    h = ~ m_shift,
    ic_expr = ~ h + alpha * (Y - f),
    targeting_steps = targeting
  )
}

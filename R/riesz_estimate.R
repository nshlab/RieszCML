
#' Estimate a Statistical Parameter using a RieszCurve
#'
#' @export
#' @examples
#'
#' rc_ate <- RieszCurve$new(
#'   nuis = function(data) {
#'
#'     m <- nadir::lnr_glm(data,
#'                         formula = Y ~ L + A)
#'     g <- nadir::lnr_logistic(
#'       data, formula = A ~ L)
#'
#'     datacf0 <- datacf1 <- data
#'     datacf0$A <- 0; datacf1$A <- 1
#'
#'     return(list(
#'       m = m(data),
#'       m1 = m(datacf1),
#'       m0 = m(datacf0),
#'       g = g(data)))
#'   },
#'   alpha = ~ A/g + (1-A)/(1-g),
#'   f = ~ m,
#'   h = ~ m1 - m0,
#'   ic_expr = ~ h + alpha * (Y - f)
#' )
#'
#' df <- tibble::tibble(
#'   L = rnorm(n = 50),
#'   A = rbinom(
#'     n = 50,
#'     size = 1,
#'     prob = plogis(L)),
#'   Y = L + rnorm(n = 50, mean = 5, sd = 1) * A)
#'
#' riesz_estimate(data = df, rc = rc_ate)
#'
#'
riesz_estimate <- function(data, rc, significance_alpha = 0.05) {
  if (!inherits(rc, "RieszCurve") && !inherits(rc, "ComposedRieszCurve")) {
    stop("`rc` must be a RieszCurve or ComposedRieszCurve.")
  }

  rc$fit(data)
  phi <- rc$fit_ic

  if (is.null(phi)) {
    stop("After `rc$fit(data)`, `rc$fit_ic` is NULL.")
  }
  if (!is.numeric(phi) || length(phi) != nrow(data)) {
    stop("`rc$fit_ic` must be a numeric vector of length nrow(data).")
  }

  estimate <- mean(phi)
  var_estimate <- stats::var(phi) / nrow(data)
  se <- sqrt(var_estimate)
  ci_low <- estimate + se * stats::qnorm(significance_alpha / 2)
  ci_high <- estimate + se * stats::qnorm(1 - significance_alpha / 2)

  RieszFit$new(
    estimate = estimate,
    var = var_estimate,
    se = se,
    ci_low = ci_low,
    ci_high = ci_high,
    significance_alpha = significance_alpha,
    estimator = "one-step",
    n = nrow(data),
    ic = phi
  )
}

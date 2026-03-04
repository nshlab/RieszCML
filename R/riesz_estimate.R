
#' Estimate a Statistical Parameter using a RieszRepresenter
#'
#'@examples
#'
#' rr_ate <- RieszRepresenter$new(
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
#' riesz_estimate(data = df, rr = rr_ate)
#'
#'
riesz_estimate <- function(data, rr) {
  if (! inherits(rr, "RieszRepresenter") & ! inherits(rr, "ComposedRieszRepresenter")) {
    stop("rr must be a RieszRepresenter or ComposedRieszRepresenter.")
  }

  if (inherits(rr, 'RieszRepresenter') || inherits(rr, 'ComposedRieszRepresenter')) {
    rr$fit(data)
    phi <- rr$fit_ic

    estimate <- mean(phi)
    var_estimate <- var(phi)/nrow(data)
    se <- sqrt(var_estimate)
    ci_low <- estimate + se * qnorm(0.025)
    ci_high <- estimate + se * qnorm(0.975)

    return(list(
      estimate = estimate, var = var_estimate, se = se, ci_low = ci_low, ci_high = ci_high
    ))
  }
}

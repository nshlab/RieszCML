
#' @export
#' @examples
#' rc_cfmean_1 <- riesz_curve_catalog$cfmean_a(
#'   L = 'L',
#'   A = 'A',
#'   a = 1,
#'   Y = 'Y',
#'
#'   nuis = function(data) {
#'     m <- nadir::lnr_glm(data, Y ~ L + A)
#'     dataAeq_a <- data
#'     dataAeq_a$A <- 1
#'
#'     g <- nadir::lnr_logistic(data, I(A == 1) ~ L)
#'     return(list(
#'       m = m(data),
#'       ma = m(dataAeq_a),
#'       g = g(data)))
#'   })
#'
#' rc_cfmean_0 <- riesz_curve_catalog$cfmean_a(
#'   A = 'A',
#'   a = 0,
#'
#'   nuis = function(data) {
#'     m <- nadir::lnr_glm(data, Y ~ L + A)
#'     dataAeq_a <- data
#'     dataAeq_a$A <- 0
#'
#'     g <- nadir::lnr_logistic(data, I(A == 0) ~ L)
#'     return(list(
#'       m = m(data),
#'       ma = m(dataAeq_a),
#'       g = g(data)))
#'   }
#' )
#'
#' n <- 10000
#' df <- tibble::tibble(
#'   L = rnorm(n = n),
#'   A = rbinom(n = n, size = 1, prob = plogis(L)),
#'   Y = L + A + rnorm(n = n, sd = .1))
#'
#' riesz_estimate(data = df, rc = rc_cfmean_0)
#' riesz_estimate(data = df, rc = rc_cfmean_1)
riesz_curve_catalog <- list()

riesz_curve_catalog$cfmean_a <- function(
    a,
    nuis) {

  alpha_formula <- as.formula(paste0('~ I(A ==', a, ')/g'))

  return(RieszCurve$new(
    nuis = nuis,
    alpha = alpha_formula,
    f = ~ m,
    h = ~ ma,
    ic_expr = ~ h + alpha * (Y - f)
  ))
}


riesz_curve_catalog$cfmean_a1 <- function(nuis) {
  riesz_curve_catalog$cfmean_a(a = 1, nuis = nuis)
}

riesz_curve_catalog$cfmean_a0 <- function(nuis) {
  riesz_curve_catalog$cfmean_a(a = 0, nuis = nuis)
}


riesz_curve_catalog$ate <- function(nuis) {

  RieszCurve$new(
    nuis = nuis,
    alpha = ~ A/g + (1 - A)/(1 - g),
    f = ~ m,
    h = ~ m1 - m0,
    ic_expr = ~ h + alpha * (Y - f)
  )
}

riesz_curve_catalog$att <- function(nuis) {

  RieszCurve$new(
    nuis = nuis,
    alpha = ~ A + (A - g)/(g * (1 - g)),
    f = ~ m,
    h = ~ m1 - m0,
    ic_expr = ~ A * (Y - m0) / mean(A) + (m1 - m0) - h
  )
}

riesz_curve_catalog$subgroup_mean <- function(v, nuis) {

  alpha_formula <- as.formula(
    paste0("~ I(V == ", v, ") / pv")
  )

  RieszCurve$new(
    nuis = nuis,
    alpha = alpha_formula,
    f = ~ m,
    h = ~ mv,
    ic_expr = ~ h + alpha * (Y - f)
  )
}

riesz_curve_catalog$missing_mean <- function(nuis) {

  RieszCurve$new(
    nuis = nuis,
    alpha = ~ R/pi,
    f = ~ m,
    h = ~ m,
    ic_expr = ~ h + alpha * (Y - f)
  )
}

riesz_curve_catalog$mtp <- function(nuis) {

  RieszCurve$new(
    nuis = nuis,
    alpha = ~ g_shift / g,
    f = ~ m,
    h = ~ m_shift,
    ic_expr = ~ h + alpha * (Y - f)
  )
}





#' Logit and expit (inverse logit) transformations
#'
#' `logit(p)` computes \eqn{\log\{p / (1 - p)\}} and `expit(x)` computes
#' \eqn{1 / (1 + e^{-x})}. Thin wrappers around [stats::qlogis()] and
#' [stats::plogis()], used by the logistic fluctuation in [riesz_tmle()].
#'
#' @param p Numeric vector of probabilities in (0, 1).
#' @param x Numeric vector on the logit scale.
#' @return A numeric vector.
#' @name logit
#' @export
#' @examples
#' expit(logit(0.3))
logit <- function(p) stats::qlogis(p)

#' @rdname logit
#' @export
expit <- function(x) stats::plogis(x)

to01 <- function(x, bounds) {
  if (is.null(bounds)) stop("`bounds` must be provided for logistic fluctuation.")
  if (!is.numeric(bounds) || length(bounds) != 2L) stop("`bounds` must be a numeric vector of length 2.")
  a <- bounds[[1]]; b <- bounds[[2]]
  if (!is.finite(a) || !is.finite(b)) stop("`bounds` must be finite.")
  if (b <= a) stop("`bounds[2]` must be > `bounds[1]`.")

  (x - a) / (b - a)
}

from01 <- function(p, bounds) {
  if (is.null(bounds)) stop("`bounds` must be provided.")
  if (!is.numeric(bounds) || length(bounds) != 2L) stop("`bounds` must be a numeric vector of length 2.")
  a <- bounds[[1]]; b <- bounds[[2]]
  if (b <= a) stop("`bounds[2]` must be > `bounds[1]`.")

  a + p * (b - a)
}

clip01 <- function(p, clip = 1e-6) {
  if (!is.numeric(clip) || length(clip) != 1L || is.na(clip) || clip <= 0 || clip >= 0.5) {
    stop("`clip` must be a numeric scalar in (0, 0.5).")
  }
  pmin(pmax(p, clip), 1 - clip)
}

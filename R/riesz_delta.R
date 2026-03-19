.numeric_gradient <- function(g, theta, eps = 1e-6) {
  if (!is.function(g)) stop("`g` must be a function.")
  if (!is.numeric(theta) || length(theta) < 1L) {
    stop("`theta` must be a numeric vector of length >= 1.")
  }
  if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) || eps <= 0) {
    stop("`eps` must be a positive numeric scalar.")
  }

  k <- length(theta)
  grad <- numeric(k)

  for (j in seq_len(k)) {
    step <- rep(0, k)
    step[j] <- eps
    grad[j] <- (g(theta + step) - g(theta - step)) / (2 * eps)
  }

  grad
}

.validate_riesz_fits <- function(fits) {
  if (!is.list(fits) || length(fits) < 1L) {
    stop("`fits` must be a non-empty list of `RieszFit` objects.")
  }

  ok <- vapply(fits, function(x) inherits(x, "RieszFit"), logical(1))
  if (!all(ok)) {
    stop("All elements of `fits` must inherit from `RieszFit`.")
  }

  n_vec <- vapply(fits, function(x) length(x$ic), integer(1))
  if (length(unique(n_vec)) != 1L) {
    stop("All fits must have influence curves of the same length.")
  }

  invisible(NULL)
}

.extract_component_estimates <- function(fits) {
  vapply(fits, function(x) x$estimate, numeric(1))
}

.extract_component_eif_matrix <- function(fits) {
  # Use centered EIFs. Since `ic` is stored as an uncentered IC-like quantity,
  # center it at the reported estimate.
  eif_mat <- do.call(
    cbind,
    lapply(fits, function(x) {
      if (is.null(x$ic)) {
        stop("Each `RieszFit` must contain an `ic` component.")
      }
      as.numeric(x$ic - x$estimate)
    })
  )

  colnames(eif_mat) <- paste0("theta", seq_len(ncol(eif_mat)))
  eif_mat
}

#' @export
riesz_delta <- function(...,
                        fits = NULL,
                        g,
                        grad = NULL,
                        significance_alpha = 0.05,
                        numgrad_eps = 1e-6,
                        parameter = NULL,
                        estimator = "delta-method") {
  if (is.null(fits)) {
    fits <- list(...)
  } else {
    if (length(list(...)) > 0L) {
      stop("Supply fits either through `...` or through `fits =`, not both.")
    }
  }

  .validate_riesz_fits(fits)

  if (!is.function(g)) {
    stop("`g` must be a function.")
  }

  if (!is.null(grad) && !is.function(grad)) {
    stop("`grad` must be NULL or a function.")
  }

  if (!is.numeric(significance_alpha) || length(significance_alpha) != 1L ||
      is.na(significance_alpha) || significance_alpha <= 0 || significance_alpha >= 1) {
    stop("`significance_alpha` must be a numeric scalar in (0, 1).")
  }

  theta_hat <- .extract_component_estimates(fits)
  eif_mat <- .extract_component_eif_matrix(fits)

  psi_hat <- g(theta_hat)
  if (!is.numeric(psi_hat) || length(psi_hat) != 1L || is.na(psi_hat)) {
    stop("`g(theta_hat)` must return a non-missing numeric scalar.")
  }

  grad_hat <- if (is.null(grad)) {
    .numeric_gradient(g, theta_hat, eps = numgrad_eps)
  } else {
    grad(theta_hat)
  }

  if (!is.numeric(grad_hat) || length(grad_hat) != length(theta_hat)) {
    stop("Gradient must be a numeric vector of length equal to the number of fits.")
  }

  eif_psi <- as.numeric(eif_mat %*% grad_hat)
  ic_psi <- psi_hat + eif_psi

  var_hat <- stats::var(eif_psi) / nrow(eif_mat)
  se_hat <- sqrt(var_hat)
  ci_low <- psi_hat + stats::qnorm(significance_alpha / 2) * se_hat
  ci_high <- psi_hat + stats::qnorm(1 - significance_alpha / 2) * se_hat

  RieszFit$new(
    estimate = psi_hat,
    var = var_hat,
    se = se_hat,
    ci_low = ci_low,
    ci_high = ci_high,
    significance_alpha = significance_alpha,
    estimator = estimator,
    parameter = parameter,
    n = nrow(eif_mat),
    ic = ic_psi
  )
}

#' @export
riesz_delta_difference <- function(fit1,
                                   fit2,
                                   significance_alpha = 0.05,
                                   parameter = "difference") {
  riesz_delta(
    fit1, fit2,
    g = function(theta) theta[1] - theta[2],
    grad = function(theta) c(1, -1),
    significance_alpha = significance_alpha,
    parameter = parameter,
    estimator = "delta-method"
  )
}

#' @export
riesz_delta_ratio <- function(fit1,
                              fit2,
                              significance_alpha = 0.05,
                              parameter = "ratio",
                              tol = 1e-8) {
  if (abs(fit2$estimate) < tol) {
    stop("Cannot form ratio: denominator estimate is too close to zero.")
  }

  riesz_delta(
    fit1, fit2,
    g = function(theta) theta[1] / theta[2],
    grad = function(theta) c(1 / theta[2], -theta[1] / theta[2]^2),
    significance_alpha = significance_alpha,
    parameter = parameter,
    estimator = "delta-method"
  )
}

#' @export
riesz_delta_log_ratio <- function(fit1,
                                  fit2,
                                  significance_alpha = 0.05,
                                  parameter = "log-ratio",
                                  tol = 1e-8) {
  if (fit1$estimate <= tol || fit2$estimate <= tol) {
    stop("Log-ratio requires both estimates to be strictly positive.")
  }

  riesz_delta(
    fit1, fit2,
    g = function(theta) log(theta[1] / theta[2]),
    grad = function(theta) c(1 / theta[1], -1 / theta[2]),
    significance_alpha = significance_alpha,
    parameter = parameter,
    estimator = "delta-method"
  )
}

#' @export
riesz_delta_odds_ratio <- function(fit1,
                                   fit2,
                                   significance_alpha = 0.05,
                                   parameter = "odds-ratio",
                                   tol = 1e-8) {
  p1 <- fit1$estimate
  p0 <- fit2$estimate

  if (p1 <= tol || p1 >= 1 - tol || p0 <= tol || p0 >= 1 - tol) {
    stop("Odds ratio requires both estimates to lie strictly inside (0, 1).")
  }

  riesz_delta(
    fit1, fit2,
    g = function(theta) {
      odds1 <- theta[1] / (1 - theta[1])
      odds0 <- theta[2] / (1 - theta[2])
      odds1 / odds0
    },
    grad = function(theta) {
      or_hat <- (theta[1] / (1 - theta[1])) / (theta[2] / (1 - theta[2]))
      c(
        or_hat / (theta[1] * (1 - theta[1])),
        -or_hat / (theta[2] * (1 - theta[2]))
      )
    },
    significance_alpha = significance_alpha,
    parameter = parameter,
    estimator = "delta-method"
  )
}

#' @export
riesz_delta_log_odds_ratio <- function(fit1,
                                       fit2,
                                       significance_alpha = 0.05,
                                       parameter = "log-odds-ratio",
                                       tol = 1e-8) {
  p1 <- fit1$estimate
  p0 <- fit2$estimate

  if (p1 <= tol || p1 >= 1 - tol || p0 <= tol || p0 >= 1 - tol) {
    stop("Log odds ratio requires both estimates to lie strictly inside (0, 1).")
  }

  riesz_delta(
    fit1, fit2,
    g = function(theta) {
      log(theta[1] / (1 - theta[1])) - log(theta[2] / (1 - theta[2]))
    },
    grad = function(theta) {
      c(
        1 / (theta[1] * (1 - theta[1])),
        -1 / (theta[2] * (1 - theta[2]))
      )
    },
    significance_alpha = significance_alpha,
    parameter = parameter,
    estimator = "delta-method"
  )
}

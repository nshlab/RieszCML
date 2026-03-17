riesz_tmle <- function(data,
                       rc,
                       fluctuation = "logistic",
                       bounds = NULL,
                       outcome_col,
                       use_intercept_and_weights = TRUE,
                       clip = 1e-6,
                       significance_alpha = 0.05) {

  .validate_riesz_tmle_inputs(
    data = data,
    rc = rc,
    fluctuation = fluctuation,
    bounds = bounds,
    outcome_col = outcome_col,
    clip = clip,
    significance_alpha = significance_alpha
  )

  if (inherits(rc, "RieszCurve")) {
    return(.riesz_tmle_single(
      data = data,
      rc = rc,
      fluctuation = fluctuation,
      bounds = bounds,
      outcome_col = outcome_col,
      use_intercept_and_weights = use_intercept_and_weights,
      clip = clip,
      significance_alpha = significance_alpha
    ))
  }

  if (inherits(rc, "ComposedRieszCurve")) {
    return(.riesz_tmle_composed(
      data = data,
      rc = rc,
      fluctuation = fluctuation,
      bounds = bounds,
      outcome_col = outcome_col,
      use_intercept_and_weights = use_intercept_and_weights,
      clip = clip,
      significance_alpha = significance_alpha
    ))
  }

  stop("`rc` must inherit from `RieszCurve` or `ComposedRieszCurve`.")
}


.validate_riesz_tmle_inputs <- function(data,
                                        rc,
                                        fluctuation,
                                        bounds,
                                        outcome_col,
                                        clip,
                                        significance_alpha) {

  if (!inherits(rc, "RieszCurve") && !inherits(rc, "ComposedRieszCurve")) {
    stop("`rc` must be a RieszCurve or ComposedRieszCurve.")
  }

  if (!is.character(outcome_col) || length(outcome_col) != 1L ||
      is.na(outcome_col) || outcome_col == "") {
    stop("`outcome_col` must be a non-empty character scalar.")
  }

  if (!outcome_col %in% names(data)) {
    stop("`outcome_col` is not a column in `data`.")
  }

  if (!fluctuation %in% c("identity", "logistic")) {
    stop("`fluctuation` must be one of {'identity', 'logistic'}.")
  }

  if (!is.numeric(clip) || length(clip) != 1L || is.na(clip) ||
      clip <= 0 || clip >= 0.5) {
    stop("`clip` must be a numeric scalar in (0, 0.5).")
  }

  if (!is.numeric(significance_alpha) || length(significance_alpha) != 1L ||
      is.na(significance_alpha) || significance_alpha <= 0 || significance_alpha >= 1) {
    stop("`significance_alpha` must be a numeric scalar in (0, 1).")
  }

  if (fluctuation == "logistic") {
    if (is.null(bounds)) {
      stop("For `fluctuation = 'logistic'`, `bounds` must be supplied.")
    }
    if (!is.numeric(bounds) || length(bounds) != 2L) {
      stop("`bounds` must be a numeric vector of length 2.")
    }
    if (!is.finite(bounds[[1]]) || !is.finite(bounds[[2]]) || bounds[[2]] <= bounds[[1]]) {
      stop("`bounds` must satisfy bounds[1] < bounds[2], both finite.")
    }
  }

  invisible(NULL)
}

.fit_tmle_fluctuation <- function(target,
                                  offset,
                                  clever_covariate,
                                  fluctuation,
                                  bounds = NULL,
                                  clip = 1e-6,
                                  use_intercept_and_weights = TRUE) {

  n <- length(target)

  if (!is.numeric(target) || length(target) != n) stop("`target` must be numeric.")
  if (!is.numeric(offset) || length(offset) != n) stop("`offset` must be numeric.")
  if (!is.numeric(clever_covariate) || length(clever_covariate) != n) {
    stop("`clever_covariate` must be numeric.")
  }

  if (fluctuation == "identity") {
    dat <- data.frame(
      target = target,
      H = clever_covariate,
      offset = offset
    )

    if (use_intercept_and_weights) {
      fit <- stats::glm(
        target ~ 1,
        family = stats::gaussian(),
        data = dat,
        weights = dat$H,
        offset = offset
      )
      coefs <- stats::coef(fit)
      intercept <- 0
      eps <- unname(coefs["(Intercept)"])
      if (is.na(intercept)) intercept <- 0
    } else {
      fit <- stats::glm(
        target ~ -1 + H,
        family = stats::gaussian(),
        data = dat,
        offset = offset
      )
      coefs <- stats::coef(fit)
      intercept <- 0
      eps <- unname(coefs["H"])
    }

    updated <- offset + intercept + eps * clever_covariate

    return(list(
      eps = eps,
      intercept = intercept,
      updated = updated,
      fit = fit
    ))
  }

  if (fluctuation == "logistic") {
    target01 <- clip01(to01(target, bounds), clip = clip)
    offset01 <- clip01(to01(offset, bounds), clip = clip)
    offset_logit <- logit(offset01)

    dat <- data.frame(
      target01 = target01,
      H = clever_covariate,
      offset_logit = offset_logit
    )

    if (use_intercept_and_weights) {
      fit <- suppressWarnings(stats::glm(
        target01 ~ 1 + H,
        family = stats::binomial(),
        data = dat,
        offset = offset_logit
      ))
      coefs <- stats::coef(fit)
      intercept <- 0
      eps <- unname(coefs["(Intercept)"])
      if (is.na(intercept)) intercept <- 0
    } else {
      fit <- suppressWarnings(stats::glm(
        target01 ~ -1 + H,
        family = stats::binomial(),
        data = dat,
        offset = offset_logit
      ))
      coefs <- stats::coef(fit)
      intercept <- 0
      eps <- unname(coefs["H"])
    }

    updated01 <- expit(offset_logit + intercept + eps * clever_covariate)
    updated01 <- clip01(updated01, clip = clip)
    updated <- from01(updated01, bounds)

    return(list(
      eps = eps,
      intercept = intercept,
      updated = updated,
      fit = fit
    ))
  }

  stop("Unsupported fluctuation.")
}

.summarize_tmle_fit <- function(estimate, ic_for_inference, significance_alpha = 0.05) {
  if (!is.numeric(estimate) || length(estimate) != 1L || is.na(estimate)) {
    stop("`estimate` must be a non-missing numeric scalar.")
  }
  if (!is.numeric(ic_for_inference) || length(ic_for_inference) < 2L) {
    stop("`ic_for_inference` must be a numeric vector of length at least 2.")
  }

  n <- length(ic_for_inference)
  var_estimate <- stats::var(ic_for_inference) / n
  se <- sqrt(var_estimate)
  ci_low <- estimate + stats::qnorm(significance_alpha / 2) * se
  ci_high <- estimate + stats::qnorm(1 - significance_alpha / 2) * se

  list(
    estimate = estimate,
    var = var_estimate,
    se = se,
    ci_low = ci_low,
    ci_high = ci_high
  )
}


.apply_fluctuation_map <- function(x,
                                   clever_covariate,
                                   eps,
                                   intercept = 0,
                                   fluctuation,
                                   bounds = NULL,
                                   clip = 1e-6) {
  if (!is.numeric(x) || !is.numeric(clever_covariate)) {
    stop("`x` and `clever_covariate` must be numeric.")
  }
  if (length(x) != length(clever_covariate)) {
    stop("`x` and `clever_covariate` must have the same length.")
  }

  if (fluctuation == "identity") {
    return(x + intercept + eps * clever_covariate)
  }

  if (fluctuation == "logistic") {
    x01 <- clip01(to01(x, bounds), clip = clip)
    x01_star <- expit(logit(x01) + intercept + eps * clever_covariate)
    x01_star <- clip01(x01_star, clip = clip)
    return(from01(x01_star, bounds))
  }

  stop("Unsupported `fluctuation`.")
}


.riesz_tmle_single <- function(data,
                               rc,
                               fluctuation,
                               bounds,
                               outcome_col,
                               use_intercept_and_weights,
                               clip,
                               significance_alpha) {

  y <- data[[outcome_col]]
  rc$fit(data)

  alpha <- rc$fit_alpha
  h <- rc$fit_h
  f <- rc$fit_f

  fluc <- .fit_tmle_fluctuation(
    target = y,
    offset = f,
    clever_covariate = alpha,
    fluctuation = fluctuation,
    bounds = bounds,
    clip = clip,
    use_intercept_and_weights = use_intercept_and_weights
  )

  f_star <- fluc$updated

  h_star <- .apply_fluctuation_map(
    x = h,
    clever_covariate = alpha,
    eps = fluc$eps,
    intercept = fluc$intercept,
    fluctuation = fluctuation,
    bounds = bounds,
    clip = clip
  )

  tmle_estimate <- mean(h_star)

  # inference IC: untargeted nuisance parts plugged in,
  # but inference centered around the TMLE plug-in estimate
  ic_for_inference <- h_star + alpha * (y - f)

  out <- .summarize_tmle_fit(
    estimate = tmle_estimate,
    ic_for_inference = ic_for_inference,
    significance_alpha = significance_alpha
  )

  out$eps <- fluc$eps
  out$intercept <- fluc$intercept
  out$ic <- ic_for_inference
  out$f_star <- f_star
  out$h_star <- h_star
  out$fluctuation_model <- fluc$fit
  out
}


.riesz_tmle_composed <- function(data,
                                 rc,
                                 fluctuation,
                                 bounds,
                                 outcome_col,
                                 use_intercept_and_weights,
                                 clip,
                                 significance_alpha) {

  y <- data[[outcome_col]]
  n <- nrow(data)

  rc$fit(data)
  J <- length(rc$rc_list)

  alpha_list <- vector("list", J)
  f_list <- vector("list", J)
  h_list <- vector("list", J)
  omega_list <- vector("list", J)

  omega_running <- rep(1, n)

  for (j in seq_len(J)) {
    rj <- rc$rc_list[[j]]
    alpha_list[[j]] <- rj$fit_alpha
    f_list[[j]] <- rj$fit_f
    h_list[[j]] <- rj$fit_h

    omega_running <- omega_running * alpha_list[[j]]
    omega_list[[j]] <- omega_running
  }

  target_list <- vector("list", J)
  f_star_list <- vector("list", J)
  h_star_list <- vector("list", J)
  eps_vec <- numeric(J)
  intercept_vec <- numeric(J)
  fluctuation_models <- vector("list", J)

  current_target <- y

  for (j in seq_len(J)) {
    target_list[[j]] <- current_target

    fluc <- .fit_tmle_fluctuation(
      target = current_target,
      offset = f_list[[j]],
      clever_covariate = omega_list[[j]],
      fluctuation = fluctuation,
      bounds = bounds,
      clip = clip,
      use_intercept_and_weights = use_intercept_and_weights
    )

    f_star_list[[j]] <- fluc$updated
    eps_vec[[j]] <- fluc$eps
    intercept_vec[[j]] <- fluc$intercept
    fluctuation_models[[j]] <- fluc$fit

    h_star_list[[j]] <- .apply_fluctuation_map(
      x = h_list[[j]],
      clever_covariate = omega_list[[j]],
      eps = fluc$eps,
      intercept = fluc$intercept,
      fluctuation = fluctuation,
      bounds = bounds,
      clip = clip
    )

    current_target <- h_star_list[[j]]
  }

  tmle_estimate <- mean(h_star_list[[J]])

  # inference IC: untargeted residual pieces, targeted final plug-in
  ic_for_inference <- rep(0, n)
  for (j in seq_len(J)) {
    ic_for_inference <- ic_for_inference +
      omega_list[[j]] * (target_list[[j]] - f_list[[j]])
  }
  ic_for_inference <- ic_for_inference + h_star_list[[J]]

  out <- .summarize_tmle_fit(
    estimate = tmle_estimate,
    ic_for_inference = ic_for_inference,
    significance_alpha = significance_alpha
  )

  out$eps <- eps_vec
  out$intercept <- intercept_vec
  out$ic <- ic_for_inference
  out$omega <- omega_list
  out$targets <- target_list
  out$f_star <- f_star_list
  out$h_star <- h_star_list
  out$fluctuation_models <- fluctuation_models
  out
}



# old 2026-03-17 ---------------------------------------------------------------------
#
# riesz_tmle <- function(data,
#                        rc,
#                        fluctuation = 'logistic',
#                        bounds,
#                        outcome_col,
#                        use_intercept_and_weights = TRUE,
#                        clip = 1e-6,
#                        significance_alpha = .05) {
#
#   # checks
#   if (!inherits(rc, "RieszCurve")) stop("`rc` must be a RieszCurve.")
#   if (!is.character(outcome_col) || length(outcome_col) != 1L || outcome_col == "" || is.na(outcome_col)) {
#     stop("`outcome_col` must be a non-empty character scalar.")
#   }
#   if (!outcome_col %in% names(data)) stop("`outcome_col` is not a column in `data`.")
#
#
#   # if identity, defer to riesz_estimate
#   if (fluctuation == 'identity') {
#     return(riesz_estimate(data, rc))
#   }
#
#   # fluctuation on the logistic scale
#   if (fluctuation == 'logistic') {
#
#     # bounds checks
#     if (is.null(bounds)) stop("For `fluctuation = 'logistic'`, you must provide `bounds = c(a, b)`.")
#     if (!is.numeric(bounds) || length(bounds) != 2L) stop("`bounds` must be numeric length 2: c(a, b).")
#     a <- bounds[[1]]; b <- bounds[[2]]
#     if (!is.finite(a) || !is.finite(b)) stop("`bounds` must be finite.")
#     if (b <= a) stop("`bounds[2]` must be > `bounds[1]`.")
#
#     # fit rc pieces
#     rc$fit(data)
#
#     if (is.null(rc$fit_alpha)) stop("After rc$fit(data), rc$fit_alpha is NULL.")
#     if (is.null(rc$fit_h)) stop("After rc$fit(data), rc$fit_h is NULL.")
#     if (is.null(rc$fit_f)) stop("After rc$fit(data), rc$fit_f is NULL (need the regression f).")
#
#     alpha <- rc$fit_alpha
#     h <- rc$fit_h
#     f <- rc$fit_f
#     y <- data[[outcome_col]]
#     n <- nrow(data)
#
#     if (!is.numeric(alpha) || length(alpha) != n) stop("`alpha` must be numeric length nrow(data).")
#     if (!is.numeric(h) || length(h) != n) stop("`h` must be numeric length nrow(data).")
#     if (!is.numeric(f) || length(f) != n) stop("`f` must be numeric length nrow(data).")
#     if (!is.numeric(y) || length(y) != n) stop("Outcome must be numeric length nrow(data).")
#
#     # map to [0,1]
#     y01 <- to01(y, bounds)
#     f01 <- to01(f, bounds)
#
#     # range checks
#     if (any(!is.finite(y01))) stop("Outcome mapped to [0,1] produced non-finite values.")
#     if (any(!is.finite(f01))) stop("f mapped to [0,1] produced non-finite values.")
#
#     # clip away from 0/1 for logit
#     y01 <- clip01(y01, clip = clip)
#     f01 <- clip01(f01, clip = clip)
#
#     # Fit fluctuation: logit(f*) = logit(f) + eps * alpha
#     tmle_data <- data.frame(
#       y01 = y01,
#       alpha = alpha,
#       offset = logit(f01)
#     )
#
#     if (use_intercept_and_weights) {
#       tmle_formula <- y01 ~ 1
#     } else {
#       tmle_formula <- y01 ~ -1 + alpha
#     }
#
#     fluctuation_model <- suppressWarnings(stats::glm(
#       formula = tmle_formula,
#       family = stats::binomial(),
#       data = tmle_data,
#       offset = tmle_data$offset,
#       weights = if (use_intercept_and_weights) tmle_data[['alpha']] else NULL
#     )) # tends to message that non-integer #successes were given in a
#     # binomial glm or that glm.fit fitted probabilities numerically 0 or 1
#
#     if (use_intercept_and_weights) {
#       eps <- as.numeric(stats::coef(fluctuation_model)[['(Intercept)']])
#     } else {
#       eps <- as.numeric(stats::coef(fluctuation_model)[["alpha"]])
#     }
#     if (!is.finite(eps)) stop("Fluctuation coefficient `eps` is non-finite.")
#
#     f01_star <- expit(logit(f01) + eps * alpha)
#     f01_star <- clip01(f01_star, clip = clip)
#
#     # map back to [a,b]
#     f_star <- from01(f01_star, bounds)
#
#     # updated IC contribution (uncentered)
#     ic_star <- h + alpha * (y - f_star)
#
#     estimate <- mean(ic_star)
#     var_estimate <- stats::var(ic_star) / n
#     se <- sqrt(var_estimate)
#     ci_low <- estimate + stats::qnorm(significance_alpha / 2) * se
#     ci_high <- estimate + stats::qnorm(1 - significance_alpha / 2) * se
#
#     return(list(
#       estimate = estimate,
#       var = var_estimate,
#       se = se,
#       ci_low = ci_low,
#       ci_high = ci_high,
#       eps = eps,
#       ic = ic_star,
#       f_star = f_star
#     ))
#   }
# }

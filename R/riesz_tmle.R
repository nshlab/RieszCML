#' @export
riesz_tmle <- function(data,
                       rc,
                       fluctuation = "logistic",
                       bounds = NULL,
                       outcome_col,
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
      clip = clip,
      significance_alpha = significance_alpha
    ))
  }

  stop("`rc` must be a RieszCurve or ComposedRieszCurve.")
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
                                  clip = 1e-6) {

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


    fit <- stats::glm(
      target ~ -1 + H,
      family = stats::gaussian(),
      data = dat,
      offset = offset
    )
    coefs <- stats::coef(fit)
    intercept <- 0
    eps <- unname(coefs["H"])

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

    fit <- suppressWarnings(stats::glm(
        target01 ~ -1 + H,
        family = stats::binomial(),
        data = dat,
        offset = offset_logit
      ))
    coefs <- stats::coef(fit)
    eps <- unname(coefs["H"])

    updated01 <- expit(offset_logit + eps * clever_covariate)
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
    ci_high = ci_high,
    n = n
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
                               clip,
                               significance_alpha) {

  y <- data[[outcome_col]]
  rc$fit(data)

  f <- rc$fit_f
  h <- rc$fit_h
  alpha <- rc$fit_alpha


  # untargeted procedure for EIF --------------------------------------------------------------

  # Untargeted plug-in
  psi_init <- mean(h)

  # Untargeted IC-like quantity
  ic_init <- rc$fit_ic

  # Centered EIF for variance
  eif_init <- ic_init - psi_init

  # targeted procedure for estimate -----------------------------------------

  alpha <- rc$fit_alpha
  h <- rc$fit_h
  f <- rc$fit_f

  fluc <- .fit_tmle_fluctuation(
    target = y,
    offset = f,
    clever_covariate = alpha,
    fluctuation = fluctuation,
    bounds = bounds,
    clip = clip
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

  summarized_tmle_fit <- .summarize_tmle_fit(
    estimate = tmle_estimate,
    ic_for_inference = eif_init,
    significance_alpha = significance_alpha
  )

  RieszFit$new(
    estimate = summarized_tmle_fit$estimate,
    var = summarized_tmle_fit$var,
    se = summarized_tmle_fit$se,
    ci_low = summarized_tmle_fit$ci_low,
    ci_high = summarized_tmle_fit$ci_high,
    significance_alpha = significance_alpha,
    estimator = "TMLE",
    n = summarized_tmle_fit$n,
    ic = eif_init,
    eps = fluc$eps,
    intercept = fluc$intercept,
    f_star = f_star,
    h_star = h_star,
    fluctuation_model = fluc$fit
  )
}

.riesz_tmle_composed <- function(data,
                                 rc,
                                 fluctuation,
                                 bounds,
                                 outcome_col,
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
      clip = clip
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

  untargeted_ic <- .compute_untargeted_ic_composed(
    y = y,
    omega_list = omega_list,
    f_list = f_list,
    h_list = h_list
  )
  ic_for_inference <- untargeted_ic$ic
  plugin_estimate <- untargeted_ic$psi

  summarized_tmle_fit <- .summarize_tmle_fit(
    estimate = tmle_estimate,
    ic_for_inference = ic_for_inference,
    significance_alpha = significance_alpha
  )

  RieszFit$new(
    estimate = summarized_tmle_fit$estimate,
    var = summarized_tmle_fit$var,
    se = summarized_tmle_fit$se,
    ci_low = summarized_tmle_fit$ci_low,
    ci_high = summarized_tmle_fit$ci_high,
    significance_alpha = significance_alpha,
    estimator = "TMLE",
    n = summarized_tmle_fit$n,
    eps = eps_vec,
    intercept = intercept_vec,
    ic = ic_for_inference,
    f_star = f_star_list,
    h_star = h_star_list,
    fluctuation_models = fluctuation_models,
    omega = omega_list,
    targets = target_list
  )
}


.compute_untargeted_ic_composed <- function(y, omega_list, f_list, h_list) {
  n <- length(y)
  J <- length(f_list)

  current_target_init <- y
  target_list_init <- vector("list", J)

  for (j in seq_len(J)) {
    target_list_init[[j]] <- current_target_init
    current_target_init <- h_list[[j]]
  }

  ic_init <- rep(0, n)
  for (j in seq_len(J)) {
    ic_init <- ic_init + omega_list[[j]] * (target_list_init[[j]] - f_list[[j]])
  }
  ic_init <- ic_init + h_list[[J]]

  list(
    ic = ic_init,
    psi = mean(h_list[[J]]),
    targets = target_list_init
  )
}



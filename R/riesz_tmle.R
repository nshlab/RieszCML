#' Targeted Minimum Loss-Based Estimation for Riesz-Represented Functionals
#'
#' Constructs a TMLE for a linear functional represented by a [RieszCurve]
#' (Algorithm 1 of the manuscript with a single stage) or a
#' [ComposedRieszCurve] (sequential-regression TMLE, Algorithm 1).
#'
#' @param data A `data.frame`-like object containing the observed data.
#' @param rc A [RieszCurve] or [ComposedRieszCurve]. For valid TMLE updates,
#'   each curve must specify how its counterfactual predictions are updated:
#'   via `targeting_steps` (single curves) and via `alpha_star` (each stage
#'   of a composed curve whose `h` differs from its `f`).
#' @param fluctuation_type One of `"logistic"` (recommended; requires
#'   `bounds`) or `"identity"`.
#' @param bounds Numeric vector of length 2 bounding the outcome scale, used
#'   by the logistic fluctuation. In the composed case the same bounds are
#'   applied to all intermediate targets, which are guaranteed to respect
#'   them by construction of the logistic update.
#' @param outcome_col Name of the outcome column in `data`.
#' @param fluctuation_weights Optional nonnegative observation weights for
#'   the fluctuation regressions.
#' @param clip Clipping constant for the logistic fluctuation, in (0, 0.5).
#' @param significance_alpha Significance level for Wald confidence
#'   intervals, based on the variance of the estimated (untargeted) EIF.
#' @param score_tol Multiple of the estimated standard error beyond which a
#'   nonzero empirical mean of the targeted EIF estimating equation triggers
#'   a warning. The targeted score should be numerically zero after a
#'   successful TMLE update; a large value indicates a failed or
#'   mis-specified targeting step.
#'
#' @return A [RieszFit]. `ic` holds the untargeted (inference) influence
#'   curve; `ic_star` holds the targeted influence curve, whose mean equals
#'   the TMLE estimate up to numerical error when the targeting step
#'   succeeded.
#' @export
riesz_tmle <- function(data,
                       rc,
                       fluctuation_type = "logistic",
                       bounds = NULL,
                       outcome_col,
                       fluctuation_weights = NULL,
                       clip = 1e-6,
                       significance_alpha = 0.05,
                       score_tol = 0.05) {

  .validate_riesz_tmle_inputs(
    data = data,
    rc = rc,
    fluctuation_type = fluctuation_type,
    bounds = bounds,
    outcome_col = outcome_col,
    fluctuation_weights = fluctuation_weights,
    clip = clip,
    significance_alpha = significance_alpha
  )

  if (inherits(rc, "RieszCurve")) {
    return(.riesz_tmle_single(
      data = data,
      rc = rc,
      fluctuation_type = fluctuation_type,
      bounds = bounds,
      outcome_col = outcome_col,
      fluctuation_weights = fluctuation_weights,
      clip = clip,
      significance_alpha = significance_alpha,
      score_tol = score_tol
    ))
  }

  if (inherits(rc, "ComposedRieszCurve")) {
    return(.riesz_tmle_composed(
      data = data,
      rc = rc,
      fluctuation_type = fluctuation_type,
      bounds = bounds,
      outcome_col = outcome_col,
      fluctuation_weights = fluctuation_weights,
      clip = clip,
      significance_alpha = significance_alpha,
      score_tol = score_tol
    ))
  }

  stop("`rc` must be a RieszCurve or ComposedRieszCurve.")
}


.validate_riesz_tmle_inputs <- function(data,
                                        rc,
                                        fluctuation_type,
                                        bounds,
                                        outcome_col,
                                        fluctuation_weights,
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

  if (!fluctuation_type %in% c("identity", "logistic")) {
    stop("`fluctuation_type` must be one of {'identity', 'logistic'}.")
  }

  if (!is.numeric(clip) || length(clip) != 1L || is.na(clip) ||
      clip <= 0 || clip >= 0.5) {
    stop("`clip` must be a numeric scalar in (0, 0.5).")
  }

  if (!is.numeric(significance_alpha) || length(significance_alpha) != 1L ||
      is.na(significance_alpha) || significance_alpha <= 0 || significance_alpha >= 1) {
    stop("`significance_alpha` must be a numeric scalar in (0, 1).")
  }

  if (fluctuation_type == "logistic") {
    if (is.null(bounds)) {
      stop("For `fluctuation_type = 'logistic'`, `bounds` must be supplied.")
    }
    if (!is.numeric(bounds) || length(bounds) != 2L) {
      stop("`bounds` must be a numeric vector of length 2.")
    }
    if (!is.finite(bounds[[1]]) || !is.finite(bounds[[2]]) || bounds[[2]] <= bounds[[1]]) {
      stop("`bounds` must satisfy bounds[1] < bounds[2], both finite.")
    }
  }
  if (!is.null(fluctuation_weights)) {
    if (!is.numeric(fluctuation_weights)) {
      stop("`fluctuation_weights` must be NULL or numeric.")
    }
    if (length(fluctuation_weights) != nrow(data)) {
      stop("`fluctuation_weights` must have length nrow(data).")
    }
    if (anyNA(fluctuation_weights)) {
      stop("`fluctuation_weights` cannot contain NA values.")
    }
    if (any(fluctuation_weights < 0)) {
      stop("`fluctuation_weights` must be nonnegative.")
    }
  }

  invisible(NULL)
}

.fit_tmle_fluctuation <- function(target,
                                  offset,
                                  clever_covariate,
                                  fluctuation_type,
                                  fluctuation_weights = NULL,
                                  bounds = NULL,
                                  clip = 1e-6) {

  n <- length(target)

  if (!is.numeric(target) || length(target) != n) stop("`target` must be numeric.")
  if (!is.numeric(offset) || length(offset) != n) stop("`offset` must be numeric.")
  if (!is.numeric(clever_covariate) || length(clever_covariate) != n) {
    stop("`clever_covariate` must be numeric.")
  }
  if (anyNA(target) || anyNA(offset) || anyNA(clever_covariate)) {
    stop("`target`, `offset`, and `clever_covariate` cannot contain NA values.")
  }
  if (!is.null(fluctuation_weights)) {
    if (!is.numeric(fluctuation_weights) || length(fluctuation_weights) != n) {
      stop("`fluctuation_weights` must be NULL or a numeric vector of length `target`.")
    }
    if (anyNA(fluctuation_weights)) {
      stop("`fluctuation_weights` cannot contain NA values.")
    }
    if (any(fluctuation_weights < 0)) {
      stop("`fluctuation_weights` must be nonnegative.")
    }
  }

  # A degenerate clever covariate cannot identify a fluctuation direction.
  if (all(clever_covariate == 0)) {
    warning("Clever covariate is identically zero; no fluctuation applied ",
            "(eps set to 0).")
    return(list(eps = 0, intercept = 0, updated = offset, fit = NULL))
  }

  finish <- function(fit, updated_fun) {
    coefs <- stats::coef(fit)
    eps <- unname(coefs["H"])
    if (!is.finite(eps)) {
      warning("Fluctuation coefficient is not finite (NA/NaN/Inf); ",
              "no fluctuation applied (eps set to 0). This usually ",
              "indicates a degenerate clever covariate or separation.")
      eps <- 0
    }
    list(eps = eps, intercept = 0, updated = updated_fun(eps), fit = fit)
  }

  if (fluctuation_type == "identity") {
    dat <- data.frame(
      target = target,
      H = clever_covariate,
      offset = offset
    )

    fit <- stats::glm(
      target ~ -1 + H,
      family = stats::gaussian(),
      data = dat,
      offset = offset,
      weights = fluctuation_weights
    )

    return(finish(fit, function(eps) offset + eps * clever_covariate))
  }

  if (fluctuation_type == "logistic") {
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
      offset = offset_logit,
      weights = fluctuation_weights
    ))

    return(finish(fit, function(eps) {
      updated01 <- expit(offset_logit + eps * clever_covariate)
      updated01 <- clip01(updated01, clip = clip)
      from01(updated01, bounds)
    }))
  }

  stop("Unsupported fluctuation_type.")
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

.check_targeted_score <- function(ic_star, estimate, se, score_tol) {
  # After a successful TMLE update, mean(ic_star) == estimate up to
  # numerical error (each fluctuation solves its score equation exactly).
  score <- mean(ic_star) - estimate
  if (is.finite(se) && se > 0 && abs(score) > score_tol * se) {
    warning(
      "The targeted EIF estimating equation is not (numerically) solved: ",
      "mean(targeted IC) - estimate = ", formatC(score, digits = 4, format = "g"),
      " (", formatC(abs(score) / se, digits = 3, format = "f"),
      " estimated standard errors). ",
      "Check `targeting_steps` / `alpha_star` specifications and bounds."
    )
  }
  invisible(score)
}

.apply_fluctuation_map <- function(x,
                                   clever_covariate,
                                   eps,
                                   intercept = 0,
                                   fluctuation_type,
                                   bounds = NULL,
                                   clip = 1e-6) {
  if (!is.numeric(x) || !is.numeric(clever_covariate)) {
    stop("`x` and `clever_covariate` must be numeric.")
  }
  if (length(x) != length(clever_covariate)) {
    stop("`x` and `clever_covariate` must have the same length.")
  }

  if (fluctuation_type == "identity") {
    return(x + intercept + eps * clever_covariate)
  }

  if (fluctuation_type == "logistic") {
    x01 <- clip01(to01(x, bounds), clip = clip)
    x01_star <- expit(logit(x01) + intercept + eps * clever_covariate)
    x01_star <- clip01(x01_star, clip = clip)
    return(from01(x01_star, bounds))
  }

  stop("Unsupported `fluctuation_type`.")
}

.riesz_tmle_single <- function(data,
                               rc,
                               fluctuation_type,
                               bounds,
                               outcome_col,
                               fluctuation_weights,
                               clip,
                               significance_alpha,
                               score_tol) {

  y <- data[[outcome_col]]
  rc$fit(data)

  f <- rc$fit_f
  h <- rc$fit_h
  alpha <- rc$fit_alpha

  # ---- untargeted quantities used for inference ----------------------------

  # Untargeted plug-in
  psi_init <- mean(h)

  # Untargeted IC-like quantity
  ic_init <- rc$fit_ic

  # Centered EIF for variance
  eif_init <- ic_init - psi_init

  # ---- targeted procedure for the point estimate ---------------------------

  fluc_fit <- .fit_tmle_fluctuation(
    target = y,
    offset = f,
    clever_covariate = alpha,
    fluctuation_type = fluctuation_type,
    fluctuation_weights = fluctuation_weights,
    bounds = bounds,
    clip = clip
  )

  nuis_star <- .build_targeted_nuisance_list(
    rc = rc,
    data = data,
    eps = fluc_fit$eps,
    intercept = fluc_fit$intercept,
    fluctuation_type = fluctuation_type,
    bounds = bounds,
    clip = clip
  )

  f_star <- rc$eval_f(data, nuis_star)
  h_star <- rc$eval_h(data, nuis_star, f_value = f_star)

  tmle_estimate <- mean(h_star)

  # Targeted influence curve; its empirical mean should equal the TMLE
  # estimate (up to numerical error) when the targeting step succeeded.
  ic_star <- NULL
  if (!is.null(rc$ic_expr)) {
    alpha_at_star <- rc$eval_alpha(data, nuis_star)
    ic_star <- unname(rc$eval_ic(
      data, nuis_star,
      alpha = alpha_at_star, f = f_star, h = h_star
    ))
  }

  summarized_tmle_fit <- .summarize_tmle_fit(
    estimate = tmle_estimate,
    ic_for_inference = eif_init,
    significance_alpha = significance_alpha
  )

  if (!is.null(ic_star)) {
    .check_targeted_score(ic_star, tmle_estimate,
                          summarized_tmle_fit$se, score_tol)
  }

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
    ic_star = ic_star,
    eps = fluc_fit$eps,
    intercept = fluc_fit$intercept,
    f_star = f_star,
    h_star = h_star,
    fluctuation_model = fluc_fit$fit
  )
}

.riesz_tmle_composed <- function(data,
                                 rc,
                                 fluctuation_type,
                                 bounds,
                                 outcome_col,
                                 fluctuation_weights,
                                 clip,
                                 significance_alpha,
                                 score_tol) {

  y <- data[[outcome_col]]
  n <- nrow(data)

  rc$fit(data)
  J <- length(rc$rc_list)

  alpha_list <- vector("list", J)
  alpha_star_list <- vector("list", J)
  f_list <- vector("list", J)
  h_list <- vector("list", J)

  for (j in seq_len(J)) {
    rj <- rc$rc_list[[j]]
    alpha_list[[j]] <- rj$fit_alpha
    f_list[[j]] <- rj$fit_f
    h_list[[j]] <- rj$fit_h

    a_star <- rj$fit_alpha_star
    if (is.null(a_star)) {
      warning(
        "Stage ", j, " has no `alpha_star`. Falling back to the ",
        "observed-data `alpha` to update its predictions in the plug-in; ",
        "this is NOT a valid TMLE update in general, because units off ",
        "the support of the intervention (or unobserved units) receive no ",
        "update. Supply `alpha_star` on every stage (e.g. `~ 1/g` for ",
        "indicator weights, `~ 1/pi` for missingness weights)."
      )
      a_star <- alpha_list[[j]]
    }
    alpha_star_list[[j]] <- a_star
  }

  # ---- cumulative products of representers (suffix products) ---------------
  # omega_j = prod_{k = j}^{J} alpha_k       : clever covariate for the
  #                                            fluctuation of stage j's f
  # omega_star_j = alpha*_j * prod_{k>j} alpha_k : clever covariate for the
  #                                            counterfactual update of
  #                                            stage j's h (only stage j's
  #                                            own treatment is set when
  #                                            h_j is evaluated; earlier
  #                                            time points stay natural)
  omega_list <- vector("list", J)
  omega_star_list <- vector("list", J)
  running <- rep(1, n)
  for (j in rev(seq_len(J))) {
    # h_j evaluates stage j's fluctuated regression with ONLY stage j's
    # treatment intervened (the package's stagewise convention: see the
    # ComposedRieszCurve examples and Diaz et al. 2021's lmtp recursion,
    # where each pseudo-outcome shifts the current treatment only). The
    # clever covariate of that evaluation is therefore alpha*_j times the
    # NATURAL-value representers of the later stages (earlier time points),
    #   omega^d_j = alpha*_j * prod_{k > j} alpha_k,
    # not prod_{k >= j} alpha*_k. For indicator representers the two
    # coincide (1(A=a)^2 = 1(A=a) = 1(A=a) * 1), which is why binary tests
    # cannot tell them apart; for density-ratio representers (MTP stages)
    # the all-starred product breaks double robustness.
    omega_star_list[[j]] <- alpha_star_list[[j]] * running
    running <- running * alpha_list[[j]]
    omega_list[[j]] <- running
  }

  # ---- sequential fluctuations, innermost-first -----------------------------
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
      fluctuation_type = fluctuation_type,
      fluctuation_weights = fluctuation_weights,
      bounds = bounds,
      clip = clip
    )

    f_star_list[[j]] <- fluc$updated
    eps_vec[[j]] <- fluc$eps
    intercept_vec[[j]] <- fluc$intercept
    fluctuation_models[[j]] <- fluc$fit

    # Counterfactual predictions are updated with the INTERVENED clever
    # covariate (omega_star), not the observed one: h_j is the evaluation of
    # the fluctuated regression f*_j under the intervention.
    h_star_list[[j]] <- .apply_fluctuation_map(
      x = h_list[[j]],
      clever_covariate = omega_star_list[[j]],
      eps = fluc$eps,
      intercept = fluc$intercept,
      fluctuation_type = fluctuation_type,
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

  targeted_ic <- .compute_targeted_ic_composed(
    y = y,
    omega_list = omega_list,
    f_star_list = f_star_list,
    h_star_list = h_star_list
  )

  summarized_tmle_fit <- .summarize_tmle_fit(
    estimate = tmle_estimate,
    ic_for_inference = ic_for_inference,
    significance_alpha = significance_alpha
  )

  .check_targeted_score(targeted_ic, tmle_estimate,
                        summarized_tmle_fit$se, score_tol)

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
    ic = ic_for_inference - untargeted_ic$psi,
    ic_star = targeted_ic,
    f_star = f_star_list,
    h_star = h_star_list,
    fluctuation_models = fluctuation_models,
    omega = omega_list,
    targets = target_list
  )
}

.compute_targeted_ic_composed <- function(y,
                                          omega_list,
                                          f_star_list,
                                          h_star_list) {
  J <- length(omega_list)
  n <- length(y)

  ic <- rep(0, n)

  current_target <- y
  for (j in seq_len(J)) {
    ic <- ic + omega_list[[j]] * (current_target - f_star_list[[j]])
    current_target <- h_star_list[[j]]
  }

  ic <- ic + h_star_list[[J]]

  ic
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


# helper to apply intervention for fluctuating nuisance models ------------------------

.apply_intervention_values <- function(data, set = NULL) {
  data_new <- data
  if (!is.null(set)) {
    for (nm in names(set)) {
      data_new[[nm]] <- set[[nm]]
    }
  }
  data_new
}


# produce targeted nuisances ----------------------------------------------

.build_targeted_nuisance_list <- function(rc, data, eps, intercept,
                                          fluctuation_type, bounds, clip) {
  nuis_star <- rc$fit_nuis

  if (is.null(rc$targeting_steps)) {
    stop("`rc$targeting_steps` must be specified for TMLE. For each nuisance ",
         "appearing in `f` or `h`, provide `list()` (updated with the ",
         "observed-data `alpha`), `list(alpha_star = ~ ...)` (updated with ",
         "an explicit counterfactual clever covariate), or ",
         "`list(set = list(...))` (updated with `alpha` re-evaluated on ",
         "intervened data).")
  }

  for (nm in names(rc$targeting_steps)) {
    if (!nm %in% names(nuis_star)) {
      stop("`targeting_steps` refers to nuisance `", nm,
           "`, which is not in the fitted nuisance list.")
    }

    spec <- rc$targeting_steps[[nm]]

    if (!is.null(spec$alpha_star)) {
      # explicit counterfactual clever covariate for this nuisance
      if (is.numeric(spec$alpha_star)) {
        H_nm <- spec$alpha_star
      } else {
        H_nm <- rc$.eval_formula(spec$alpha_star, data, rc$fit_nuis)
      }
    } else if (!is.null(spec$set)) {
      # legacy mechanism: re-evaluate alpha on intervened data. Only valid
      # when `alpha` depends on the intervened columns directly.
      data_fluc <- .apply_intervention_values(data, spec$set)
      H_nm <- rc$eval_alpha(data_fluc, rc$fit_nuis)
    } else {
      # empty spec: update with the observed-data representer
      H_nm <- rc$fit_alpha
    }

    if (length(H_nm) == 1L) H_nm <- rep(H_nm, length(nuis_star[[nm]]))
    H_nm <- as.vector(unname(H_nm))

    nuis_star[[nm]] <- .apply_fluctuation_map(
      x = rc$fit_nuis[[nm]],
      clever_covariate = H_nm,
      eps = eps,
      intercept = intercept,
      fluctuation_type = fluctuation_type,
      bounds = bounds,
      clip = clip
    )
  }

  nuis_star
}

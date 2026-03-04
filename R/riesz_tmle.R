
riesz_tmle <- function(data,
                       rr,
                       fluctuation = 'logistic',
                       bounds,
                       outcome_col,
                       clip = 1e-6) {

  # checks
  if (!inherits(rr, "RieszRepresenter")) stop("`rr` must be a RieszRepresenter.")
  if (!is.character(outcome_col) || length(outcome_col) != 1L || outcome_col == "" || is.na(outcome_col)) {
    stop("`outcome_col` must be a non-empty character scalar.")
  }
  if (!outcome_col %in% names(data)) stop("`outcome_col` is not a column in `data`.")


  # if identity, defer to riesz_estimate
  if (fluctuation == 'identity') {
    return(riesz_estimate(data, rr))
  }

  # fluctuation on the logistic scale
  if (fluctuation == 'logistic') {

    # bounds checks
    if (is.null(bounds)) stop("For `fluctuation = 'logistic'`, you must provide `bounds = c(a, b)`.")
    if (!is.numeric(bounds) || length(bounds) != 2L) stop("`bounds` must be numeric length 2: c(a, b).")
    a <- bounds[[1]]; b <- bounds[[2]]
    if (!is.finite(a) || !is.finite(b)) stop("`bounds` must be finite.")
    if (b <= a) stop("`bounds[2]` must be > `bounds[1]`.")

    # fit rr pieces
    rr$fit(data)

    if (is.null(rr$fit_alpha)) stop("After rr$fit(data), rr$fit_alpha is NULL.")
    if (is.null(rr$fit_h)) stop("After rr$fit(data), rr$fit_h is NULL.")
    if (is.null(rr$fit_f)) stop("After rr$fit(data), rr$fit_f is NULL (need the regression f).")

    alpha <- rr$fit_alpha
    h <- rr$fit_h
    f <- rr$fit_f
    y <- data[[outcome_col]]
    n <- nrow(data)

    if (!is.numeric(alpha) || length(alpha) != n) stop("`alpha` must be numeric length nrow(data).")
    if (!is.numeric(h) || length(h) != n) stop("`h` must be numeric length nrow(data).")
    if (!is.numeric(f) || length(f) != n) stop("`f` must be numeric length nrow(data).")
    if (!is.numeric(y) || length(y) != n) stop("Outcome must be numeric length nrow(data).")

    # map to [0,1]
    y01 <- to01(y, bounds)
    f01 <- to01(f, bounds)

    # (optional) range checks; you can relax these if you prefer
    if (any(!is.finite(y01))) stop("Outcome mapped to [0,1] produced non-finite values.")
    if (any(!is.finite(f01))) stop("f mapped to [0,1] produced non-finite values.")

    # clip away from 0/1 for logit
    y01 <- clip01(y01, clip = clip)
    f01 <- clip01(f01, clip = clip)

    # Fit fluctuation: logit(f*) = logit(f) + eps * alpha
    tmle_data <- data.frame(
      y01 = y01,
      alpha = alpha,
      offset = logit(f01)
    )

    fluctuation_model <- suppressMessages(stats::glm(
      y01 ~ -1 + alpha,
      family = stats::binomial(),
      data = tmle_data,
      offset = tmle_data$offset
    )) # tends to message that non-integer #successes were given in a
    # binomial glm or that glm.fit fitted probabilities numerically 0 or 1

    eps <- as.numeric(stats::coef(fluctuation_model)[["alpha"]])
    if (!is.finite(eps)) stop("Fluctuation coefficient `eps` is non-finite.")

    f01_star <- expit(logit(f01) + eps * alpha)
    f01_star <- clip01(f01_star, clip = clip)

    # map back to [a,b]
    f_star <- from01(f01_star, bounds)

    # updated IC contribution (uncentered)
    ic_star <- h + alpha * (y - f_star)

    estimate <- mean(ic_star)
    var_estimate <- stats::var(ic_star) / n
    se <- sqrt(var_estimate)
    ci_low <- estimate + stats::qnorm(0.025) * se
    ci_high <- estimate + stats::qnorm(0.975) * se

    return(list(
      estimate = estimate,
      var = var_estimate,
      se = se,
      ci_low = ci_low,
      ci_high = ci_high,
      eps = eps,
      ic = ic_star,
      f_star = f_star
    ))
  }
}

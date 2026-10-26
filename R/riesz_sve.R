#' Stochastic-Interventional Vaccine Efficacy under Two-Phase Sampling
#' and Right Censoring
#'
#' Estimates SVE(delta) = 1 - psi_1(delta) / psi_0 (Hejazi et al. 2023),
#' where psi_1(delta) = E[E(Y | S + delta, L, C = 1, A = 1) | A = 1] and
#' psi_0 = E[E(Y | L, C = 1, A = 0) | A = 0], with the marker S measured
#' only in a phase-two sample of the vaccine arm.
#'
#' @param data data.frame with columns `A`, `S`, `R`, `C`, `Y`, covariates.
#'   `NA` values of `S` (unmeasured) and `Y` (censored) are filled with 0;
#'   those rows receive zero weight in the influence curve.
#' @param delta Numeric scalar shift of the marker (SD units if S is
#'   standardized).
#' @param nuis_censoring,nuis_phase2,nuis_placebo Functions of `data`
#'   returning the nuisance lists documented in `riesz_curve_catalog`.
#' @param nuis_shift Function of `(data, delta)` returning the shift-stage
#'   nuisance list.
#' @param estimator `"tmle"` or `"onestep"` for the two numerators.
#' @param bounds Outcome bounds for the logistic fluctuation.
#' @param significance_alpha Significance level for Wald intervals.
#' @param ... Further arguments passed to [riesz_tmle()] when
#'   `estimator = "tmle"` (e.g. `fluctuation_method`, `clip`, `score_tol`).
#' @return A list with `RieszFit`s `sve`, `risk_vaccine`, `risk_placebo`.
#' @export
riesz_sve <- function(data, delta,
                      nuis_censoring, nuis_shift, nuis_phase2, nuis_placebo,
                      estimator = c("tmle", "onestep"),
                      bounds = c(0, 1),
                      significance_alpha = 0.05, ...) {
  estimator <- match.arg(estimator)
  data <- .prepare_sve_data(data)

  rc_th1 <- riesz_curve_catalog$sve_vaccine_numerator(
    nuis_censoring = nuis_censoring,
    nuis_shift = function(d) nuis_shift(d, delta),
    nuis_phase2 = nuis_phase2
  )
  rc_th0 <- riesz_curve_catalog$sve_placebo_numerator(nuis_placebo)

  fit_num <- function(rc) {
    if (estimator == "tmle") {
      riesz_tmle(data, rc, fluctuation_type = "logistic", bounds = bounds,
                 outcome_col = "Y", significance_alpha = significance_alpha,
                 ...)
    } else {
      riesz_estimate(data, rc, significance_alpha = significance_alpha)
    }
  }
  fit_th1 <- fit_num(rc_th1)
  fit_th0 <- fit_num(rc_th0)
  fit_p1 <- riesz_estimate(data, riesz_curve_catalog$arm_proportion(1))
  fit_p0 <- riesz_estimate(data, riesz_curve_catalog$arm_proportion(0))

  # theta = (theta_1, p_1, theta_0, p_0); SVE = 1 - (t1 / p1) / (t0 / p0)
  sve <- riesz_delta(
    fit_th1, fit_p1, fit_th0, fit_p0,
    g = function(t) 1 - (t[1] * t[4]) / (t[2] * t[3]),
    grad = function(t) {
      r <- (t[1] * t[4]) / (t[2] * t[3])
      c(-r / t[1], r / t[2], r / t[3], -r / t[4])
    },
    significance_alpha = significance_alpha,
    parameter = sprintf("SVE(delta = %g)", delta),
    estimator = paste0("delta-method (", estimator, ")")
  )

  list(
    sve = sve,
    risk_vaccine = riesz_delta_ratio(fit_th1, fit_p1, significance_alpha,
                                     parameter = "psi_1(delta)"),
    risk_placebo = riesz_delta_ratio(fit_th0, fit_p0, significance_alpha,
                                     parameter = "psi_0")
  )
}

.prepare_sve_data <- function(data) {
  need <- c("A", "S", "R", "C", "Y")
  miss <- setdiff(need, names(data))
  if (length(miss)) stop("`data` is missing column(s): ",
                         paste(miss, collapse = ", "))
  for (v in c("A", "R", "C")) {
    if (anyNA(data[[v]]) || !all(data[[v]] %in% c(0, 1))) {
      stop("`", v, "` must be a 0/1 indicator without NA.")
    }
  }
  if (anyNA(data$S[data$A == 1 & data$R == 1])) {
    stop("`S` must be observed for all vaccine-arm phase-two units.")
  }
  if (anyNA(data$Y[data$C == 1])) {
    stop("`Y` must be observed for all uncensored units.")
  }
  data$S[is.na(data$S)] <- 0
  data$Y[is.na(data$Y)] <- 0
  data
}

#' Reference GLM nuisance estimators for [riesz_sve()]
#'
#' Parametric working models (logistic outcome and censoring regressions,
#' Gaussian marker density), IPW-weighted by the known phase-two sampling
#' probabilities in column `pi_col`. Swap in crossfit super learners for
#' real analyses; the contracts are documented in `riesz_curve_catalog`.
#' @param covariates Character vector of baseline covariate column names
#'   (the `L` in the nuisance contracts).
#' @param pi_col Name of the column holding the known phase-two sampling
#'   probabilities \eqn{P(R = 1 | V, A = 1)}.
#' @return A named list of nuisance functions: `censoring`, `shift`,
#'   `phase2`, `placebo`, suitable for the matching `nuis_*` arguments of
#'   [riesz_sve()].
#' @export
sve_glm_nuisances <- function(covariates, pi_col = "pi_R") {
  rhs_L <- paste(covariates, collapse = " + ")
  vac2 <- function(d) d$A == 1 & d$R == 1
  ipw <- function(d) ifelse(vac2(d), 1 / d[[pi_col]], 0)

  censoring <- function(data) {
    d <- data; d$w_ipw <- ipw(d)
    d2 <- d[vac2(d), ]
    d2c <- d[vac2(d) & d$C == 1, ]
    fitC <- stats::glm(stats::as.formula(paste("C ~ S +", rhs_L)),
                       family = stats::quasibinomial(),
                       data = d2, weights = d2$w_ipw)
    fitY <- stats::glm(stats::as.formula(paste("Y ~ S +", rhs_L)),
                       family = stats::quasibinomial(),
                       data = d2c, weights = d2c$w_ipw)
    list(mC = stats::predict(fitY, newdata = d, type = "response"),
         pC = stats::predict(fitC, newdata = d, type = "response"))
  }

  shift <- function(data, delta) {
    d <- data; d$w_ipw <- ipw(d); d$h1 <- data[["h[j=1]"]]
    d2 <- d[vac2(d), ]
    fitM <- stats::glm(stats::as.formula(paste("h1 ~ S +", rhs_L)),
                       family = stats::quasibinomial(),
                       data = d2, weights = d2$w_ipw)
    d_up <- d; d_up$S <- d$S + delta

    fitS <- stats::lm(stats::as.formula(paste("S ~", rhs_L)),
                      data = d2, weights = d2$w_ipw)
    mu <- stats::predict(fitS, newdata = d)
    w <- d2$w_ipw
    s <- sqrt(sum(w * stats::residuals(fitS)^2) / sum(w))

    list(m       = stats::predict(fitM, newdata = d,    type = "response"),
         m_shift = stats::predict(fitM, newdata = d_up, type = "response"),
         g       = stats::dnorm(d$S,         mu, s),
         g_shift = stats::dnorm(d$S - delta, mu, s),
         g_up    = stats::dnorm(d$S + delta, mu, s))
  }

  phase2 <- function(data) {
    d <- data; d$h2 <- data[["h[j=2]"]]; d$CY <- d$C * d$Y
    fit3 <- stats::glm(stats::as.formula(paste("h2 ~ C + CY +", rhs_L)),
                       family = stats::quasibinomial(), data = d[vac2(d), ])
    list(m  = ifelse(d$A == 1,
                     stats::predict(fit3, newdata = d, type = "response"), 0),
         pi = ifelse(d$A == 1, d[[pi_col]], 1))
  }

  placebo <- function(data) {
    pl <- data$A == 0
    fitC0 <- stats::glm(stats::as.formula(paste("C ~", rhs_L)),
                        family = stats::binomial(), data = data[pl, ])
    fitY0 <- stats::glm(stats::as.formula(paste("Y ~", rhs_L)),
                        family = stats::binomial(),
                        data = data[pl & data$C == 1, ])
    list(m0  = stats::predict(fitY0, newdata = data, type = "response"),
         pC0 = stats::predict(fitC0, newdata = data, type = "response"))
  }

  list(censoring = censoring, shift = shift,
       phase2 = phase2, placebo = placebo)
}

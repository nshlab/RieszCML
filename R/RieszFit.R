#' Fitted Riesz Estimator
#'
#' R6 class holding the output of [riesz_estimate()], [riesz_tmle()], and
#' [riesz_delta()]: a point estimate, influence-curve-based variance and
#' Wald interval, and (for TMLE) the targeting internals.
#'
#' @field estimate Point estimate.
#' @field var,se Estimated variance and standard error.
#' @field ci_low,ci_high Wald confidence interval limits.
#' @field significance_alpha Significance level of the interval.
#' @field estimator,parameter Labels for the estimator and target parameter.
#' @field n Sample size.
#' @field ic Influence curve used for inference.
#' @field ic_star Targeted influence curve (TMLE only).
#' @field eps,intercept Fluctuation coefficients (TMLE only).
#' @field f_star,h_star Targeted predictions (TMLE only).
#' @field fluctuation_model,fluctuation_models Fitted fluctuation GLM(s).
#' @field omega,targets Cumulative representer products and stagewise
#'   targets (composed TMLE only).
#' @export
RieszFit <- R6::R6Class(
  classname = "RieszFit",
  portable = TRUE,
  public = list(
    estimate = NULL,
    var = NULL,
    se = NULL,
    ci_low = NULL,
    ci_high = NULL,
    significance_alpha = NULL,
    estimator = NULL,
    parameter = NULL,
    n = NULL,

    ic = NULL,
    ic_star = NULL,
    eps = NULL,
    intercept = NULL,
    f_star = NULL,
    h_star = NULL,
    fluctuation_model = NULL,
    fluctuation_models = NULL,
    omega = NULL,
    targets = NULL,

    initialize = function(estimate,
                          var,
                          se,
                          ci_low,
                          ci_high,
                          significance_alpha = 0.05,
                          estimator = "one-step",
                          parameter = NULL,
                          n = NULL,
                          ic = NULL,
                          ic_star = NULL,
                          eps = NULL,
                          intercept = NULL,
                          f_star = NULL,
                          h_star = NULL,
                          fluctuation_model = NULL,
                          fluctuation_models = NULL,
                          omega = NULL,
                          targets = NULL) {
      self$estimate <- estimate
      self$var <- var
      self$se <- se
      self$ci_low <- ci_low
      self$ci_high <- ci_high
      self$significance_alpha <- significance_alpha
      self$estimator <- estimator
      self$parameter <- parameter
      self$n <- n

      self$ic <- ic
      self$ic_star <- ic_star
      self$eps <- eps
      self$intercept <- intercept
      self$f_star <- f_star
      self$h_star <- h_star
      self$fluctuation_model <- fluctuation_model
      self$fluctuation_models <- fluctuation_models
      self$omega <- omega
      self$targets <- targets

      invisible(self)
    },

    print = function(...) {
      digits <- max(3L, getOption("digits") - 2L)
      ci_level <- 100 * (1 - self$significance_alpha)

      cat("Riesz fit\n")

      if (!is.null(self$parameter) && !is.na(self$parameter)) {
        cat("Parameter: ", self$parameter, "\n", sep = "")
      }

      if (!is.null(self$estimator) && !is.na(self$estimator)) {
        cat("Estimator: ", self$estimator, "\n", sep = "")
      }

      if (!is.null(self$n) && !is.na(self$n)) {
        cat("n: ", self$n, "\n", sep = "")
      }

      cat("\n")
      cat("Estimate: ", formatC(self$estimate, digits = digits, format = "f"), "\n", sep = "")
      cat("Std. Error: ", formatC(self$se, digits = digits, format = "f"), "\n", sep = "")
      cat(formatC(ci_level, digits = 4, format = "fg"), "% CI: [",
          formatC(self$ci_low, digits = digits, format = "f"), ", ",
          formatC(self$ci_high, digits = digits, format = "f"), "]\n", sep = "")

      if (!is.null(self$eps)) {
        cat("\n")
        if (length(self$eps) == 1L) {
          cat("Fluctuation epsilon: ", formatC(self$eps, digits = digits, format = "f"), "\n", sep = "")
        } else {
          cat("Fluctuation epsilons:\n")
          eps_fmt <- formatC(self$eps, digits = digits, format = "f")
          for (j in seq_along(eps_fmt)) {
            cat("  [", j, "] ", eps_fmt[[j]], "\n", sep = "")
          }
        }
      }

      stored <- c()
      if (!is.null(self$ic)) stored <- c(stored, "ic")
      if (!is.null(self$ic_star)) stored <- c(stored, "ic_star")
      if (!is.null(self$f_star)) stored <- c(stored, "f_star")
      if (!is.null(self$h_star)) stored <- c(stored, "h_star")
      if (!is.null(self$fluctuation_model)) stored <- c(stored, "fluctuation_model")
      if (!is.null(self$fluctuation_models)) stored <- c(stored, "fluctuation_models")
      if (!is.null(self$omega)) stored <- c(stored, "omega")
      if (!is.null(self$targets)) stored <- c(stored, "targets")

      if (length(stored) > 0L) {
        cat("\nStored internals: ", paste(stored, collapse = ", "), "\n", sep = "")
      }

      invisible(self)
    },

    as_list = function() {
      list(
        estimate = self$estimate,
        var = self$var,
        se = self$se,
        ci_low = self$ci_low,
        ci_high = self$ci_high,
        significance_alpha = self$significance_alpha,
        estimator = self$estimator,
        parameter = self$parameter,
        n = self$n,
        ic = self$ic,
        ic_star = self$ic_star,
        eps = self$eps,
        intercept = self$intercept,
        f_star = self$f_star,
        h_star = self$h_star,
        fluctuation_model = self$fluctuation_model,
        fluctuation_models = self$fluctuation_models,
        omega = self$omega,
        targets = self$targets
      )
    }
  )
)

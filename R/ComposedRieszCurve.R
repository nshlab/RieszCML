#' Composed (Sequential) Riesz Curves for Nested Linear Functionals
#'
#' A `ComposedRieszCurve` represents a nested sequence of linear functionals
#' \eqn{\Psi(P) = E[h_1(\cdot\,; Q_1)]} with
#' \eqn{Q_t = E[h_{t+1}(\cdot\,; Q_{t+1}) | \bar A_t, \bar L_t]}, as in
#' Theorem 2 of the accompanying manuscript. Stages are supplied
#' *innermost-first*: `rc_list[[1]]` is the regression acting directly on the
#' outcome `Y` (the last time point), and each later stage acts on the
#' previous stage's `h`, which is made available to it as a data column named
#' `h[j=<stage>]`.
#'
#' @section Influence curve:
#' With stages indexed innermost-first, the (uncentered) efficient influence
#' curve is
#' \deqn{\phi = h_J + \sum_{j=1}^{J} \omega_j (T_j - f_j), \qquad
#'       \omega_j = \prod_{k=j}^{J} \alpha_k,}
#' where \eqn{T_1 = Y} and \eqn{T_j = h_{j-1}} for \eqn{j > 1}. Note that the
#' cumulative products \eqn{\omega_j} accumulate **from the outermost stage
#' (first time point) inward**: the residual against the innermost regression
#' (\eqn{Y - f_1}) is weighted by the *full* product of all Riesz
#' representers, and the residual against the outermost regression is
#' weighted by that stage's representer alone. This matches the EIF
#' \eqn{\sum_t \prod_{k \le t} \alpha_k [h_{t+1} - Q_t] + h_1 - \psi}
#' once the reversed stage indexing (stage 1 = final time point) is
#' accounted for.
#'
#' @export
#' @examples
#' # E[Y^{A1=1, A2=1}] in a two-timepoint longitudinal design
#' rc1 <- RieszCurve$new(  # innermost stage: time point 2
#'   nuis = function(data) {
#'     m_fit <- lm(Y ~ L1 + A1 + L2 + A2, data = data)
#'     dataA2eq1 <- data
#'     dataA2eq1$A2 <- 1
#'     g_fit <- glm(A2 ~ L1 + A1 + L2, family = binomial(), data = data)
#'
#'     list(
#'       m = predict(m_fit, newdata = data),
#'       m1 = predict(m_fit, newdata = dataA2eq1),
#'       g = predict(g_fit, newdata = data, type = "response")
#'     )
#'   },
#'   h = ~ m1,
#'   f = ~ m,
#'   alpha = ~ A2/g,
#'   alpha_star = ~ 1/g,
#'   ic_expr = ~ h + alpha * (Y - f),
#'   targeting_steps = list(m = list(), m1 = list(alpha_star = ~ 1/g))
#' )
#' rc2 <- RieszCurve$new(  # outermost stage: time point 1
#'   nuis = function(data) {
#'     m_fit <- lm(`h[j=1]` ~ L1 + A1, data = data)
#'     dataA1eq1 <- data
#'     dataA1eq1$A1 <- 1
#'     g_fit <- glm(A1 ~ L1, family = binomial(), data = data)
#'
#'     list(
#'       m = predict(m_fit, newdata = data),
#'       m1 = predict(m_fit, newdata = dataA1eq1),
#'       g = predict(g_fit, newdata = data, type = "response"))
#'   },
#'   h = ~ m1,
#'   f = ~ m,
#'   alpha = ~ A1/g,
#'   alpha_star = ~ 1/g,
#'   ic_expr = ~ h + alpha * (`h[j=1]` - f),
#'   targeting_steps = list(m = list(), m1 = list(alpha_star = ~ 1/g))
#' )
#'
#' set.seed(1)
#' sim_dgp <- function(n) {
#'   L1 <- rnorm(n)
#'   A1 <- rbinom(n, size = 1, prob = plogis(-0.2 + 1.0 * L1))
#'   L2 <- 0.5 + 0.7 * L1 + 2.0 * A1 + rnorm(n)
#'   A2 <- rbinom(n, size = 1, prob = plogis(-0.3 + 0.6 * A1 + 1.2 * L2))
#'   Y <- 1.0 + 0.4 * L1 + 1.5 * L2 + 3.0 * A2 + rnorm(n)
#'   data.frame(L1 = L1, A1 = A1, L2 = L2, A2 = A2, Y = Y)
#' }
#' # E[Y^{1,1}] = 1.0 + 1.5 * (0.5 + 2.0) + 3.0 = 7.75
#' df <- sim_dgp(n = 5000)
#'
#' # The composed influence curve, written out by hand
#' # (innermost residual carries the FULL product of representers):
#' rc1$fit(df)
#' df$`h[j=1]` <- rc1$fit_h
#' rc2$fit(df)
#' ic_manual <- rc2$fit_h +
#'   rc2$fit_alpha * (rc1$fit_h - rc2$fit_f) +
#'   rc2$fit_alpha * rc1$fit_alpha * (df$Y - rc1$fit_f)
#' mean(ic_manual)
#'
#' # Same thing via ComposedRieszCurve:
#' rc_composed <- ComposedRieszCurve$new(rc_list = list(rc1, rc2))
#' rc_composed$fit(df)
#' all.equal(rc_composed$fit_ic, ic_manual)
#' riesz_estimate(df, rc_composed)
ComposedRieszCurve <- R6::R6Class(
  classname = 'ComposedRieszCurve',
  portable = TRUE,
  public = list(
    rc_list = NULL,
    outcome_col = NULL,

    fit_ic = NULL,
    fit_omega = NULL,   # fit_omega[[j]] = prod_{k = j}^{J} alpha_k

    initialize = function(
      rc_list,
      outcome_col = "Y"
    ) {
      if (!is.list(rc_list) || length(rc_list) < 1L ||
          !all(vapply(rc_list, function(x) inherits(x, "RieszCurve"), logical(1)))) {
        stop("`rc_list` must be a non-empty list of `RieszCurve` objects.")
      }
      if (!is.character(outcome_col) || length(outcome_col) != 1L) {
        stop("`outcome_col` must be a character scalar.")
      }
      self$rc_list <- rc_list
      self$outcome_col <- outcome_col
    },

    fit = function(data) {
      J <- length(self$rc_list)
      n <- nrow(data)

      if (!self$outcome_col %in% names(data)) {
        stop("`", self$outcome_col, "` is not a column in `data`.")
      }
      y <- data[[self$outcome_col]]

      # ---- 1. fit all stages, innermost-first, chaining h columns ----
      data_modified <- data
      for (j in seq_len(J)) {
        self$rc_list[[j]]$fit(data_modified)

        for (comp in c("fit_alpha", "fit_f", "fit_h")) {
          v <- self$rc_list[[j]][[comp]]
          if (is.null(v) || length(v) != n) {
            stop("Stage ", j, ": `", comp, "` has length ",
                 length(v), " but `data` has ", n, " rows.")
          }
        }

        data_modified[[paste0('h[j=', j, ']')]] <- self$rc_list[[j]]$fit_h
      }

      # ---- 2. cumulative products of Riesz representers ----
      # omega_j = prod_{k = j}^{J} alpha_k: products accumulate from the
      # OUTERMOST stage (first time point) inward, so the innermost residual
      # (Y - f_1) is weighted by the full product. (Accumulating from stage 1
      # upward attaches the products to the wrong residuals and destroys
      # double robustness; see tests/testthat/test-composed-ordering.R.)
      omega <- vector("list", J)
      running <- rep(1, n)
      for (j in rev(seq_len(J))) {
        running <- running * self$rc_list[[j]]$fit_alpha
        omega[[j]] <- running
      }
      self$fit_omega <- omega

      # ---- 3. assemble the uncentered influence curve ----
      ic <- self$rc_list[[J]]$fit_h
      prev_target <- y
      for (j in seq_len(J)) {
        ic <- ic + omega[[j]] * (prev_target - self$rc_list[[j]]$fit_f)
        prev_target <- self$rc_list[[j]]$fit_h
      }
      self$fit_ic <- ic

      return(invisible(NULL))
    },

    print = function(...) {

      cat("ComposedRieszCurve object\n")
      cat("------------------------\n")

      J <- length(self$rc_list)

      cat("Stages: ", J, "\n", sep = "")
      cat("Status: ",
          if (!is.null(self$fit_ic)) "fitted" else "not yet fitted",
          "\n", sep = "")

      fmt_formula <- function(x) {
        if (inherits(x, "formula")) {
          paste(deparse(x), collapse = "")
        } else if (is.character(x)) {
          x
        } else if (is.numeric(x)) {
          "<numeric>"
        } else if (is.null(x)) {
          "<NULL>"
        } else {
          paste0("<", class(x)[1], ">")
        }
      }

      cat("\nPipeline (stage 1 = innermost, acting on ",
          self$outcome_col, "):\n", sep = "")

      for (j in seq_len(J)) {
        rcj <- self$rc_list[[j]]

        cat("\n")
        cat("  Stage ", j, "\n", sep = "")
        cat("  ", strrep("~", 8), "\n", sep = "")

        cat("    alpha:      ", fmt_formula(rcj$alpha), "\n", sep = "")
        cat("    alpha_star: ", fmt_formula(rcj$alpha_star), "\n", sep = "")
        cat("    f:          ", fmt_formula(rcj$f), "\n", sep = "")
        cat("    h:          ", fmt_formula(rcj$h), "\n", sep = "")

        if (!is.null(rcj$ic_expr)) {
          cat("    ic:         ", fmt_formula(rcj$ic_expr), "\n", sep = "")
        }

        if (is.function(rcj$nuis)) {
          cat("    nuis:  <function>\n")
        } else if (is.list(rcj$nuis)) {
          cat("    nuis:  <list>")
          if (!is.null(names(rcj$nuis)) && length(names(rcj$nuis)) > 0) {
            cat(" [", paste(names(rcj$nuis), collapse = ", "), "]", sep = "")
          }
          cat("\n")
        } else {
          cat("    nuis:  <", class(rcj$nuis)[1], ">\n", sep = "")
        }

        stage_status <- !is.null(rcj$fit_h) || !is.null(rcj$fit_f) || !is.null(rcj$fit_alpha)
        cat("    status: ", if (stage_status) "fit" else "not yet fit", "\n", sep = "")
      }

      cat("\nComposition rule:\n")
      cat("  stage 1 acts on ", self$outcome_col, " directly;\n", sep = "")
      cat("  each later stage acts on the previous stage's h (column `h[j=<stage>]`).\n")
      cat("  Residual weights are omega_j = prod_{k >= j} alpha_k, so the\n")
      cat("  innermost residual carries the full product of representers.\n")

      if (!is.null(self$fit_ic)) {
        cat("\nFit summary:\n")
        cat("  n: ", length(self$fit_ic), "\n", sep = "")
        cat("  mean(IC): ", formatC(mean(self$fit_ic), digits = 4, format = "f"), "\n", sep = "")
        cat("  sd(IC):   ", formatC(stats::sd(self$fit_ic), digits = 4, format = "f"), "\n", sep = "")
      }

      invisible(self)
    }

  )
)

#' @export
#' @examples
#' rc1 <- RieszCurve$new(
#'   nuis = function(data) {
#'     m <- nadir::lnr_glm(data, formula = Y ~ L1 + A1 + L2 + A2)
#'     dataA2eq1 <- data
#'     dataA2eq1$A2 <- 1
#'     g <- nadir::lnr_logistic(data, formula = A2 ~ L1 + A1 + L2)
#'
#'     return(list(
#'       m = m(data),
#'       m1 = m(dataA2eq1),
#'       g = g(data)
#'     ))
#'   },
#'   h = ~ m1,
#'   f = ~ m,
#'   alpha = ~ A2/g,
#'   ic_expr = ~ h + alpha * (Y - f)
#' )
#' rc2 <- RieszCurve$new(
#'   nuis = function(data) {
#'     m <- nadir::lnr_glm(data, formula = `h[j=1]` ~ L1 + A1)
#'     dataA1eq1 <- data
#'     dataA1eq1$A1 <- 1
#'     g <- nadir::lnr_logistic(data, formula = A1 ~ L1)
#'
#'     return(list(
#'       m = m(data),
#'       m1 = m(dataA1eq1),
#'       g = g(data)))
#'       },
#'    h = ~ m1,
#'    f = ~ m,
#'    alpha = ~ A1/g,
#'    ic_expr = ~ h + alpha * (h - f)
#' )
#'
#' df <- tibble::tibble(
#'   L1 = rnorm(n = 50, sd = 2),
#'   A1 = rbinom(n = 50, prob = plogis(L1), size = 1),
#'   L2 = rnorm(n = 50, sd = 1) + L1 + (A1*2-1) * rnorm(n = 50, mean = 1),
#'   A2 = rbinom(n = 50, prob = plogis(A1 - 0.5 + L2), size = 1),
#'   Y = L2 + A2)
#' # E[Y^\bar{1}] ==
#' rc1$fit(df)
#' df$`h[j=1]` <- rc1$fit_h
#' rc2$fit(df)
#'
#' # This is the uncentered influence curve
#' # that estimates psi=E[Y^\bar{1}] in mean:
#'
#' rc2$fit_h + rc2$fit_alpha * rc1$fit_alpha * (rc1$fit_h - rc2$fit_f) +
#'   rc1$fit_alpha * (df$Y - rc1$fit_h)
#'
#' library(tibble)
#'
#' sim_dgp <- function(n) {
#'   L1 <- rnorm(n, mean = 0, sd = 1)
#'
#'   # Baseline confounded treatment
#'   A1 <- rbinom(n, size = 1, prob = plogis(-0.2 + 1.0 * L1))
#'
#'   # Time-varying confounder affected by A1
#'   U2 <- rnorm(n, mean = 0, sd = 1)
#'   L2 <- 0.5 + 0.7 * L1 + 2.0 * A1 + U2
#'
#'   # Second treatment confounded by L2 (and depends on A1 too)
#'   A2 <- rbinom(n, size = 1, prob = plogis(-0.3 + 0.6 * A1 + 1.2 * L2))
#'
#'   # Outcome depends on L2 and A2 (so L2 is a time-varying confounder for A2 -> Y)
#'   eps <- rnorm(n, mean = 0, sd = 1)
#'   Y <- 1.0 + 0.4 * L1 + 1.5 * L2 + 3.0 * A2 + eps
#'
#'   tibble(L1 = L1, A1 = A1, L2 = L2, A2 = A2, Y = Y)
#' }
#'
#'
#' # simulation parameters
#' theta0 <- 1.0
#' theta2 <- 1.5
#' theta3 <- 3.0
#' beta0  <- 0.5
#' beta2  <- 2.0
#'
#' EY_bar1 <- theta0 + theta2 * (beta0 + beta2) + theta3
#' EY_bar1
#' df <- sim_dgp(n=10000)
#' rc1$fit(df)
#' df$`h[j=1]` <- rc1$fit_h
#' rc2$fit(df)
#'
#' ic_curve <- rc2$fit_h + rc2$fit_alpha * rc1$fit_alpha * (rc1$fit_f - rc2$fit_f) +
#'   rc1$fit_alpha * (df$Y - rc1$fit_h)
#' mean(ic_curve)
#' var(ic_curve)/nrow(df)
#'
#' # now using the ComposedRieszCurve:
#'
#' rc_list <- list(rc1, rc2)
#' rc_composed <- ComposedRieszCurve$new( rc_list = rc_list )
#' rc_composed$fit(df)
#' mean(rc_composed$fit_ic)
#' var(rc_composed$fit_ic)/nrow(df)
#' riesz_estimate(df, rc_composed)
ComposedRieszCurve <- R6::R6Class(
  classname = 'ComposedRieszCurve',
  portable = TRUE,
  public = list(
    rc_list = NULL,

    fit_ic = NULL,

    initialize = function(
      rc_list
    ) {
      self$rc_list <- rc_list
    },

    fit = function(data) {
      data_modified <- data

      product_alpha <- 1

      for (j in 1:length(self$rc_list)) {
        # fit the jth riesz representer
        self$rc_list[[j]]$fit(data_modified)

        h_str <- paste0('h[j=', j,']')

        data_modified[[h_str]] <- self$rc_list[[j]]$fit_h

        product_alpha <- product_alpha * self$rc_list[[j]]$fit_alpha

        if (j == 1) {

          self$fit_ic <-
            product_alpha * (data$Y - self$rc_list[[j]]$fit_f)
        } else {

          self$fit_ic <- self$fit_ic +
            product_alpha * (self$rc_list[[j-1]]$fit_h - self$rc_list[[j]]$fit_f)
        }

        if (j == length(self$rc_list)) {
          self$fit_ic <- self$fit_ic + self$rc_list[[j]]$fit_h
        }
      }

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

      cat("\nPipeline:\n")

      for (j in seq_len(J)) {
        rcj <- self$rc_list[[j]]

        cat("\n")
        cat("  Stage ", j, "\n", sep = "")
        cat("  ", strrep("~", 8), "\n", sep = "")

        cat("    alpha: ", fmt_formula(rcj$alpha), "\n", sep = "")
        cat("    f:     ", fmt_formula(rcj$f), "\n", sep = "")
        cat("    h:     ", fmt_formula(rcj$h), "\n", sep = "")

        if (!is.null(rcj$ic_expr)) {
          cat("    ic:    ", fmt_formula(rcj$ic_expr), "\n", sep = "")
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
      cat("  stage 1 acts on Y directly;\n")
      cat("  each later stage acts on the previous stage's h.\n")

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

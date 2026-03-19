
#' An R6 Class with Components for Influence Function Based Estimation Using Riesz representers as Weights
#'
#' The title is a bit tongue in cheek, referencing the idea of a "Rieszfluence Curve."
#'
#' These objects can be thought of as storing the necessary components to
#' construct the efficient influence function of the estimator and
#' to compute a one-step estimator or TMLE based on that.
#'
#' Each 'rc' or RieszCurve object should have a \code{nuis}, \code{alpha},
#' \code{f}, \code{h}, \code{ic_expr}, which respectively define:
#'
#'   - \code{nuis} a named list of nuisance function estimates
#'   - \code{alpha} an one-sided formula expression for the weights which can
#'   use data columns and the names of nuisance estimates
#'   - \code{f} ...
#'   - \code{h} ...
#'   - \code{ic_expr} A one-sided formula expression for an
#'     uncentered influence curve to use in IF-based estimation.
#'     Often \code{ ~ h - alpha * (Y - f) }.
#'
#' @export
#' @examples
#' rc_ate <- RieszCurve$new(
#'   nuis = function(data) {
#'
#'     m <- nadir::lnr_glm(data,
#'                         formula = Y ~ L + A)
#'     g <- nadir::lnr_logistic(
#'       data, formula = A ~ L)
#'
#'     datacf0 <- datacf1 <- data
#'     datacf0$A <- 0; datacf1$A <- 1
#'
#'     return(list(
#'       m = m(data),
#'       m1 = m(datacf1),
#'       m0 = m(datacf0),
#'       g = g(data)))
#'   },
#'   alpha = ~ A/g + (1-A)/(1-g),
#'   f = ~ m,
#'   h = ~ m1 - m0,
#'   ic_expr = ~ h + alpha * (Y - f)
#' )
#' rc_ate
#'
#' df <- tibble::tibble(
#'   L = rnorm(n = 50),
#'   A = rbinom(
#'     n = 50,
#'     size = 1,
#'     prob = plogis(L)),
#'   Y = L + rnorm(n = 50, mean = 5, sd = 1) * A)
#'
#' rc_ate$fit(data = df)
#' mean(rc_ate$fit_ic)
#' rc_ate
RieszCurve <- R6::R6Class(
  classname = 'RieszCurve',
  portable = TRUE,
  public = list(
    nuis = NULL,   # returns a named list with nuis vectors
    alpha = NULL,  # a formula expression defining weights using nuis
    f = NULL,      # a formula expression defining f using nuis
    h = NULL,      # a formula expression defining h

    fit_nuis = NULL, # for storing the fit vectors
    fit_alpha = NULL,
    fit_f = NULL,
    fit_h = NULL,
    fit_ic = NULL,

    ic_expr = NULL, # an uncentered influence curve

    initialize = function(
      nuis,
      alpha,
      f,
      h,
      ic_expr = NULL) {

      # preflight checks for nuis
      if (is.list(nuis)) {
        if ((length(nuis) > 0 & is.null(names(nuis))) | any(names(nuis) == "" | is.na(names(nuis)))) {
          stop("if nuis is a list it must be named.")
        }
      }
      if (! class(nuis) %in% c("list", "function")) {
        stop("nuis must be either a named list or function.")
      }

      self$nuis <- nuis
      self$alpha <- alpha
      self$f <- f
      self$h <- h

      if (! is.null(ic_expr)) {
        self$ic_expr <- ic_expr
      }

      return(invisible(self))
    },

    fit = function(data) {
      if (is.function(self$nuis)) {
        self$fit_nuis <- self$nuis(data)
      } else if (is.list(self$nuis)) {
        # TODO: check that every list element is a
        # numeric vector here

        self$fit_nuis <- self$nuis
      }

      # evaluating alpha
      if (is.numeric(self$alpha)) {
        self$fit_alpha <- self$alpha
      } else if (inherits(self$alpha, 'formula')) {
        alpha <- self$alpha[[2]]
        # evaluate alpha in the nuis and data environment
        alpha_eval_env <- list2env(c(as.list(data), self$fit_nuis), parent = baseenv())
        self$fit_alpha <- eval(alpha, envir = alpha_eval_env)
      }

      # evaluating f
      if (is.numeric(self$f)) {
        self$fit_f <- self$f
      } else if (inherits(self$f, 'formula')) {
        f <- self$f[[2]]
        # evaluate f in the nuis and data environment
        f_eval_env <- list2env(c(as.list(data), self$fit_nuis), parent = baseenv())
        self$fit_f <- eval(f, envir = f_eval_env)
      }

      # evaluating h
       if (inherits(self$h, 'formula')) {
        h <- self$h[[2]]
        # evaluate h in the nuis and data environment
        h_eval_env <- list2env(c(as.list(data), self$fit_nuis, list(f = self$fit_f)), parent = baseenv())
        self$fit_h <- eval(h, envir = h_eval_env)
      } else if (is.numeric(self$h)) {
        self$fit_h <- self$h
      }

      if (inherits(self$ic_expr, 'formula')) {
        ic <- self$ic_expr[[2]]

        ic_eval_env <- list2env(
          c(as.list(data), self$fit_nuis,
            list(alpha = self$fit_alpha,
                 f = self$fit_f,
                 h = self$fit_h)),
          parent = baseenv())

        fit_ic <- eval(ic, envir = ic_eval_env)
        fit_ic <- unname(fit_ic)
        self$fit_ic <- fit_ic
      }

      return(invisible(NULL))
    },

    print = function(...) {

      cat("RieszCurve object\n")
      cat("-----------------\n")

      # Helper to format formulas nicely
      fmt_formula <- function(x) {
        if (inherits(x, "formula")) {
          paste(deparse(x), collapse = "")
        } else if (is.character(x)) {
          x
        } else {
          paste0("<", class(x)[1], ">")
        }
      }

      cat("Components:\n")

      cat("  alpha:  ", fmt_formula(self$alpha), "\n", sep = "")
      cat("  f:      ", fmt_formula(self$f), "\n", sep = "")
      cat("  h:      ", fmt_formula(self$h), "\n", sep = "")
      cat("  ic:     ", fmt_formula(self$ic_expr), "\n", sep = "")

      # Nuisance info (don’t print function body)
      cat("\nNuisance specification: \n")
      if (is.function(self$nuis)) {
        cat("  nuis: <function>\n")
      } else if (is.list(self$nuis)) {
        cat("  nuis: <list of precomputed vectors>\n")
        cat("    names: ", paste(names(self$nuis), collapse = ", "), "\n", sep = "")
      } else {
        cat("  nuis: <", class(self$nuis)[1], ">\n", sep = "")
      }

      # If fitted, show summary of IC
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





#' Components for Influence-Function-Based Estimation Using Riesz Representers as Weights
#'
#' A `RieszCurve` stores the components needed to construct the (uncentered)
#' efficient influence function (EIF) of a linear functional
#' \eqn{\Psi(\eta) = E[h(O; \eta)]}, and to compute a one-step estimator
#' ([riesz_estimate()]) or a TMLE ([riesz_tmle()]) based on it.
#'
#' @details
#' Each `RieszCurve` is defined by:
#' \describe{
#'   \item{`nuis`}{Either a function of `data` returning a *named list* of
#'     numeric nuisance vectors (each of length `nrow(data)`, or scalar), or
#'     such a named list directly (for precomputed nuisances).}
#'   \item{`alpha`}{A one-sided formula (or numeric vector) giving the Riesz
#'     representer evaluated at the *observed* data, e.g. `~ I(A == 1)/g`.
#'     May reference data columns and names of nuisance estimates.}
#'   \item{`f`}{A one-sided formula (or numeric vector) giving the nuisance
#'     regression evaluated at the *observed* data (the "offset" that the
#'     residual in the EIF is taken against), e.g. `~ m` for
#'     \eqn{m(A, L) = E[Y | A, L]}.}
#'   \item{`h`}{A one-sided formula (or numeric vector) giving the
#'     transformation \eqn{h(O; \eta)} whose empirical mean is the plug-in
#'     estimator, typically the regression evaluated *under the intervention*,
#'     e.g. `~ m1` for \eqn{m(1, L)}.}
#'   \item{`alpha_star`}{Optional one-sided formula (or numeric vector) giving
#'     the Riesz representer evaluated *under the intervention defining `h`*,
#'     e.g. `~ 1/g` when `alpha = ~ I(A == 1)/g`. This is the clever covariate
#'     with which counterfactual predictions (`h`, and targeted nuisances that
#'     enter `h`) must be updated during TMLE. It is required for valid TMLE
#'     updates in [riesz_tmle()] whenever `h` differs from `f`, and for
#'     composed (sequential) TMLE. For stages where `h` is an observed-data
#'     quantity equal to `f`, it may be omitted.}
#'   \item{`ic_expr`}{A one-sided formula for the *uncentered* influence
#'     curve, most often `~ h + alpha * (Y - f)`. Symbols `alpha`, `f`, `h`
#'     refer to the evaluated components above.}
#'   \item{`targeting_steps`}{A named list describing how each nuisance is
#'     fluctuated during TMLE. Each element must itself be a list containing
#'     at most one of:
#'     \itemize{
#'       \item `set`: a named list of data-column interventions (e.g.
#'         `list(A = 1)`); the clever covariate for that nuisance is `alpha`
#'         re-evaluated on the intervened data. Only appropriate when `alpha`
#'         depends on the intervened columns directly (e.g. indicator
#'         weights), *not* when `alpha` is built purely from precomputed
#'         nuisance vectors (e.g. density ratios for shift interventions).
#'       \item `alpha_star`: a one-sided formula giving the clever covariate
#'         for that nuisance directly (e.g. `~ 1/g`), evaluated in the
#'         (data, nuisances) environment. This is the recommended, general
#'         mechanism.
#'     }
#'     An empty list (`list()`) means the nuisance is updated with the
#'     observed-data `alpha` (appropriate for the natural-value regression
#'     `m` appearing in `f`).}
#' }
#'
#' @export
#' @examples
#' rc_ate <- RieszCurve$new(
#'   nuis = function(data) {
#'
#'     m_fit <- lm(Y ~ L + A, data = data)
#'     g_fit <- glm(A ~ L, family = binomial(), data = data)
#'
#'     datacf0 <- datacf1 <- data
#'     datacf0$A <- 0; datacf1$A <- 1
#'
#'     return(list(
#'       m = predict(m_fit, newdata = data),
#'       m1 = predict(m_fit, newdata = datacf1),
#'       m0 = predict(m_fit, newdata = datacf0),
#'       g = predict(g_fit, newdata = data, type = "response")))
#'   },
#'   alpha = ~ A/g - (1 - A)/(1 - g),
#'   f = ~ m,
#'   h = ~ m1 - m0,
#'   ic_expr = ~ h + alpha * (Y - f),
#'   targeting_steps = list(
#'     m  = list(),
#'     m1 = list(alpha_star = ~ 1/g),
#'     m0 = list(alpha_star = ~ -1/(1 - g))
#'   )
#' )
#' rc_ate
#'
#' set.seed(1)
#' df <- data.frame(L = rnorm(n = 200))
#' df$A <- rbinom(n = 200, size = 1, prob = plogis(df$L))
#' df$Y <- df$L + rnorm(n = 200, mean = 5, sd = 1) * df$A
#'
#' rc_ate$fit(data = df)
#' mean(rc_ate$fit_ic)
#' rc_ate
RieszCurve <- R6::R6Class(
  classname = 'RieszCurve',
  portable = TRUE,
  public = list(
    nuis = NULL,        # returns a named list with nuisance vectors
    alpha = NULL,       # formula/numeric: Riesz representer at observed data
    alpha_star = NULL,  # formula/numeric: Riesz representer under the intervention defining h
    f = NULL,           # formula/numeric: regression at observed data
    h = NULL,           # formula/numeric: plug-in transformation (counterfactual evaluation)
    targeting_steps = NULL, # named list specifying TMLE updates per nuisance

    fit_nuis = NULL,    # fitted nuisance list
    fit_alpha = NULL,
    fit_alpha_star = NULL,
    fit_f = NULL,
    fit_h = NULL,
    fit_ic = NULL,

    ic_expr = NULL,     # uncentered influence curve

    initialize = function(
      nuis,
      alpha,
      f,
      h,
      ic_expr = NULL,
      alpha_star = NULL,
      targeting_steps = NULL) {

      # ---- preflight checks: nuis ----
      if (is.list(nuis)) {
        if ((length(nuis) > 0 & is.null(names(nuis))) ||
            any(names(nuis) == "" | is.na(names(nuis)))) {
          stop("if `nuis` is a list it must be fully named.")
        }
      }
      if (!class(nuis)[1] %in% c("list", "function")) {
        stop("`nuis` must be either a named list or a function of `data`.")
      }

      # ---- preflight checks: formula-or-numeric components ----
      .check_form <- function(x, nm, allow_null = FALSE) {
        if (is.null(x)) {
          if (allow_null) return(invisible(NULL))
          stop("`", nm, "` must be supplied.")
        }
        if (!(inherits(x, "formula") || is.numeric(x))) {
          stop("`", nm, "` must be a one-sided formula or a numeric vector.")
        }
        if (inherits(x, "formula") && length(x) != 2L) {
          stop("`", nm, "` must be a one-sided formula, e.g. `~ m1`.")
        }
        invisible(NULL)
      }
      .check_form(alpha, "alpha")
      .check_form(f, "f")
      .check_form(h, "h")
      .check_form(alpha_star, "alpha_star", allow_null = TRUE)
      .check_form(ic_expr, "ic_expr", allow_null = TRUE)

      # ---- preflight checks: targeting_steps ----
      if (!is.null(targeting_steps)) {
        if (!is.list(targeting_steps)) {
          stop("`targeting_steps` must be a (fully named) list.")
        }
        if (length(targeting_steps) > 0 &&
            (is.null(names(targeting_steps)) || any(names(targeting_steps) == ""))) {
          stop("`targeting_steps` must be a fully named list.")
        }
        for (nm in names(targeting_steps)) {
          spec <- targeting_steps[[nm]]
          if (!is.list(spec)) {
            stop("`targeting_steps[['", nm, "']]` must be a list ",
                 "(possibly empty), e.g. `list(set = list(A = 1))` or ",
                 "`list(alpha_star = ~ 1/g)`.")
          }
          bad <- setdiff(names(spec), c("set", "alpha_star"))
          if (length(bad) > 0) {
            stop("Unknown field(s) in `targeting_steps[['", nm, "']]`: ",
                 paste(bad, collapse = ", "),
                 ". Allowed fields are `set` and `alpha_star`. ",
                 "Did you write `list(A = 1)` instead of `list(set = list(A = 1))`?")
          }
          if (!is.null(spec$set) && !is.null(spec$alpha_star)) {
            stop("`targeting_steps[['", nm, "']]` may specify `set` or ",
                 "`alpha_star`, not both.")
          }
          if (!is.null(spec$set) &&
              (!is.list(spec$set) || is.null(names(spec$set)))) {
            stop("`targeting_steps[['", nm, "']]$set` must be a named list ",
                 "of data-column interventions, e.g. `list(A = 1)`.")
          }
          if (!is.null(spec$alpha_star) &&
              !(inherits(spec$alpha_star, "formula") || is.numeric(spec$alpha_star))) {
            stop("`targeting_steps[['", nm, "']]$alpha_star` must be a ",
                 "one-sided formula or a numeric vector.")
          }
        }
      }

      self$nuis <- nuis
      self$alpha <- alpha
      self$alpha_star <- alpha_star
      self$f <- f
      self$h <- h
      self$ic_expr <- ic_expr
      self$targeting_steps <- targeting_steps

      return(invisible(self))
    },

    fit = function(data) {
      n <- nrow(data)

      if (is.function(self$nuis)) {
        self$fit_nuis <- self$nuis(data)
      } else if (is.list(self$nuis)) {
        self$fit_nuis <- self$nuis
      }

      # validate fitted nuisances: named numeric vectors of length n or 1
      if (is.list(self$fit_nuis) && length(self$fit_nuis) > 0) {
        if (is.null(names(self$fit_nuis)) || any(names(self$fit_nuis) == "")) {
          stop("The fitted nuisance list must be fully named.")
        }
        for (nm in names(self$fit_nuis)) {
          v <- self$fit_nuis[[nm]]
          if (!is.numeric(v)) {
            stop("Nuisance `", nm, "` is not numeric.")
          }
          if (!(length(v) %in% c(1L, n))) {
            stop("Nuisance `", nm, "` has length ", length(v),
                 " but `data` has ", n, " rows. Nuisance functions must ",
                 "return predictions for every row of `data` (do not ",
                 "drop rows, e.g. via `drop_na()`, inside `nuis`).")
          }
          if (anyNA(v)) {
            stop("Nuisance `", nm, "` contains NA values.")
          }
        }
      }

      recycle_n <- function(x, what) {
        if (is.null(x)) return(NULL)
        if (length(x) == 1L) x <- rep(x, n)
        if (length(x) != n) {
          stop("Evaluated `", what, "` has length ", length(x),
               " but `data` has ", n, " rows.")
        }
        # plain numeric vector: drop names and classes such as `AsIs`
        # picked up from I(...) terms in formulas
        as.vector(unname(x))
      }

      eval_component <- function(component) {
        if (is.null(component)) return(NULL)
        if (is.numeric(component)) return(component)
        expr <- component[[2]]
        env <- list2env(c(as.list(data), self$fit_nuis), parent = baseenv())
        eval(expr, envir = env)
      }

      self$fit_alpha <- recycle_n(eval_component(self$alpha), "alpha")
      self$fit_alpha_star <- recycle_n(eval_component(self$alpha_star), "alpha_star")
      self$fit_f <- recycle_n(eval_component(self$f), "f")

      # h may additionally reference the evaluated f
      if (inherits(self$h, 'formula')) {
        h_expr <- self$h[[2]]
        h_eval_env <- list2env(
          c(as.list(data), self$fit_nuis, list(f = self$fit_f)),
          parent = baseenv())
        self$fit_h <- recycle_n(eval(h_expr, envir = h_eval_env), "h")
      } else if (is.numeric(self$h)) {
        self$fit_h <- recycle_n(self$h, "h")
      }

      if (inherits(self$ic_expr, 'formula')) {
        ic <- self$ic_expr[[2]]
        ic_eval_env <- list2env(
          c(as.list(data), self$fit_nuis,
            list(alpha = self$fit_alpha,
                 f = self$fit_f,
                 h = self$fit_h)),
          parent = baseenv())
        self$fit_ic <- recycle_n(eval(ic, envir = ic_eval_env), "ic_expr")
      }

      return(invisible(NULL))
    },

    .eval_formula = function(formula, data, nuis_list, extra = list()) {
      expr <- formula[[2]]
      env <- list2env(c(as.list(data), nuis_list, extra), parent = baseenv())
      eval(expr, envir = env)
    },

    eval_alpha = function(data, nuis_list) {
      if (is.numeric(self$alpha)) return(self$alpha)
      self$.eval_formula(self$alpha, data, nuis_list)
    },

    eval_alpha_star = function(data, nuis_list) {
      if (is.null(self$alpha_star)) return(NULL)
      if (is.numeric(self$alpha_star)) return(self$alpha_star)
      self$.eval_formula(self$alpha_star, data, nuis_list)
    },

    eval_f = function(data, nuis_list) {
      if (is.numeric(self$f)) return(self$f)
      self$.eval_formula(self$f, data, nuis_list)
    },

    eval_h = function(data, nuis_list, f_value = NULL) {
      if (is.numeric(self$h)) return(self$h)
      self$.eval_formula(self$h, data, nuis_list, extra = list(f = f_value))
    },

    eval_ic = function(data, nuis_list, alpha, f, h) {
      self$.eval_formula(
        self$ic_expr,
        data,
        nuis_list,
        extra = list(alpha = alpha, f = f, h = h)
      )
    },

    print = function(...) {

      cat("RieszCurve object\n")
      cat("-----------------\n")

      fmt_formula <- function(x) {
        if (inherits(x, "formula")) {
          paste(deparse(x), collapse = "")
        } else if (is.character(x)) {
          x
        } else if (is.null(x)) {
          "<not supplied>"
        } else {
          paste0("<", class(x)[1], ">")
        }
      }

      cat("Components:\n")

      cat("  alpha:      ", fmt_formula(self$alpha), "\n", sep = "")
      cat("  alpha_star: ", fmt_formula(self$alpha_star), "\n", sep = "")
      cat("  f:          ", fmt_formula(self$f), "\n", sep = "")
      cat("  h:          ", fmt_formula(self$h), "\n", sep = "")
      cat("  ic:         ", fmt_formula(self$ic_expr), "\n", sep = "")

      cat("\nNuisance specification: \n")
      if (is.function(self$nuis)) {
        cat("  nuis: <function>\n")
      } else if (is.list(self$nuis)) {
        cat("  nuis: <list of precomputed vectors>\n")
        cat("    names: ", paste(names(self$nuis), collapse = ", "), "\n", sep = "")
      } else {
        cat("  nuis: <", class(self$nuis)[1], ">\n", sep = "")
      }

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

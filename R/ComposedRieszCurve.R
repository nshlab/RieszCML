
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
#'
#' riesz_estimate(df, rc_composed)
#'
#' # What did we learn:
#' # riesz representers have to be fit in calculation order
#' # 1st is a regression on Y but all after are on h_{j-1}
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
      }

  )
)


#' @examples
#' rr1 <- RieszRepresenter$new(
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
#' rr2 <- RieszRepresenter$new(
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
#' rr1$fit(df)
#' df$`h[j=1]` <- rr1$fit_h
#' rr2$fit(df)
#'
#' # This is the uncentered influence curve
#' # that estimates psi=E[Y^\bar{1}] in mean:
#'
#' rr2$fit_h + rr2$fit_alpha * rr1$fit_alpha * (rr1$fit_h - rr2$fit_f) +
#'   rr1$fit_alpha * (df$Y - rr1$fit_h)
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
#' rr1$fit(df)
#' df$`h[j=1]` <- rr1$fit_h
#' rr2$fit(df)
#'
#' ic_curve <- rr2$fit_h + rr2$fit_alpha * rr1$fit_alpha * (rr1$fit_f - rr2$fit_f) +
#'   rr1$fit_alpha * (df$Y - rr1$fit_h)
#' mean(ic_curve)
#' var(ic_curve)/nrow(df)
#'
#' # now using the ComposedRieszRepresenter:
#'
#' rr_list <- list(rr1, rr2)
#' rr_composed <- ComposedRieszRepresenter$new( rr_list = rr_list )
#' rr_composed$fit(df)
#' mean(rr_composed$fit_ic)
#' var(rr_composed$fit_ic)/nrow(df)
#'
#' riesz_estimate(df, rr_composed)
#'
#' What did we learn:
#' # riesz representers have to be fit in calculation order
#' # 1st is a regression on Y but all after are on h_{j-1}
ComposedRieszRepresenter <- R6::R6Class(
  classname = 'ComposedRieszRepresenter',
  portable = TRUE,
  public = list(
    rr_list = NULL,

    fit_ic = NULL,

    initialize = function(
      rr_list
    ) {
      self$rr_list <- rr_list
    },

    fit = function(data) {
      data_modified <- data

      product_alpha <- 1

      for (i in 1:length(self$rr_list)) {
        # fit the ith riesz representer
        self$rr_list[[i]]$fit(data_modified)

        h_str <- paste0('h[j=', i,']')

        data_modified[[h_str]] <- self$rr_list[[i]]$fit_h

        product_alpha <- product_alpha * self$rr_list[[i]]$fit_alpha

        if (i == 1) {

          self$fit_ic <-
            product_alpha * (data$Y - self$rr_list[[i]]$fit_f)
        } else {

          self$fit_ic <- self$fit_ic +
            product_alpha * (self$rr_list[[i-1]]$fit_h - self$rr_list[[i]]$fit_f)
        }

        if (i == length(self$rr_list)) {
          self$fit_ic <- self$fit_ic + self$rr_list[[i]]$fit_h
        }
      }

      return(invisible(NULL))
      }

  )
)

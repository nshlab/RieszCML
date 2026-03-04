
riesz_representer_catalog <- list()

riesz_representer_catalog$cfmean_a <- function(
    L,
    A,
    a,
    Y,
    nuis) {

  alpha_formula <- as.formula(paste0('~ I(', A, '==', a, ')/g'))

  rr <- RieszRepresenter$new(
    nuis = nuis,
    alpha = alpha_formula,
    f = ~ m,
    h = ~ ma,
    ic_expr = ~ h + alpha * (Y - f)
  )
}

# # example
# rr_cfmean_1 <- riesz_representer_catalog$cfmean_a(
#   L = 'L',
#   A = 'A',
#   a = 1,
#   Y = 'Y',
#
#   nuis = function(data) {
#     m <- nadir::lnr_glm(data, Y ~ L + A)
#     dataAeq_a <- data
#     dataAeq_a$A <- 1
#
#     g <- nadir::lnr_logistic(data, I(A == 1) ~ L)
#     return(list(
#       m = m(data),
#       ma = m(dataAeq_a),
#       g = g(data)))
#   })
#
# rr_cfmean_0 <- riesz_representer_catalog$cfmean_a(
#   L = 'L',
#   A = 'A',
#   a = 0,
#   Y = 'Y',
#
#   nuis = function(data) {
#     m <- nadir::lnr_glm(data, Y ~ L + A)
#     dataAeq_a <- data
#     dataAeq_a$A <- 0
#
#     g <- nadir::lnr_logistic(data, I(A == 0) ~ L)
#     return(list(
#       m = m(data),
#       ma = m(dataAeq_a),
#       g = g(data)))
#   }
# )
#
# n <- 10000
# df <- tibble::tibble(
#   L = rnorm(n = n),
#   A = rbinom(n = n, size = 1, prob = plogis(L)),
#   Y = L + A + rnorm(n = n, sd = .1))
#
# riesz_estimate(data = df, rr = rr_cfmean_0)
# riesz_estimate(data = df, rr = rr_cfmean_1)


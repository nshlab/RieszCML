
# # Old / Incomplete
# composed_riesz_tmle <- function(data,
#                                 rc,
#                                 fluctuation = 'logistic',
#                                 bounds,
#                                 outcome_col,
#                                 use_intercept_and_weights = TRUE,
#                                 clip = 1e-6,
#                                 significance_alpha = 0.05) {
#   # checks
#   # ... to fill in later
#
#   # if identity, defer to riesz_estimate
#   if (fluctuation == 'identity') {
#     return(riesz_estimate(data, rc))
#   }
#
#   if (fluctuation == 'logistic') {
#
#     # bounds check
#
#     # fit rc pieces
#     rc$fit(data)
#
#     omega <- matrix(data = NA,
#                     nrow = nrow(data),
#                     ncol = length(rc$rc_list))
#
#     for (j in 1:length(rc$rc_list)) {
#       omega[, j] <- apply(sapply(1:j, \(.x) {
#         rc$rc_list[[.x]]$fit_alpha
#       }), 1, prod)
#     }
#
#     # checks for fit nuisances
#     for (j in length(rc$rc_list):1) {
#       tmle_data <- data.frame(
#         `h_j_next` = if (j == length(rc$rc_list)) data$Y else rc$rc_list[[j+1]]$fit_h,
#         f_j = rc$rc_list[[j]]$fit_f,
#         omega_j = omega[,j]
#       )
#
#       tmle_formula <- if (use_intercept_and_weights) {
#         h_j_next ~ 1
#       } else {
#         h_j_next ~ -1 + omega_j
#       }
#       weights <- if (use_intercept_and_weights) logit(tmle_data$omega_j) else NULL
#
#       fit <- glm(
#         formula = tmle_formula,
#         data = tmle_data,
#         weights = weights,
#         family = binomial(link='logit'),
#         offset = logit(tmle_data$f_j)
#       )
#
#       eps <- if (use_intercept_and_weights) {
#         coef(fit)[['(Intercept)']]
#       } else {
#         coef(fit)[['omega_j']]
#       }
#
#
#     }
#   }
# }

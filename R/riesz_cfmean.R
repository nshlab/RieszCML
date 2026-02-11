riesz_cfmean <- function(a = 1) {

  rr <- RieszRepresenter$new(
    name = 'cfmean_a1',
    description = sprintf("Counterfactual mean E[E(Y|A=%s,L)] with Riesz EIF", a),
    required_nuisances = c('Q_a'),
    psi_hat = function(data, nuis) {
      mean(nuis$Q_a)
    },
    eif = function(data, nuis, psi = NULL) {

    }
  )
}

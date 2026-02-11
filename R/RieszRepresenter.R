
#' @importFrom assertthat assert_that
RieszRepresenter <- R6::R6Class(
  "RieszRepresenter",
  lock_objects = FALSE,
  public = list(
    name = NULL,
    description = NULL,

    required_nuisances = NULL,
    psi_hat = NULL,
    eif = NULL,
    eval = NULL,
    print = NULL,

    initialize = function(name,
                          description = '',
                          required_nuisances,
                          psi_hat,
                          eif,
                          eval,
                          print) {
      self$name <- name
      self$description <- description

      if (! missing(required_nuisances)) {
        # we expect a character vector here
        self$required_nuisances <- required_nuisances
      }

      # accept psi_hat function argument
      if (! missing(psi_hat)) {
        assert_that(is.function(psi_hat), msg = "psi_hat(data, nuis) must be a function.")
        assert_that(
          all(c('data', 'nuis') %in% formalArgs(psi_hat)),
          msg = 'psi_hat must be a function taking data and nuis as arguments')

        self$psi_hat <- psi_hat
        environment(self$psi_hat) <- environment(self$initialize)
      }

      # accept eif function argument
      if (! missing(eif)) {
        assert_that(is.function(eif), msg = "eif(data, nuis, psi) must be a function." )
        assert_that(
          all(c('data', 'nuis', 'psi') %in% formalArgs(eif)),
          msg = 'eif must be a function taking data, nuis, and psi as arguments.')

        self$eif <- eif
        environment(self$eif) <- environment(self$initialize)
      }

      # accept eval function argument
      if (! missing(eval)) {
        assert_that(is.function(eval), msg = "eval(data, nuis) must be a function." )
        assert_that(
          all(c('data', 'nuis') %in% formalArgs(eval)),
          msg = 'eval must be a function taking data and nuis as arguments.')

        self$eval <- eval
        environment(self$eval) <- environment(self$initialize)
      }

      # accept print function argument
      if (! missing(print)) {
        assert_that(is.function(print), msg = "print must be a function." )

        self$print <- print
        environment(self$print) <- environment(self$initialize)
      }

    }

  )
)

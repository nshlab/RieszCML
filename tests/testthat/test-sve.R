sim_sve <- function(n) {
  L <- rnorm(n)
  A <- rbinom(n, 1, 0.5)
  S <- ifelse(A == 1, rnorm(n, 0.5 * L, 1), 0)
  pY <- ifelse(A == 1, plogis(-3 - 0.8 * S + 0.3 * L), plogis(-2 + 0.3 * L))
  Y <- rbinom(n, 1, pY)
  C <- rbinom(n, 1, plogis(2 + 0.3 * L - 0.3 * A * S))
  pi_R <- ifelse(A == 1, ifelse(C == 1 & Y == 1, 1, 0.3), 1)
  R <- ifelse(A == 1, rbinom(n, 1, pi_R), 0)
  data.frame(L = L, A = A, R = R, C = C, pi_R = pi_R,
             S = ifelse(R == 1, S, NA), Y = ifelse(C == 1, Y, NA))
}

sve_truth <- function(delta, n_mc = 2e6) {
  L <- rnorm(n_mc); S <- rnorm(n_mc, 0.5 * L, 1)
  1 - mean(plogis(-3 - 0.8 * (S + delta) + 0.3 * L)) /
      mean(plogis(-2 + 0.3 * L))
}

test_that("riesz_sve recovers the truth (one-step and TMLE)", {
  set.seed(1)
  df <- sim_sve(20000)
  nu <- sve_glm_nuisances(covariates = "L")
  for (delta in c(-0.5, 0, 0.5)) {
    truth <- sve_truth(delta)
    for (est in c("onestep", "tmle")) {
      fit <- riesz_sve(df, delta,
                       nuis_censoring = nu$censoring, nuis_shift = nu$shift,
                       nuis_phase2 = nu$phase2, nuis_placebo = nu$placebo,
                       estimator = est)$sve
      expect_lt(abs(fit$estimate - truth), 4 * fit$se)
    }
  }
})

test_that("SVE(0) agrees with the ordinary censoring-adjusted VE", {
  set.seed(2)
  df <- sim_sve(20000)
  nu <- sve_glm_nuisances("L")
  fit <- riesz_sve(df, 0, nu$censoring, nu$shift, nu$phase2, nu$placebo,
                   estimator = "onestep")$sve
  expect_lt(abs(fit$estimate - sve_truth(0)), 4 * fit$se)
})

#' @title Simulate from the model.
#' @export
#' @description Draw a set of parameters from the prior of the model,
#'   then draw a set of data from the likelihood conditional
#'   on those parameter values.
#' @details The model matrix `x` is also random.
#' @return A Stan-compatible list of data and hyperparameters.
#'   The `.join_data` element is a list of parameter values that
#'   `stantargets` uses to compare the prior to the posterior:
#'   <https://docs.ropensci.org/stantargets/articles/simulation.html>
#' @param n_patient Positive integer of length 1, number of patients.
#' @param n_visit Positive integer of length 1,
#'   number of scheduled study visits.
#' @param n_beta Positive integer of length 1, number of fixed effects.
#' @param n_missing Positive integer of length 1,
#'   number of missing observations.
#' @param s_beta Positive numeric of length 1, prior standard deviation
#'   of the fixed effect parameters beta.
#' @param s_sigma Positive numeric of length 1, prior uniform upper bound
#'   of the residual standard deviation parameters sigma.
#' @param s_lambda Positive numeric of length 1, prior LKJ shape parameter
#'   of the Cholesky factor lambda of the within-patient residual
#'   correlation matrix.
#' @examples
#' data <- simulate_data()
#' str(data)
#' model <- cmdstanr::cmdstan_model("model.stan")
simulate_data <- function(
  n_patient = 50,
  n_visit = 4,
  n_beta = 2,
  n_missing = as.integer(0.1 * n_patient * n_visit),
  s_beta = 1,
  s_sigma = 1,
  s_lambda = 1
) {
  x <- matrix(rnorm(n = n_patient * n_visit * n_beta), nrow = n_beta)
  beta <- rnorm(n = n_beta, mean = 0, sd = s_beta)
  sigma <- runif(n = n_visit, min = 0, max = s_sigma)
  correlation <- trialr::rlkjcorr(n = 1, K = n_visit, eta = s_lambda)
  lambda <- t(chol(correlation))
  covariance <- diag(sigma) %*% correlation %*% diag(sigma)
  mu <- matrix(t(x) %*% beta, nrow = n_visit)
  y <- apply(X = mu, MARGIN = 2, FUN = function(mu_column) {
    MASS::mvrnorm(n = 1, mu = mu_column, Sigma = covariance)
  })
  missing <- matrix(0L, nrow = nrow(y), ncol = ncol(y))
  index <- sample.int(n = prod(dim(y)), size = n_missing)
  y[index] <- -1e6
  missing[index] <- 1L
  list(
    n_patient = n_patient,
    n_visit = n_visit,
    n_beta = n_beta,
    n_missing = n_missing,
    missing = missing,
    s_beta = s_beta,
    s_sigma = s_sigma,
    s_lambda = s_lambda,
    x = x,
    y = y,
    .join_data = list(
      beta = beta,
      sigma = sigma,
      lambda = lambda
    )
  )
}

#' @title Simulate from the model.
#' @export
#' @description Draw a set of parameters from the prior of the model,
#'   then draw a set of data from the likelihood conditional
#'   on those parameter values.
#' @details The model matrix `x` is also random.
#' @return A Stan-compatible list of data and hyperparameters.
#'   The model matrix x is given as the transpose of the model matrix
#'   in the slides in order to speed up matrix multiplication in the Stan code.
#'   The data records (columns) in x are ordered
#'   by visit within patient, i.e. same as the row order of
#'   tidyr::expand_grid(patient = seq_len(n_patient), visit = seq_len(n_visit))
#'   The `.join_data` element is a list of parameter values that
#'   `stantargets` uses to compare the prior to the posterior:
#'   <https://docs.ropensci.org/stantargets/articles/simulation.html>.
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
#' data <- simulate_data(n_missing = 0)
#' str(data)
#' data$.join_data <- NULL
#' model <- cmdstanr::cmdstan_model("model.stan")
#' fit <- model$sample(data = data, iter_sampling = 2e3, iter_warmup = 2e3)
#' fit$summary()
simulate_data <- function(
  n_patient = 50,
  n_visit = 4,
  n_beta = 2,
  n_missing = 0,
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

#' @title Get simulation-based calibration rank statistics.
#' @export
#' @description Run `SBC::calculate_ranks_draws_matrix()` to calculate
#'   simulation-based calibration rank statistics on each model parameter.
#' @references
#'   Kim, Shinyoung, Hyunji Moon, Martin Modrák, and Teemu Säilynoja. 2022.
#'     SBC: Simulation Based Calibration for Rstan/Cmdstanr Models.
#'   Talts, Sean, Michael Betancourt, Daniel Simpson, Aki Vehtari,
#'     and Andrew Gelman. 2020. "Validating Bayesian Inference Algorithms with
#'     Simulation-Based Calibration." https://arxiv.org/abs/1804.06788.
#' @return A data frame with one column per model parameter and a single
#'   row with the rank statistics.
#' @param data List of data for the Stan model.
#' @param draws A data frame with one row per MCMC draw and one column
#'   per model parameter.
#' @examples
#' library(dplyr)
#' library(posterior)
#' library(purrr)
#' library(SBC)
#' library(tibble)
#' data <- simulate_data() # See above for the function definition.
#' draws <- tibble(`beta[1]` = rnorm(100), `beta[2]` = rnorm(100))
#' get_ranks(data = data, draws = draws)
get_ranks <- function(data, draws) {
  draws <- select(
    draws,
    starts_with(names(data$.join_data))
  )
  truth <- map_dbl(
    .x = names(draws),
    .f = ~eval(parse(text = .x), envir = data$.join_data)
  )
  out <- calculate_ranks_draws_matrix(
    variables = truth,
    dm = as_draws_matrix(draws)
  )
  as_tibble(as.list(out))
}

#' @title Plot simulation-based calibration rank statistics.
#' @export
#' @description Plot simulation-based calibration rank statistics.
#' @details For most models, if the model is implemented correctly,
#'   then the SBC rank statistics should be uniformly distributed
#'   for each parameter.
#' @references
#'   Kim, Shinyoung, Hyunji Moon, Martin Modrák, and Teemu Säilynoja. 2022.
#'     SBC: Simulation Based Calibration for Rstan/Cmdstanr Models.
#'   Talts, Sean, Michael Betancourt, Daniel Simpson, Aki Vehtari,
#'     and Andrew Gelman. 2020. "Validating Bayesian Inference Algorithms with
#'     Simulation-Based Calibration." https://arxiv.org/abs/1804.06788.
#' @return A data frame with one column per model parameter and a single
#'   row with the rank statistics.
#' @param ranks Data frame of SBC ranks from get_ranks()
#' @param draws A data frame with one row per MCMC draw and one column
#'   per model parameter.
#' @examples
#' library(dplyr)
#' library(ggplot2)
#' library(tidyr)
#' ranks <- tibble(
#'   `beta[1]` = runif(1e5, min = 0, max = 1e3),
#'   `beta[2]` = runif(1e5, min = 0, max = 1e3)
#' )
#' plot_ranks(ranks = ranks, prefix = "beta")
plot_ranks <- function(ranks, prefix) {
  long <- ranks %>%
    select(starts_with(prefix)) %>%
    pivot_longer(
      cols = everything(),
      names_to = "parameter",
      values_to = "rank"
    )
  ggplot(long) +
    geom_histogram(aes(x = rank), bins = 15) +
    facet_wrap(~parameter) +
    theme_gray(12)
}

# Goal: conduct a simulation study to check that the model
# is implemented correctly in Stan. The pipeline
# simulates from the prior predictive distribution
# and assesses the ability of the model to capture the
# underlying parameter draws in posterior intervals.
library(targets)
library(stantargets)
tar_option_set(
  packages = c(
    "dplyr",
    "ggplot2",
    "posterior",
    "purrr",
    "SBC",
    "tibble",
    "tidyr"
  )
)
tar_source()
# For clustermq configuration with Slurm, see
# https://books.ropensci.org/targets/hpc.html#clustermq-remote-configuration
options(clustermq.scheduler = "multiprocess")
list(
  tar_stan_mcmc_rep_draws(
    name = analysis,
    stan_files = "model.stan",
    data = simulate_data(),
    batches = 100,
    reps = 10,
    chains = 4,
    stdout = nullfile(),
    iter_warmup = 2e3,
    iter_sampling = 2e3,
    transform = get_ranks, # Supply the transform to get SBC ranks.
    format = "feather",
    memory = "transient",
    garbage_collection = TRUE,
    storage = "worker",
    retrieval = "worker"
  ),
  tar_target(beta, plot_ranks(analysis_model, "beta")),
  tar_target(sigma, plot_ranks(analysis_model, "sigma")),
  tar_target(lambda, plot_ranks(analysis_model, "lambda"))
)

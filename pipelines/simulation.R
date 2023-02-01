# Goal: conduct a simulation study to check that the model
# is implemented correctly in Stan. The pipeline
# simulates from the prior predictive distribution
# and assesses the ability of the model to capture the
# underlying parameter draws in posterior intervals.
library(targets)
library(stantargets)
tar_option_set(packages = "dplyr", deployment = "main")
tar_source()
# For clustermq configuration with Slurm, see
# https://books.ropensci.org/targets/hpc.html#clustermq-remote-configuration
options(clustermq.scheduler = "multiprocess")
list(
  tar_stan_mcmc_rep_summary(
    name = analysis,
    stan_files = "model.stan",
    data = simulate_data(),
    batches = 100,
    reps = 10,
    chains = 4,
    stdout = nullfile(),
    iter_warmup = 2e3,
    iter_sampling = 2e3,
    summaries = list(
      ~posterior::quantile2(.x, probs = c(0.025, 0.25, 0.75, 0.975)),
      rhat = posterior::rhat,
      ess_bulk = posterior::ess_bulk,
      ess_tail = posterior::ess_tail
    ),
    memory = "transient",
    garbage_collection = TRUE,
    storage = "worker",
    retrieval = "worker",
    deployment = "worker"
  ),
  tar_target(
    name = convergence,
    command = tibble::tibble(
      max_rhat = max(analysis_model$rhat, na.rm = TRUE),
      min_ess_bulk = min(analysis_model$ess_bulk, na.rm = TRUE),
      min_ess_tail = min(analysis_model$ess_tail, na.rm = TRUE)
    ),
    deployment = "main"
  ),
  tar_target(
    name = coverage,
    command = analysis_model %>%
      filter(grepl("^beta|^sigma|^lambda", variable), !is.na(rhat)) %>%
      group_by(variable) %>%
      summarize(
        coverage_50 = mean(q25 < .join_data & .join_data < q75),
        coverage_95 = mean(q2.5 < .join_data & .join_data < q97.5)
      ),
    deployment = "main"
  )
)

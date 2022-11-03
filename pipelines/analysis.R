library(targets)
library(stantargets)
tar_source()
# For clustermq configuration with Slurm, see
# https://books.ropensci.org/targets/hpc.html#clustermq-remote-configuration
options(clustermq.scheduler = "multisession")
list(
  tar_stan_mcmc(
    name = analysis,
    stan_files = "model.stan",
    data = simulate_data(),
    stdout = nullfile(),
    chains = 4,
    iter_warmup = 2e3,
    iter_sampling = 2e3,
    memory = "transient",
    garbage_collection = TRUE,
    storage = "worker",
    retrieval = "worker"
  )
)

# Goal: analyze a simulated dataset and inspect the results.
rstudioapi::restartSession()
library(targets)
tar_destroy()

# Write the target script.
file.copy("pipelines/analysis.R", "_targets.R", overwrite = TRUE)

# Look at the target script.
tar_edit()

# Inspect the pipeline.
tar_visnetwork()

# Run the pipeline with 2 persistent workers.
tar_make_clustermq(workers = 2)

# Read the model summary.
View(tar_read(analysis_summary_model))

# Read the MCMC draws.
View(tar_read(analysis_draws_model))

# Read the HMC diagnostics
View(tar_read(analysis_diagnostics_model))

# Access the cmdstanr CmdStanMCMC output object.
tar_read(analysis_mcmc_model)

# Access the simulated dataset.
str(tar_read(analysis_data))

# {targets} watches the model file for changes.
tar_read(analysis_file_model)
tar_meta(names = analysis_file_model, fields = data)

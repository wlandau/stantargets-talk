# Goal: conduct a simulation study to check that the model
# is implemented correctly in Stan. The pipeline
# simulates from the prior predictive distribution
# and assesses the ability of the model to capture the
# underlying parameter draws in posterior intervals.
rstudioapi::restartSession()
library(targets)
tar_destroy()

# Write the target script.
file.copy("pipelines/simulation.R", "_targets.R", overwrite = TRUE)

# Look at the target script.
tar_edit()

# Inspect the pipeline.
tar_visnetwork()

# Run the pipeline with 2 persistent workers. First run takes a long time.
# To increase computational power, visit
# https://books.ropensci.org/targets/hpc.html
tar_make_clustermq(workers = 2)

# Read the coverage results. If the model is implemented correctly,
# then cover_95 should be about 0.95 on average, and cover_50
# should be about 0.5 on average.
tar_read(coverage)

# Read the convergence results.
tar_read(convergence)

# The underlying summary output is stored in separate
# files and aggregated when the need arises.
# You can load one or more branches (batches).
tar_read(analysis_model, branches = c(2, 3))

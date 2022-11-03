rstudioapi::restartSession()
library(targets)
tar_destroy()

# Write the target script.
file.copy("pipelines/simulation.R", "_targets.R", overwrite = TRUE)

# Look at the target script.
tar_edit()

# Inspect the pipeline.
tar_visnetwork()

# Run the pipeline with 2 persistent workers.
tar_make_clustermq(workers = 2)

# Read the coverage results. Is the model well calibrated?
tar_read(coverage)

# Read the convergence results.
tar_read(convergence)

# The underlying summary output is stored in separate
# files and aggregated when the need arises.
# You can load one or more branches (batches).
tar_read(analysis_model, branches = c(2, 3))

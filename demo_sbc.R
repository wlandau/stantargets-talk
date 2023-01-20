# Goal: conduct simulation-based calibration (SBC)
# to check that the model is implemented correctly in Stan.
# The pipeline simulates from the prior predictive distribution
# and assesses the alignment of prior parameter draws
# with posterior MCMC draws.
rstudioapi::restartSession()
library(targets)
tar_destroy()

# Write the target script.
file.copy("pipelines/sbc.R", "_targets.R", overwrite = TRUE)

# Look at the target script.
tar_edit()

# Inspect the pipeline.
tar_visnetwork()

# Run the pipeline with 2 persistent workers.
# Takes a long time.
tar_make_clustermq(workers = 2)

# Plot the SBC rank statistics of the fixed effect parameters.
# Should be uniformly distributed if the model is impelemented correctly.
tar_read(beta)

# Same for the residual standard deviation parameters.
tar_read(sigma)

# Same for the correlation matrix parameters.
tar_read(lambda)

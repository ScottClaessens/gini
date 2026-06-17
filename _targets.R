library(stantargets)
library(targets)
library(tarchetypes)

tar_option_set(packages = c("cmdstanr", "tidyverse"))
tar_source()

list(
  
  # ─────────────────────────────────────────
  # Load data
  # ─────────────────────────────────────────
  
  tar_target(data_file, "data/SiteGiniLevel.csv", format = "file"),
  tar_target(data, load_data(data_file)),
  
  # ─────────────────────────────────────────
  # Fit Stan model
  # ─────────────────────────────────────────
  
  tar_stan_mcmc(
    fit,
    stan_files = "stan/model.stan",
    data = get_data_list(data),
    seed = 1
  )
  
)

library(stantargets)
library(targets)
library(tarchetypes)

tar_option_set(packages = c("bayesplot", "cmdstanr", "patchwork", "tidyverse"))
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
    parallel_chains = 4,
    seed = 1
  ),
  
  # ─────────────────────────────────────────
  # Plot model predictions and checks
  # ─────────────────────────────────────────
  
  tar_target(plot_predictions, plot_model_predictions(data, fit_draws_model)),
  tar_target(plot_pp_check, plot_predictive_check(data, fit_draws_model))

)

library(stantargets)
library(targets)
library(tarchetypes)

tar_option_set(
  packages = c("bayesplot", "cmdstanr", "patchwork", "tidybayes", "tidyverse")
)
tar_source()

list(
  
  # ─────────────────────────────────────────
  # Load data
  # ─────────────────────────────────────────
  
  tar_target(file_gini, "data/gini/SiteGiniLevel.csv", format = "file"),
  tar_target(file_hyde, "data/hyde/hyde.csv", format = "file"),
  tar_target(data, load_data(file_gini, file_hyde)),
  
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
  tar_target(plot_pp_check, plot_predictive_check(data, fit_draws_model)),
  
  # ─────────────────────────────────────────
  # Fit to synthetic data and plot
  # ─────────────────────────────────────────
  
  tar_stan_mcmc(
    sim,
    stan_files = "stan/model.stan",
    data = get_data_list(generate_synthetic_data()),
    parallel_chains = 4,
    seed = 1
  ),
  tar_target(plot_sim, plot_simulation(sim_draws_model))
  

)

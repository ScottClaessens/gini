library(stantargets)
library(targets)
library(tarchetypes)

tar_option_set(packages = c("cmdstanr", "loo", "patchwork", "tidyverse"))
tar_source()

list(
  
  # ─────────────────────────────────────────
  # Load data
  # ─────────────────────────────────────────
  
  tar_target(data_file, "data/SiteGiniLevel.csv", format = "file"),
  tar_target(data, read_csv(data_file, show_col_types = FALSE)),
  
  # ─────────────────────────────────────────
  # Plot Gini and population size over time
  # ─────────────────────────────────────────
  
  tar_target(plot_gini_global, plot_gini_time(data)),
  tar_target(plot_gini_by_region, plot_gini_time(data, split_by_region = TRUE)),
  tar_target(plot_population_global, plot_population_time(data)),
  tar_target(
    plot_population_by_region,
    plot_population_time(data, split_by_region = TRUE)
  ),
  
  # ─────────────────────────────────────────
  # Baseline model
  # ─────────────────────────────────────────
  
  # prior predictive check
  tar_stan_mcmc(
    name = prior_baseline,
    stan_files = "stan/model_baseline.stan",
    data = get_stan_data_list(data, prior_only = TRUE),
    parallel_chains = 4,
    seed = 1
  ),
  tar_target(
    prior_check_baseline,
    plot_predictive_check(data, prior_baseline_draws_model_baseline,
                          filename = "plots/prior_check_baseline.pdf")
  ),
  
  # fit model
  tar_stan_mcmc(
    name = fit_baseline,
    stan_files = "stan/model_baseline.stan",
    data = get_stan_data_list(data),
    parallel_chains = 4,
    seed = 1
  ),
  
  # loo-cv
  tar_target(loo_baseline, fit_baseline_mcmc_model_baseline$loo()),
  
  # posterior predictive check
  tar_target(
    pred_check_baseline,
    plot_predictive_check(data, fit_baseline_draws_model_baseline,
                          filename = "plots/pred_check_baseline.pdf")
  ),
  
  # ─────────────────────────────────────────
  # SDE model
  # ─────────────────────────────────────────
  
  # prior predictive check
  tar_stan_mcmc(
    name = prior_sde,
    stan_files = "stan/model_sde.stan",
    data = get_stan_data_list(data, prior_only = TRUE),
    parallel_chains = 4,
    seed = 1
  ),
  tar_target(
    prior_check_sde,
    plot_predictive_check(data, prior_sde_draws_model_sde,
                          filename = "plots/prior_check_sde.pdf")
  ),
  
  # fit model
  tar_stan_mcmc(
    name = fit_sde,
    stan_files = "stan/model_sde.stan",
    data = get_stan_data_list(data),
    parallel_chains = 4,
    seed = 1
  ),
  
  # loo-cv
  tar_target(loo_sde, fit_sde_mcmc_model_sde$loo()),
  
  # posterior predictive check
  tar_target(
    pred_check_sde,
    plot_predictive_check(data, fit_sde_draws_model_sde,
                          filename = "plots/pred_check_sde.pdf")
  )
  
)

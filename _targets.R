library(stantargets)
library(targets)
library(tarchetypes)

tar_option_set(packages = c("cmdstanr", "loo", "tidyverse"))
tar_source()

list(
  
  # load data
  tar_target(data_file, "data/SiteGiniLevel.csv", format = "file"),
  tar_target(data, read_csv(data_file, show_col_types = FALSE)),
  
  # plot gini and population over time
  tar_target(plot_gini_global, plot_gini_time(data)),
  tar_target(plot_gini_by_region, plot_gini_time(data, split_by_region = TRUE)),
  tar_target(plot_population_global, plot_population_time(data)),
  tar_target(
    plot_population_by_region,
    plot_population_time(data, split_by_region = TRUE)
  ),
  
  # get stan data list
  tar_target(stan_data_list, get_stan_data_list(data)),
  
  # fit baseline model
  tar_stan_mcmc(
    name = fit_baseline,
    stan_files = "stan/model_baseline.stan",
    data = stan_data_list,
    parallel_chains = 4,
    seed = 1
  ),
  tar_target(loo_baseline, fit_baseline_mcmc_model_baseline$loo()),
  
  # fit sde model
  tar_stan_mcmc(
    name = fit_sde,
    stan_files = "stan/model_sde.stan",
    data = stan_data_list,
    parallel_chains = 4,
    seed = 1
  ),
  tar_target(loo_sde, fit_sde_mcmc_model_sde$loo())
  
)

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
  # fit stan model
  tar_stan_mcmc(
    name = fit,
    stan_files = "stan/model.stan",
    data = stan_data_list,
    parallel_chains = 4,
    seed = 1
  ),
  # calculate loo-cv
  tar_target(loo, fit_mcmc_model$loo())
  
)

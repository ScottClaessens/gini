library(stantargets)
library(targets)
library(tarchetypes)

tar_option_set(packages = c("cmdstanr", "tidyverse"))
tar_source()

list(
  # load data
  tar_target(data_file, "data/SiteGiniLevel.csv", format = "file"),
  tar_target(data, read_csv(data_file, show_col_types = FALSE)),
  # plot gini and polity population over time
  tar_target(plot, plot_gini_population(data)),
  # fit stan model
  tar_target(stan_data_list, get_stan_data_list(data)),
  tar_stan_mcmc(fit, stan_files = "stan/model.stan", data = stan_data_list),
  # plot posterior predictive check
  tar_target(plot_post, plot_post_pred_check(stan_data_list, fit_draws_model))
)

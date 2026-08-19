options(tidyverse.quiet = TRUE)
library(stantargets)
library(targets)
library(tarchetypes)
library(tidyverse)

tar_option_set(
  packages = c("bayesplot", "cmdstanr", "patchwork", "tidybayes", "tidyverse")
)
tar_source()

list(
  
  # ────────────────────────────────────────────
  # Load data
  # ────────────────────────────────────────────
  
  tar_target(file_gini, "data/gini/SiteGiniLevel.csv", format = "file"),
  tar_target(file_hyde, "data/hyde/hyde.csv", format = "file"),
  tar_target(data, load_data(file_gini, file_hyde)),
  
  # ────────────────────────────────────────────
  # Fit Stan model
  # ────────────────────────────────────────────
  
  tar_stan_mcmc(
    fit,
    stan_files = "stan/model.stan",
    data = get_data_list(data),
    iter_warmup = 500,
    iter_sampling = 500,
    parallel_chains = 4,
    seed = 1234
  ),
  
  # ────────────────────────────────────────────
  # Plot model results, predictions, and checks
  # ────────────────────────────────────────────
  
  tar_target(
    plot_pp_check,
    plot_predictive_check(data, fit_draws_model)
  ),
  tar_target(
    plot_pred_global,
    plot_global_trajectories(fit_draws_model)
  ),
  tar_map(
    values = tibble(
      variable = c("pop_size", "cropland", "urban", "gini")
    ),
    tar_target(
      plot_pred_regional,
      plot_regional_predictions(data, fit_draws_model, variable)
    )
  ),
  tar_target(plot_delta_gini, plot_gini_change(fit_draws_model)),
  tar_target(plot_rates, plot_model_rates(fit_draws_model)),
  tar_map(
    tibble(effect = paste0("beta", c("P", "C", "U"))),
    tar_target(
      plot_effects,
      plot_varying_effects(data, fit_draws_model, effect)
    )
  ),
  
  # ────────────────────────────────────────────
  # Produce report
  # ────────────────────────────────────────────
  
  tar_quarto(report, "quarto/report.qmd")

)

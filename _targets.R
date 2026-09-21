options(tidyverse.quiet = TRUE)
library(crew)
library(crew.cluster)
library(stantargets)
library(targets)
library(tarchetypes)
library(tidyverse)
tar_option_set(
  packages = c("bayesplot", "cmdstanr", "patchwork", "tidybayes", "tidyverse"),
  controller = crew_controller_slurm(
    workers = 1,
    options_metrics = crew_options_metrics(
      path = "/dev/stdout",
      seconds_interval = 60
    ),
    options_cluster = crew_options_slurm(
      script_lines = c(
        "#SBATCH --account=arch039044",
        "module load languages/R/4.5.1"
      ),
      memory_gigabytes_required = 200,
      cpus_per_task = 8,
      time_minutes = 10 * 24 * 60,
      log_output = "crew_log_%A.out",
      log_error = "crew_log_%A.err"
    )
  )
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
    iter_warmup = 1000,
    iter_sampling = 500,
    chains = 8,
    parallel_chains = 8,
    seed = 1
  ),
  
  # ────────────────────────────────────────────
  # Plot model results, predictions, and checks
  # ────────────────────────────────────────────
  
  tar_target(plot_pp_check, plot_predictive_check(data, fit_draws_model)),
  tar_target(plot_pred_global, plot_global_trajectories(fit_draws_model)),
  tar_map(
    values = tibble(
      variable = c("pop_size", "cropland", "irrigated", "urban", "gini")
    ),
    tar_target(
      plot_pred_regional,
      plot_regional_predictions(data, fit_draws_model, variable)
    )
  ),
  tar_target(plot_delta_gini, plot_gini_change(fit_draws_model)),
  tar_target(plot_rates, plot_model_rates(fit_draws_model)),
  tar_map(
    tibble(effect = paste0("beta", c("P", "C", "I", "U"))),
    tar_target(
      plot_effects,
      plot_varying_effects(data, fit_draws_model, effect)
    )
  )#,
  
  # ────────────────────────────────────────────
  # Produce report
  # ────────────────────────────────────────────
  
  #tar_quarto(report, "quarto/report.qmd")

)

#' Plot posterior predictive check for Stan model
#'
#' @param data Tibble of raw data
#' @param fit_draws_model Tibble of posterior draws from the model
#'
#' @returns A ggplot object
#'
plot_predictive_check <- function(data, fit_draws_model) {
  
  # y and yrep for pop size
  y <- data$pop_size[!is.na(data$pop_size)] + 0.001
  yrep <-
    fit_draws_model |>
    dplyr::select(starts_with("pop_size_rep")) |>
    as.matrix()
  
  # plot for pop size
  out <- 
    bayesplot::ppc_dens_overlay(y, yrep[1:50, ]) +
    scale_x_continuous(
      name = "Population size (log scale)",
      transform = "log",
      breaks = c(1e-02, 1e+01, 1e+04)
    )
  
  # save
  ggsave(
    plot = out,
    filename = "plots/pp_check.pdf",
    height = 4,
    width = 4
  )
  
  # return
  out
  
}
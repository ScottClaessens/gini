#' Plot posterior predictive check for Stan model
#'
#' @param data Tibble of raw data
#' @param fit_draws_model Tibble of posterior draws from the model
#'
#' @returns A ggplot object
#'
plot_predictive_check <- function(data, fit_draws_model) {
  
  # y and yrep for pop size
  y <- data$pop_size[!is.na(data$pop_size)]
  yrep <-
    fit_draws_model |>
    dplyr::select(starts_with("pop_size_rep")) |>
    as.matrix()
  
  # plot for pop size
  pA <- 
    bayesplot::ppc_dens_overlay(y, yrep[1:50, ]) +
    scale_x_continuous(
      name = "Population size (log + 1)",
      transform = "log1p",
      breaks = c(0, 100, 10000)
    )
  
  # y and yrep for cropland
  y <- data$cropland[!is.na(data$cropland)]
  yrep <-
    fit_draws_model |>
    dplyr::select(starts_with("cropland_rep")) |>
    as.matrix()
  
  # plot for cropland
  pB <- 
    bayesplot::ppc_dens_overlay(y, yrep[1:50, ]) +
    scale_x_continuous(
      name = "Cropland (log + 1)",
      transform = "log1p",
      limits = c(0, 150),
      breaks = c(0, 10, 100)
    )
  
  # y and yrep for urban
  y <- data$urban[!is.na(data$urban)]
  yrep <-
    fit_draws_model |>
    dplyr::select(starts_with("urban_rep")) |>
    as.matrix()
  
  # plot for urban
  pC <- 
    bayesplot::ppc_dens_overlay(y, yrep[1:50, ]) +
    scale_x_continuous(
      name = "Urban area (log + 1)",
      transform = "log1p",
      limits = c(0, 5),
      breaks = c(0, 1, 2, 4)
    )
  
  # y and yrep for gini
  y <- data$gini[!is.na(data$gini)]
  yrep <-
    fit_draws_model |>
    dplyr::select(starts_with("gini_rep")) |>
    as.matrix()
  
  # plot for gini
  pD <- 
    bayesplot::ppc_dens_overlay(y, yrep[1:50, ]) +
    scale_x_continuous(
      name = "Gini",
      breaks = c(0.25, 0.5, 0.75)
    )
  
  # put together
  out <- 
    ((pA + pB) / (pC + pD)) +
    plot_layout(guides = "collect")
  
  # cleanup
  rm(data, fit_draws_model)
  
  # save
  ggsave(
    plot = out,
    filename = "plots/pp_check.pdf",
    height = 5,
    width = 5
  )
  
  # return
  out
  
}
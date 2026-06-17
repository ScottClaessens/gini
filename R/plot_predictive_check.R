#' Plot posterior predictive check for Stan model
#'
#' @param data Tibble of raw data
#' @param fit_draws_model Tibble of posterior draws from the model
#'
#' @returns A ggplot object
#'
plot_predictive_check <- function(data, fit_draws_model) {
  
  # y and yrep for gini
  y <- data$gini
  yrep <- 
    fit_draws_model |>
    dplyr::select(starts_with("gini_rep")) |>
    as.matrix()
  
  # plot for gini
  pA <- 
    bayesplot::ppc_dens_overlay(y, yrep[1:50, ]) +
    scale_x_continuous(
      name = "Gini",
      limits = c(0, 1)
    )
  
  # y and yrep for pop size
  y <- data$polity_pop
  yrep <-
    fit_draws_model |>
    dplyr::select(starts_with("pop_size_rep")) |>
    as.matrix()
  
  # plot for pop size
  obs <- !is.na(y)
  pB <- 
    bayesplot::ppc_dens_overlay(y[obs], yrep[1:50, obs]) +
    scale_x_continuous(
      name = "Population size (log scale)",
      transform = "log",
      breaks = c(1e-02, 1e+01, 1e+04, 1e+07, 1e+10)
    )
  
  # combine
  out <- (pA / pB) + plot_layout(guides = "collect")
  
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
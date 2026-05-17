#' Plot predictive check from a fitted model
#'
#' @param data Tibble of data
#' @param draws Draws from a fitted model
#'
#' @returns A ggplot object
#'
plot_predictive_check <- function(data, draws, filename) {
  
  # get data in stan list form
  data <- get_stan_data_list(data)
  
  # gini plot
  y <- data$gini
  
  yrep <- 
    draws |>
    dplyr::select(starts_with("gini_rep")) |>
    slice(1:25) |>
    as.matrix()
  
  p1 <- 
    bayesplot::ppc_dens_overlay(y, yrep) +
    scale_x_continuous(
      name = "Gini",
      limits = c(0, 1)
    )
  
  # population size plot
  y <- ifelse(data$pop_scaled == -9999, NA, data$pop_scaled)
  
  yrep <-
    draws |>
    dplyr::select(starts_with("pop_scaled_rep")) |>
    slice(1:25) |>
    as.matrix()
  
  p2 <- 
    bayesplot::ppc_dens_overlay(
      y = log(y[!is.na(y)]),
      yrep = log(yrep[, !is.na(y)])
    ) +
    scale_x_continuous(
      name = "Population size\n(logged and scaled)",
      breaks = c(-20, -10, 0),
      limits = c(-25, 5)
    )
  
  # put together
  out <- 
    p1 + p2 +
    patchwork::plot_layout(guides = "collect")
  
  # save
  ggsave(
    filename = filename,
    plot = out,
    height = 3,
    width = 5
  )
  
  # return
  out
  
}

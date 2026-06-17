#' Plot time-series predictions from Stan model
#'
#' @param data Tibble of raw data
#' @param fit_draws_model Tibble of posterior draws from the model
#'
#' @returns A ggplot object
#'
plot_model_predictions <- function(data, fit_draws_model) {
  
  # plot predictions for gini
  pA <-
    fit_draws_model |>
    dplyr::select(starts_with("latent_rep") & ends_with("1]")) |>
    pivot_longer(
      cols = everything(),
      names_to = "id",
      names_pattern = "latent_rep\\[(.*),1]"
    ) |>
    mutate(id = as.numeric(id)) |>
    group_by(id) |>
    summarise(
      median = median(plogis(value)),
      lower = quantile(plogis(value), 0.025),
      upper = quantile(plogis(value), 0.975)
    ) |>
    mutate(
      time = seq(0, 1, length.out = 100),
      time_before_present = time * (max(data$date) - min(data$date)) + 
        min(data$date) - 2026
    ) |>
    
    ggplot() +
    geom_point(
      data = data,
      aes(
        x = date - 2026,
        y = gini
      ),
      size = 0.2
    ) +
    geom_ribbon(
      aes(
        x = time_before_present,
        ymin = lower,
        ymax = upper
      ),
      fill = "lightblue",
      alpha = 0.8
    ) +
    geom_line(
      aes(
        x = time_before_present,
        y = median
      )
    ) +
    scale_y_continuous(
      name = "Wealth inequality (Gini)",
      limits = c(0, 1)
    ) +
    scale_x_continuous(
      name = "Time before present"
    ) +
    theme_classic()
  
  # plot predictions for population size
  pB <-
    fit_draws_model |>
    dplyr::select(starts_with("latent_rep") & ends_with("2]")) |>
    pivot_longer(
      cols = everything(),
      names_to = "id",
      names_pattern = "latent_rep\\[(.*),2]"
    ) |>
    mutate(id = as.numeric(id)) |>
    group_by(id) |>
    summarise(
      median = median(exp(value)),
      lower = quantile(exp(value), 0.025),
      upper = quantile(exp(value), 0.975)
    ) |>
    mutate(
      time = seq(0, 1, length.out = 100),
      time_before_present = time * (max(data$date) - min(data$date)) + 
        min(data$date) - 2026
    ) |>
    
    ggplot() +
    geom_point(
      data = data,
      aes(
        x = date - 2026,
        y = polity_pop
      ),
      size = 0.2
    ) +
    geom_ribbon(
      aes(
        x = time_before_present,
        ymin = lower,
        ymax = upper
      ),
      fill = "lightblue",
      alpha = 0.8
    ) +
    geom_line(
      aes(
        x = time_before_present,
        y = median
      )
    ) +
    scale_y_continuous(
      name = "Population size",
      transform = "log",
      breaks = c(1e+01, 1e+03, 1e+05, 1e+07),
      labels = scales::comma
    ) +
    scale_x_continuous(
      name = "Time before present"
    ) +
    theme_classic()
  
  # combine
  out <- (pA / pB) + plot_layout(axis_titles = "collect_x")
  
  # save
  ggsave(
    plot = out,
    filename = "plots/predictions.pdf",
    height = 5,
    width = 5
  )
  
  # return
  out

}

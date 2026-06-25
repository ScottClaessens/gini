#' Plot regional time-series predictions from Stan model
#'
#' @param data Tibble of raw data
#' @param fit_draws_model Tibble of posterior draws from the model
#'
#' @returns A ggplot object
#'
plot_regional_predictions <- function(data, fit_draws_model) {
  
  out <-
    fit_draws_model |>
    dplyr::select(starts_with("regional_latent_rep")) |>
    pivot_longer(
      cols = everything(),
      names_to = c("region", "t", "variable"),
      names_pattern = paste0("regional_latent_rep\\[(.*),(.*),(.*)]")
    ) |>
    mutate(
      region = levels(factor(data$bigregion))[as.numeric(region)],
      t = seq(0, 1, length.out = 100)[as.numeric(t)],
      variable = c("gini", "pop_size", "cropland")[as.numeric(variable)],
      value = ifelse(variable == "gini", plogis(value), log(exp(value) + 1))
    ) |>
    group_by(region, t, variable) |>
    summarise(
      median = median(value),
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975),
      .groups = "drop"
    ) |>
    mutate(
      time_before_present = t * (max(data$date) - min(data$date)) + 
        min(data$date) - 2026
    ) |>
    
    ggplot() +
    geom_point(
      data = data |>
        pivot_longer(
          cols = c("gini", "pop_size", "cropland"),
          names_to = "variable"
        ) |>
        mutate(
          region = bigregion,
          value = case_when(
            variable == "gini" ~ value,
            variable == "pop_size" ~ log(value + 1),
            variable == "cropland" ~ log(value + 1)
          )
        ),
      aes(
        x = (date - 2026) / 1000,
        y = value
      ),
      size = 0.2
    ) +
    geom_ribbon(
      aes(
        x = time_before_present / 1000,
        ymin = lower,
        ymax = upper
      ),
      fill = "lightgrey",
      alpha = 0.8
    ) +
    geom_line(
      aes(
        x = time_before_present / 1000,
        y = median
      )
    ) +
    facet_grid(
      factor(variable, levels = c("gini", "pop_size", "cropland")) ~ region,
      scales = "free_y"
    ) +
    scale_x_continuous(
      name = "Time before present (ky)",
      breaks = c(-10, -5, 0)
    ) +
    theme_classic()
  
  # save
  ggsave(
    plot = out,
    filename = "plots/regional_predictions.pdf",
    height = 5,
    width = 8
  )
  
  # return
  out

}

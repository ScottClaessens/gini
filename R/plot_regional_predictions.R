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
      names_to = c("region", "date", "variable"),
      names_pattern = paste0("regional_latent_rep\\[(.*),(.*),(.*)]")
    ) |>
    mutate(
      region = levels(factor(data$subregion))[as.numeric(region)],
      date = seq(
        min(data$date), max(data$date), length.out = 100
      )[as.numeric(date)],
      variable = c("pop_size", "cropland")[as.numeric(variable)]
    ) |>
    group_by(region, date, variable) |>
    summarise(
      median = median(value),
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975),
      .groups = "drop"
    ) |>
    mutate(time_before_present = date - 2026) |>
    filter(variable == "cropland") |>
    
    ggplot() +
    geom_point(
      data = mutate(data, region = subregion, value = log(cropland + 0.001)),
      aes(
        x = (date - 2026) / 1000,
        y = value
      ),
      size = 0.05,
      colour = "black"
    ) +
    geom_ribbon(
      aes(
        x = time_before_present / 1000,
        ymin = lower,
        ymax = upper
      ),
      fill = "lightblue",
      alpha = 0.8
    ) +
    geom_line(
      aes(
        x = time_before_present / 1000,
        y = median
      )
    ) +
    facet_wrap(~ region) +
    scale_x_continuous(
      name = "Time before present (ky)",
      breaks = c(-10, -5, 0)
    ) +
    scale_y_continuous(
      name = "Cropland (log)"
    ) +
    theme_classic() +
    theme(strip.text = element_text(size = 6))
  
  # cleanup
  rm(data, fit_draws_model)
  
  # save
  ggsave(
    plot = out,
    filename = "plots/regional_predictions.pdf",
    height = 8,
    width = 8
  )
  
  # return
  out
  
}
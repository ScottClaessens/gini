#' Plot regional time-series predictions from Stan model
#'
#' @param data Tibble of raw data
#' @param fit_draws_model Tibble of posterior draws from the model
#' @param variable Variable to predict
#'
#' @returns A ggplot object
#'
plot_regional_predictions <- function(data, fit_draws_model, 
                                      variable = "pop_size") {
  
  out <-
    fit_draws_model |>
    dplyr::select(starts_with("regional_latent_rep")) |>
    mutate(draw = 1:nrow(fit_draws_model)) |>
    pivot_longer(
      cols = !draw,
      names_to = c("region", "date", "var"),
      names_pattern = paste0("regional_latent_rep\\[(.*),(.*),(.*)]")
    ) |>
    mutate(
      region = levels(factor(data$subregion))[as.numeric(region)],
      date   = seq(-10000, 2000, length.out = 121)[as.numeric(date)],
      var    = c("logit_pop", "raw_pop", 
                 "logit_crop", "raw_crop")[as.numeric(var)]
    ) |>
    pivot_wider(
      names_from = "var",
      values_from = "value"
    ) |>
    mutate(
      pop_size = ifelse(
        rbinom(n(), 1, plogis(logit_pop)) == 0, 0, exp(raw_pop)
      ),
      cropland = ifelse(
        rbinom(n(), 1, plogis(logit_crop)) == 0, 0, exp(raw_crop)
      )
    ) |>
    dplyr::select(c(region, date, pop_size, cropland)) |>
    pivot_longer(
      cols = c(pop_size, cropland),
      names_to = "var"
    ) |>
    group_by(region, date, var) |>
    summarise(
      median = median(log(value + 1)),
      lower = quantile(log(value + 1), 0.025),
      upper = quantile(log(value + 1), 0.975),
      .groups = "drop"
    ) |>
    mutate(time_before_present = date - 2026) |>
    filter(var == variable) |>
    
    ggplot() +
    geom_point(
      data = mutate(
        data,
        region = subregion,
        value = log(!!sym(variable) + 1)
      ),
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
      name = ifelse(
        variable == "pop_size",
        "Population size (log + 1)",
        "Cropland (log + 1)"
      )
    ) +
    theme_classic() +
    theme(strip.text = element_text(size = 4))
  
  # cleanup
  rm(data, fit_draws_model)
  
  # save
  ggsave(
    plot = out,
    filename = paste0("plots/regional_predictions_", variable, ".pdf"),
    height = 8,
    width = 8
  )
  
  # return
  out
  
}
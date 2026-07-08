#' Plot regional time-series predictions from Stan model
#'
#' @param data Tibble of raw data
#' @param fit_draws_model Tibble of posterior draws from the model
#' @param variable Variable to predict
#' @param pred_transform Function to transform model predictions
#' @param data_transform Function to transform raw data
#'
#' @returns A ggplot object
#'
plot_regional_predictions <- function(data, fit_draws_model, 
                                      variable = "pop_size",
                                      pred_transform = function(x) log(x + 1),
                                      data_transform = function(x) log(x + 1)) {
  
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
                 "logit_crop", "raw_crop", "gini")[as.numeric(var)]
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
      ),
      gini = plogis(gini)
    ) |>
    dplyr::select(c(region, date, pop_size, cropland, gini)) |>
    pivot_longer(
      cols = c(pop_size, cropland, gini),
      names_to = "var"
    ) |>
    group_by(region, date, var) |>
    summarise(
      median = median(pred_transform(value)),
      lower = quantile(pred_transform(value), 0.025),
      upper = quantile(pred_transform(value), 0.975),
      .groups = "drop"
    ) |>
    mutate(time_before_present = date - 2026) |>
    filter(var == variable) |>
    
    ggplot() +
    geom_point(
      data = mutate(
        data,
        region = subregion,
        value = data_transform(!!sym(variable))
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
      name = case_when(
        variable == "pop_size" ~ "Population size (log + 1)",
        variable == "cropland" ~ "Cropland (log + 1)",
        variable == "gini" ~ "Gini"
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
#' Plot time-series predictions from Stan model
#'
#' @param data Tibble of raw data
#' @param fit_draws_model Tibble of posterior draws from the model
#'
#' @returns A ggplot object
#'
plot_model_predictions <- function(data, fit_draws_model) {
  
  # plotting function
  plot_fun <- function(var_id, var_name, pred_transform,
                       data_transform, ylab, ylim) {
    fit_draws_model |>
      dplyr::select(
        starts_with("latent_rep") & ends_with(paste0(var_id, "]"))
      ) |>
      pivot_longer(
        cols = everything(),
        names_to = "id",
        names_pattern = paste0("latent_rep\\[(.*),", var_id, "]")
      ) |>
      mutate(id = as.numeric(id)) |>
      group_by(id) |>
      summarise(
        median = median(pred_transform(value)),
        lower = quantile(pred_transform(value), 0.025),
        upper = quantile(pred_transform(value), 0.975)
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
          y = data_transform(!!sym(var_name))
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
        name = ylab,
        limits = ylim
      ) +
      scale_x_continuous(
        name = "Time before present"
      ) +
      theme_classic()
  }
  
  # plot predictions
  pA <- 
    plot_fun(
      var_id = 1,
      var_name = "gini",
      pred_transform = plogis,
      data_transform = function(x) x,
      ylab = "Wealth inequality (Gini)",
      ylim = c(0, 1)
    )
  pB <- 
    plot_fun(
      var_id = 2,
      var_name = "pop_size",
      pred_transform = log,
      data_transform = log,
      ylab = "Population size (log)",
      ylim = NULL
    )
  pC <- 
    plot_fun(
      var_id = 3,
      var_name = "cropland",
      pred_transform = log,
      data_transform = function(x) log(x + 0.001),
      ylab = "Cropland (log)",
      ylim = NULL
    )
  
  # combine
  out <- (pA / pB / pC) + plot_layout(axis_titles = "collect_x")
  
  # save
  ggsave(
    plot = out,
    filename = "plots/predictions.pdf",
    height = 7,
    width = 5
  )
  
  # return
  out

}

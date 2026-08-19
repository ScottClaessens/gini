#' Plot rates from Stan model
#' 
#' Plot rates of population growth and production of cropland, and urban areas 
#' from Stan model
#'
#' @param fit_draws_model Tibble of posterior draws from the model
#'
#' @returns A ggplot object
#'
plot_model_rates <- function(fit_draws_model) {
  
  # vector of variable names
  var_names <- c(
    "pop"  = "Average rate of population\ngrowth per capita",
    "crop" = "Average rate of cropland\nproduction per capita",
    "urb"  = "Average rate of urban area\nproduction per capita"
  )
  
  # vector of period names
  period_names <- c(
    "1" = "10,000 BCE - 0CE",
    "2" = "0 CE - 1600 CE",
    "3" = "1600 CE - Present"
  )
  
  # plot
  p <-
    tibble(
      pop_1  = exp(fit_draws_model$`theta[1,2]`),
      pop_2  = exp(fit_draws_model$`theta[2,2]`),
      pop_3  = exp(fit_draws_model$`theta[3,2]`),
      crop_1 = exp(fit_draws_model$`theta[1,4]`),
      crop_2 = exp(fit_draws_model$`theta[2,4]`),
      crop_3 = exp(fit_draws_model$`theta[3,4]`),
      urb_1  = exp(fit_draws_model$`theta[1,6]`),
      urb_2  = exp(fit_draws_model$`theta[2,6]`),
      urb_3  = exp(fit_draws_model$`theta[3,6]`),
    ) |>
    pivot_longer(
      cols = everything(),
      names_to = c("variable", "period"),
      names_sep = "_"
    ) |>
    mutate(
      variable = factor(var_names[variable], levels = var_names),
      period = factor(period_names[period], levels = period_names)
    ) |>
    
    ggplot(
      aes(
        x = value,
        y = fct_rev(period)
      )
    ) +
    ggdist::stat_pointinterval() +
    facet_wrap(
      . ~ variable,
      ncol = 2,
      nrow = 2
    ) +
    labs(
      x = NULL,
      y = "Period",
      caption = "Rates measured per century"
    ) +
    scale_x_continuous(
      breaks = c(0.001, 0.01, 0.1, 1),
      limits = c(NA, 1),
      transform = "log",
      labels = scales::label_number(drop0trailing = TRUE)
    ) +
    theme_classic()
  
  # save
  ggsave(
    plot = p,
    file = "plots/rates.pdf",
    height = 5,
    width = 5
  )
  
  # return
  p
  
}

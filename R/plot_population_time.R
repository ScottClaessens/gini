#' Plot log population size over time
#'
#' @param data Tibble of data
#' @param split_by_region Logical. Should the plot split by region?
#'
#' @returns A ggplot object
#'
plot_population_time <- function(data, split_by_region = FALSE) {
  
  p <-
    data |>
    # remove some sites for plot
    filter(PolityPop > 0 & BeginDate > -8000) |>
    # plot
    ggplot(
      mapping = aes(
        x = BeginDate - 2026,
        xend = EndDate - 2026,
        y = log(PolityPop),
        colour = Gini
      )
    ) +
    geom_segment(
      linewidth = 1.8,
      lineend = "round",
      position = position_jitter(
        height = 0.4,
        seed = 1
      )
    ) +
    scale_colour_viridis_c(
      limits = c(0, 1),
      option = "A",
      direction = -1
    ) +
    labs(
      x = "Years before present",
      y = "Log polity population"
    ) +
    theme_classic()
  
  # split by region?
  if (split_by_region) {
    p <-
      p +
      facet_wrap(. ~ Bigregion)
  }
  
  # save
  ggsave(
    filename = paste0(
      "plots/population_",
      ifelse(split_by_region, "by_region", "global"),
      ".pdf"
    ),
    plot = p,
    height = ifelse(split_by_region, 6, 4),
    width = ifelse(split_by_region, 8, 6)
  )
  
  # return
  p
  
}

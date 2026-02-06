# function to plot gini and population over time
plot_gini_population <- function(data) {
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
  # save
  ggsave(
    filename = "plots/gini_population.pdf",
    plot = p,
    height = 4,
    width = 6
  )
  # return
  p
}

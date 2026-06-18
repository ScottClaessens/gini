plot_simulation <- function(sim_draws_model) {
  
  # wrangle posterior draws
  out <-
    sim_draws_model |>
    dplyr::select(starts_with("init"), starts_with("theta"), phi, sigma) |>
    rename(
      k     = `theta[1]`,
      alpha = `theta[2]`,
      beta  = `theta[3]`,
      r     = `theta[4]`
    ) |>
    pivot_longer(
      cols = everything(),
      names_to = "Parameter",
      values_to = "Value"
    ) |>
    
    # plot
    ggplot(
      aes(
        x = Value,
        y = Parameter
      )
    ) +
    stat_pointinterval() +
    geom_point(
      data = tibble(
        Parameter = c("init_gini", "init_pop_size", "k",
                      "alpha", "beta", "r", "phi", "sigma"),
        Value = c(-1, 5, -1, 0, 1, 2, 8, 4)
      ),
      colour = "red",
      shape = 17
    ) +
    theme_classic()
  
  # save
  ggsave(
    filename = "plots/simulation.pdf",
    plot = out,
    height = 2.5,
    width = 4.5
  )
  
  # return
  out
  
}

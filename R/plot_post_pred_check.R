# function to plot posterior predictive check
plot_post_pred_check <- function(stan_data_list, fit_draws_model) {
  # get sequence of years
  year_seq <- seq(
    from = min(stan_data_list$start_year),
    to = max(stan_data_list$end_year),
    length.out = 100
  )
  # get predictions
  post <- fit_draws_model
  preds <- lapply(year_seq, function(x) post$alpha + post$beta * (x / 1000))
  # raw data
  d <- tibble(
    pop_size   = stan_data_list$pop_size,
    start_year = stan_data_list$start_year,
    end_year   = stan_data_list$end_year
  )
  # plot
  p <-
    tibble(
      year = year_seq,
      pred = preds
    ) |>
    rowwise() |>
    mutate(
      estimate = median(pred),
      lower = quantile(pred, 0.025),
      upper = quantile(pred, 0.975)
    ) |>
    # plot
    ggplot() +
    geom_segment(
      data = d,
      mapping = aes(
        x = start_year,
        xend = end_year,
        y = log(pop_size)
      ),
      colour = "grey"
    ) +
    geom_ribbon(
      mapping = aes(
        x = year,
        ymin = lower,
        ymax = upper
      ),
      fill = "grey",
      alpha = 0.4
    ) +
    geom_line(
      mapping = aes(
        x = year,
        y = estimate
      )
    ) +
    scale_x_continuous(
      name = "Years before present",
      labels = function(x) x - 2026,
      breaks = seq(-8000, 0, by = 2000) + 2026
    ) +
    ylab("Log polity population") +
    theme_classic()
  # save
  ggsave(
    filename = "plots/posterior_predictive_check.pdf",
    plot = p,
    height = 4,
    width = 6
  )
  p
}

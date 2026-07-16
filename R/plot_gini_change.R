plot_gini_change <- function(fit_draws_model) {
  
  # get parameters
  alpha <- fit_draws_model$`theta[1,5]`
  betaP <- fit_draws_model$`theta[1,6]`
  betaC <- fit_draws_model$`theta[1,7]`
  
  # plot marginal effect of population size on delta gini
  pA <-
    tibble(logP_plus1 = seq(0, 15, by = 0.1)) |>
    rowwise() |>
    mutate(
      post = list(plogis(alpha + (betaP * logP_plus1)) - 0.5),
      median = median(post),
      lower50 = quantile(post, 0.25),
      upper50 = quantile(post, 0.75),
      lower95 = quantile(post, 0.025),
      upper95 = quantile(post, 0.975)
    ) |>
    
    ggplot(
      aes(
        x = logP_plus1,
        y = median
      )
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    geom_ribbon(
      aes(
        ymin = lower95,
        ymax = upper95
      ),
      fill = "lightblue",
      alpha = 0.5
    ) +
    geom_ribbon(
      aes(
        ymin = lower50,
        ymax = upper50
      ),
      fill = "lightblue",
      alpha = 0.8
    ) +
    geom_line() +
    scale_x_continuous(
      name = "Population size",
      breaks = log(c(0, 1e+02, 1e+04, 1e+06) + 1),
      labels = function(x) scales::comma(exp(x) - 1)
    ) +
    scale_y_continuous(
      name = "Change in Gini",
      limits = c(-0.3, 0.3)
    ) +
    theme_classic()
    
  # plot marginal effect of cropland on delta gini
  pB <-
    tibble(logC_plus1 = seq(0, 5, by = 0.1)) |>
    rowwise() |>
    mutate(
      post = list(plogis(alpha + (betaC * logC_plus1)) - 0.5),
      median = median(post),
      lower50 = quantile(post, 0.25),
      upper50 = quantile(post, 0.75),
      lower95 = quantile(post, 0.025),
      upper95 = quantile(post, 0.975)
    ) |>
    
    ggplot(
      aes(
        x = logC_plus1,
        y = median
      )
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    geom_ribbon(
      aes(
        ymin = lower95,
        ymax = upper95
      ),
      fill = "lightgreen",
      alpha = 0.5
    ) +
    geom_ribbon(
      aes(
        ymin = lower50,
        ymax = upper50
      ),
      fill = "lightgreen",
      alpha = 0.8
    ) +
    geom_line() +
    scale_x_continuous(
      name = "Cropland",
      breaks = log(c(0, 1e+01, 1e+02) + 1),
      labels = function(x) scales::comma(exp(x) - 1)
    ) +
    scale_y_continuous(
      name = "Change in Gini",
      limits = c(-0.3, 0.3)
    ) +
    theme_classic()
  
  # put together
  out <- pA + pB + plot_layout(axis_titles = "collect_y")
  
  # save
  ggsave(
    plot = out,
    file = "plots/gini_change.pdf",
    height = 4,
    width = 6
  )
  
  # return
  out
  
}

#' Plot marginal effects of variables on change in Gini values
#'
#' @param fit_draws_model Tibble of posterior draws from the model.
#'
#' @returns A patchwork of ggplot objects
#'
plot_gini_change <- function(fit_draws_model) {
  
  # get parameters
  alpha <- fit_draws_model$`theta[1,7]`
  betaP <- fit_draws_model$`theta[1,8]`
  betaC <- fit_draws_model$`theta[1,9]`
  betaU <- fit_draws_model$`theta[1,10]`
  
  # function to plot marginal effect of predictor on delta gini
  plot_marginal <- function(beta, xlim, colour, xlab) {
    
    tibble(log_predictor_plus1 = seq(xlim[1], xlim[2], by = 0.1)) |>
      rowwise() |>
      mutate(
        post = list(plogis(alpha + (beta * log_predictor_plus1)) - 0.5),
        median = median(post),
        lower50 = quantile(post, 0.25),
        upper50 = quantile(post, 0.75),
        lower95 = quantile(post, 0.025),
        upper95 = quantile(post, 0.975)
      ) |>
      ggplot(
        aes(
          x = log_predictor_plus1,
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
        fill = colour,
        alpha = 0.5
      ) +
      geom_ribbon(
        aes(
          ymin = lower50,
          ymax = upper50
        ),
        fill = colour,
        alpha = 0.8
      ) +
      geom_line() +
      scale_x_continuous(
        name = xlab,
        breaks = log(c(0, 1e+02, 1e+04, 1e+06) + 1),
        labels = function(x) scales::comma(exp(x) - 1)
      ) +
      scale_y_continuous(
        name = "Change in Gini",
        limits = c(-0.6, 0.6)
      ) +
      theme_classic()
  }
  
  # plot margin effects
  pA <- plot_marginal(betaP, xlim = c(0, 15), colour = "lightblue", 
                      xlab = "Population size")
  pB <- plot_marginal(betaC, xlim = c(0, 5), colour = "lightgreen",
                      xlab = "Cropland")
  pC <- plot_marginal(betaU, xlim = c(0, 5), colour = "#CC93CC",
                      xlab = "Urban area")
  
  # put together
  out <- 
    pA + pB + pC +
    plot_layout(
      nrow = 1,
      ncol = 3,
      axis_titles = "collect_y"
    )
  
  # save
  ggsave(
    plot = out,
    file = "plots/gini_change.pdf",
    height = 3,
    width = 6
  )
  
  # return
  out
  
}

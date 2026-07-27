#' Plot varying effects from the model
#'
#' @param data Tibble of raw data.
#' @param fit_draws_model Tibble of posterior draws from the model.
#' @param effect String. Effect to plot. Either 'betaP' or 'betaC'.
#' @param prob Numeric. Width of inner credible intervals in forest plot.
#' @param prob_outer Numeric. Width of outer credible intervals in forest plot.
#'
#' @returns A ggplot object
#'
plot_varying_effects <- function(data, fit_draws_model, effect, prob = 0.50,
                                 prob_outer = 0.95) {
  
  # set up data frame with varying effects
  d <-
    tibble(
      subregion = unique(data$subregion),
      median = NA,
      lower_in = NA,
      upper_in = NA,
      lower_out = NA,
      upper_out = NA
    )
  
  # get index for effect
  j <- case_when(
    effect == "betaP" ~ 6,
    effect == "betaC" ~ 7,
    TRUE ~ NA
  )
  
  # extract varying effects
  for (r in 1:86) {
    post <- fit_draws_model[[paste0("theta_r[1,", r, ",", j, "]")]]
    d$median[r] <- median(post)
    d$lower_in[r] <- quantile(post, (1 - prob) / 2)
    d$upper_in[r] <- quantile(post, (1 + prob) / 2)
    d$lower_out[r] <- quantile(post, (1 - prob_outer) / 2)
    d$upper_out[r] <- quantile(post, (1 + prob_outer) / 2)
  }
  
  # plot
  out <-
    ggplot(
      data = d,
      mapping = aes(
        x = median,
        y = fct_reorder(subregion, median)
      )
    ) +
    geom_vline(xintercept = 0) +
    geom_linerange(
      aes(
        xmin = lower_out,
        xmax = upper_out
      )
    ) +
    geom_linerange(
      aes(
        xmin = lower_in,
        xmax = upper_in
      ),
      linewidth = 1.05
    ) +
    geom_point() +
    labs(
      x = case_when(
        effect == "betaP" ~ "Direct effect of population size on inequality",
        effect == "betaC" ~ "Direct effect of cropland on inequality",
        TRUE ~ NA
      ),
      y = NULL
    ) +
    theme_classic() +
    theme(axis.text.y = element_text(size = 5))
  
  # save
  ggsave(
    filename = paste0("plots/varying_effects_", effect, ".pdf"),
    plot = out,
    height = 6.5,
    width = 5
  )
  
  # clean up
  rm(fit_draws_model)
  
  # return
  out
  
}


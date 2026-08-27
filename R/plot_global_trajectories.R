#' Plot global trajectories implied by fitted model
#'
#' @param fit_draws_model Tibble of posterior draws from the model.
#' @param intervention_vars Integers, 1-9. Which variables to intervene on.
#' @param intervention_values Numeric. Which values to set during the 
#'   interventions.
#' @param intervention_times Numeric. When to implement the interventions.
#'
#' @returns A patchwork of ggplot objects
#'
plot_global_trajectories <- function(fit_draws_model,
                                     intervention_vars = NULL,
                                     intervention_values = NULL,
                                     intervention_times = NULL) {
  
  # rename draws
  draws <- fit_draws_model
  
  # get number of draws
  ndraws <- length(draws$init_logit_pop)
  
  # get times for each period (in centuries)
  times1 <- seq(-100, 0)
  times2 <- seq(1, 16)
  times3 <- seq(17, 20)
  
  # get indexes for each period
  i1 <- 1:length(times1)
  i2 <- (length(times1) + 1):length(c(times1, times2))
  i3 <- (length(c(times1, times2)) + 1):length(c(times1, times2, times3))
  
  # set up results array
  out <- array(NA, dim = c(ndraws, length(c(times1, times2, times3)), 9))
  
  # get posterior trajectories
  for (j in 1:ndraws) {
    
    # period 1
    out[j, i1, ] <-
      simulate_ODE(
        initial_values = c(
          draws$init_logit_pop[j],
          draws$init_pop_size[j],
          draws$init_logit_crop[j],
          draws$init_cropland[j],
          draws$init_logit_irr[j],
          draws$init_irrigated[j],
          draws$init_logit_urb[j],
          draws$init_urban[j],
          draws$init_gini[j]
        ),
        theta = c(
          draws$`theta[1,1]`[j],
          draws$`theta[1,2]`[j],
          draws$`theta[1,3]`[j],
          draws$`theta[1,4]`[j],
          draws$`theta[1,5]`[j],
          draws$`theta[1,6]`[j],
          draws$`theta[1,7]`[j],
          draws$`theta[1,8]`[j],
          draws$`theta[1,9]`[j],
          draws$`theta[1,10]`[j],
          draws$`theta[1,11]`[j],
          draws$`theta[1,12]`[j],
          draws$`theta[1,13]`[j],
          draws$`theta[1,14]`[j]
        ),
        times = times1,
        intervention_vars = intervention_vars,
        intervention_times = intervention_times,
        intervention_values = intervention_values
      )
    
    # period 2
    out[j, i2, ] <-
      simulate_ODE(
        initial_values = c(
          # use values at end of previous period
          out[j, length(times1), 1],
          out[j, length(times1), 2],
          out[j, length(times1), 3],
          out[j, length(times1), 4],
          out[j, length(times1), 5],
          out[j, length(times1), 6],
          out[j, length(times1), 7],
          out[j, length(times1), 8],
          out[j, length(times1), 9]
        ),
        theta = c(
          draws$`theta[2,1]`[j],
          draws$`theta[2,2]`[j],
          draws$`theta[2,3]`[j],
          draws$`theta[2,4]`[j],
          draws$`theta[2,5]`[j],
          draws$`theta[2,6]`[j],
          draws$`theta[2,7]`[j],
          draws$`theta[2,8]`[j],
          draws$`theta[1,9]`[j],
          draws$`theta[1,10]`[j],
          draws$`theta[1,11]`[j],
          draws$`theta[1,12]`[j],
          draws$`theta[1,13]`[j],
          draws$`theta[1,14]`[j]
        ),
        times = times2,
        intervention_vars = intervention_vars,
        intervention_times = intervention_times,
        intervention_values = intervention_values
      )
    
    # period 3
    out[j, i3, ] <-
      simulate_ODE(
        initial_values = c(
          # use values at end of previous period
          out[j, length(c(times1, times2)), 1],
          out[j, length(c(times1, times2)), 2],
          out[j, length(c(times1, times2)), 3],
          out[j, length(c(times1, times2)), 4],
          out[j, length(c(times1, times2)), 5],
          out[j, length(c(times1, times2)), 6],
          out[j, length(c(times1, times2)), 7],
          out[j, length(c(times1, times2)), 8],
          out[j, length(c(times1, times2)), 9]
        ),
        theta = c(
          draws$`theta[3,1]`[j],
          draws$`theta[3,2]`[j],
          draws$`theta[3,3]`[j],
          draws$`theta[3,4]`[j],
          draws$`theta[3,5]`[j],
          draws$`theta[3,6]`[j],
          draws$`theta[3,7]`[j],
          draws$`theta[3,8]`[j],
          draws$`theta[1,9]`[j],
          draws$`theta[1,10]`[j],
          draws$`theta[1,11]`[j],
          draws$`theta[1,12]`[j],
          draws$`theta[1,13]`[j],
          draws$`theta[1,14]`[j]
        ),
        times = times3,
        intervention_vars = intervention_vars,
        intervention_times = intervention_times,
        intervention_values = intervention_values
      )
    
  }
  
  # get trajectories for each variable
  P <- plogis(out[, , 1]) * exp(out[, , 2])
  C <- plogis(out[, , 3]) * exp(out[, , 4])
  I <- plogis(out[, , 5]) * exp(out[, , 6])
  U <- plogis(out[, , 7]) * exp(out[, , 8])
  G <- plogis(out[, , 9])
  
  # plot gini trajectory
  pA <-
    tibble(
      variable = "Gini",
      time = c(times1, times2, times3),
      median = apply(G, 2, function(x) median(x)),
      lower  = apply(G, 2, function(x) quantile(x, 0.025)),
      upper  = apply(G, 2, function(x) quantile(x, 0.975)),
    ) |>
    ggplot(
      mapping = aes(
        x = ((time * 100) - 2026) / 1000,
        y = median,
        ymin = lower,
        ymax = upper
      )
    ) +
    geom_ribbon(
      fill = "lightpink"
    ) +
    geom_line() +
    labs(
      x = "Time before present (ky)",
      y = "Gini"
    ) +
    scale_y_continuous(
      limits = c(0, 1)
    ) +
    theme_classic()
  
  # plot population size trajectory
  pB <-
    tibble(
      variable = "Population size",
      time = c(times1, times2, times3),
      median = apply(P, 2, function(x) median(log(x + 1))),
      lower  = apply(P, 2, function(x) quantile(log(x + 1), 0.025)),
      upper  = apply(P, 2, function(x) quantile(log(x + 1), 0.975)),
    ) |>
    ggplot(
      mapping = aes(
        x = ((time * 100) - 2026) / 1000,
        y = median,
        ymin = lower,
        ymax = upper
      )
    ) +
    geom_ribbon(
      fill = "lightblue"
    ) +
    labs(
      x = "Time before present (ky)",
      y = "Population size\n(log + 1)"
    ) +
    geom_line() +
    theme_classic()
  
  # plot cropland trajectory
  pC <-
    tibble(
      variable = "Cropland",
      time = c(times1, times2, times3),
      median = apply(C, 2, function(x) median(log(x + 1))),
      lower  = apply(C, 2, function(x) quantile(log(x + 1), 0.025)),
      upper  = apply(C, 2, function(x) quantile(log(x + 1), 0.975)),
    ) |>
    ggplot(
      mapping = aes(
        x = ((time * 100) - 2026) / 1000,
        y = median,
        ymin = lower,
        ymax = upper
      )
    ) +
    geom_ribbon(
      fill = "lightgreen"
    ) +
    labs(
      x = "Time before present (ky)",
      y = "Cropland\n(log + 1)"
    ) +
    geom_line() +
    theme_classic()
  
  # plot irrigated trajectory
  pD <-
    tibble(
      variable = "Irrigated",
      time = c(times1, times2, times3),
      median = apply(I, 2, function(x) median(log(x + 1))),
      lower  = apply(I, 2, function(x) quantile(log(x + 1), 0.025)),
      upper  = apply(I, 2, function(x) quantile(log(x + 1), 0.975)),
    ) |>
    ggplot(
      mapping = aes(
        x = ((time * 100) - 2026) / 1000,
        y = median,
        ymin = lower,
        ymax = upper
      )
    ) +
    geom_ribbon(
      fill = "#FFB09C"
    ) +
    labs(
      x = "Time before present (ky)",
      y = "Irrigated area\n(log + 1)"
    ) +
    geom_line() +
    theme_classic()
  
  # plot urban trajectory
  pE <-
    tibble(
      variable = "Urban",
      time = c(times1, times2, times3),
      median = apply(U, 2, function(x) median(log(x + 1))),
      lower  = apply(U, 2, function(x) quantile(log(x + 1), 0.025)),
      upper  = apply(U, 2, function(x) quantile(log(x + 1), 0.975)),
    ) |>
    ggplot(
      mapping = aes(
        x = ((time * 100) - 2026) / 1000,
        y = median,
        ymin = lower,
        ymax = upper
      )
    ) +
    geom_ribbon(
      fill = "#CC93CC"
    ) +
    labs(
      x = "Time before present (ky)",
      y = "Urban area\n(log + 1)"
    ) +
    geom_line() +
    theme_classic()
  
  # put together
  p <- 
    pA + pB + pC + pD + pE + 
    plot_layout(
      axes = "collect_x",
      design = '
      AA
      AA
      BC
      DE
      '
    )
  
  # save
  ggsave(
    plot = p,
    file = "plots/global_predictions.pdf",
    height = 6,
    width = 6
  )
  
  # cleanup
  rm(fit_draws_model, out, P, C, G, I, U)
  
  # return
  p
  
}

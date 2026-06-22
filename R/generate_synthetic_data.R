#' Generate synthetic data given true parameter values
#'
#' @param N_times Number of time steps
#' @param prob_site Probability of sampling a site in a time step
#' @param init_gini Initial value for Gini (on logit scale)
#' @param init_pop_size Initial value for population size (on log scale)
#' @param k Rate of mean reversion for Gini
#' @param alpha Intercept for equilibrium value for Gini
#' @param beta Effect of log population size on equilibrium value for Gini
#' @param r Rate of exponential growth for population size
#' @param phi Precision parameter for beta distribution
#' @param sigma Variance parameter for log-normal distribution
#'
#' @returns A data frame of synthetic data
#'
generate_synthetic_data <- function(N_times = 1000,
                                    prob_site = 0.3,
                                    init_gini = -1,
                                    init_pop_size = 5,
                                    k = -1,
                                    alpha = 0,
                                    beta = 1,
                                    r = 2,
                                    phi = 8,
                                    sigma = 4) {
  
  # sequence of times
  ts <- seq(0, 1, length.out = N_times)
  
  # get initial states
  latent <- matrix(NA, nrow = N_times, ncol = 2)
  latent[1, 1] <- init_gini
  latent[1, 2] <- init_pop_size
  
  # get remaining states
  for (i in 2:N_times) {
    t <- ts[i] - ts[i - 1]
    y <- latent[i - 1, 1]
    x <- latent[i - 1, 2]
    latent[i, 1] <- y + ((-exp(k) * (y - (alpha + beta * x))) * t)
    latent[i, 2] <- x + ((exp(r)) * t)
  }
  
  # simulate sites with measurement error
  d <- data.frame()
  for (i in 1:N_times) {
    while (runif(1) < prob_site) {
      mu <- plogis(latent[i, 1])
      nu <- latent[i, 2]
      d <- rbind(
        d,
        data.frame(
          t = ts[i],
          gini = rbeta(1, mu * phi, (1 - mu) * phi),
          pop_size = rlnorm(1, nu, sigma)
        )
      )
    }
  }
  
  # return simulated data
  d
  
}

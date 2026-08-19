#' Simulate ODE trajectory given initial values and parameters
#'
#' @param initial_values Numeric vector of length 7. Initial starting values for
#'   the following variables, in order: the probability of positive population
#'   size (logit scale), population size (log scale), the probability of 
#'   positive cropland (logit scale), cropland (log scale), the probability of 
#'   positive urban (logit scale), urban (log scale), and Gini values (logit 
#'   scale).
#' @param theta Numeric vector of length 11. ODE parameters governing the 
#'   trajectory of the variables. See Stan code for more details.
#' @param times Sequence of time steps to iterate over.
#' @param intervention_vars Integers, 1-7. Which variables to intervene on.
#' @param intervention_values Numeric. Which values to set during the 
#'   interventions.
#' @param intervention_times Numeric. When to implement the interventions.
#'
#' @returns Matrix of results with \code{length(times)} rows and five columns
#'
simulate_ODE <- function(initial_values, theta, times,
                         intervention_vars = NULL,
                         intervention_values = NULL,
                         intervention_times = NULL) {
  
  # check the intervention inputs
  if (!is.null(intervention_vars)) {
    if (!all(intervention_vars %in% 1:7)) {
      stop("intervention_vars must be from 1 to 7")
    }
    if (length(intervention_vars) != length(intervention_values)) {
      stop("intervention_vars and intervention_values must be the same length")
    }
    if (length(intervention_vars) != length(intervention_times)) {
      stop("intervention_vars and intervention_times must be the same length")
    }
  }
  
  # initial values
  logitP <- initial_values[1]  # logit( Pr(pop_size > 0) )
  logP   <- initial_values[2]  # log(pop_size)
  logitC <- initial_values[3]  # logit( Pr(cropland > 0) )
  logC   <- initial_values[4]  # log(cropland)
  logitU <- initial_values[5]  # logit( Pr(urban > 0) )
  logU   <- initial_values[6]  # log(urban)
  logitG <- initial_values[7]  # logit(gini)
  
  # parameters
  bP <- exp(theta[1])     # increase in probability of population > 0
  rP <- exp(theta[2])     # rate of population growth
  bC <- exp(theta[3])     # increase in probability of cropland > 0
  rC <- exp(theta[4])     # rate of cropland production
  bU <- exp(theta[5])     # increase in probability of urban > 0
  rU <- exp(theta[6])     # rate of urban production
  alpha <- theta[7]       # continuous-time intercept for gini
  betaP <- theta[8]       # effect of population size on gini
  betaC <- theta[9]       # effect of cropland on gini
  betaU <- theta[10]      # effect of urban on gini
  gamma <- exp(theta[11]) # damping parameter
  
  # get matrix of results
  out <- matrix(NA, nrow = length(times), ncol = 7)
  
  # record first time step
  out[1, ] <- c(logitP, logP, logitC, logC, logitU, logU, logitG)
  
  # iterate over remaining time steps
  for (i in 2:length(times)) {
    
    # time since previous time step
    dt <- times[i] - times[i - 1]
    
    # transformed parameter
    rG <- 
      alpha + 
      (betaP * (plogis(logitP) * log(exp(logP) + 1))) + 
      (betaC * (plogis(logitC) * log(exp(logC) + 1))) - 
      (betaU * (plogis(logitU) * log(exp(logU) + 1))) - 
      (gamma * logitG)
    
    # calculate change in variables over dt
    dlogitP <- bP * dt
    dlogP   <- rP * dt
    dlogitC <- bC * dt
    dlogC   <- rC * exp(logP - logC) * dt
    dlogitU <- bU * dt
    dlogU   <- rU * exp(logP - logU) * dt
    dlogitG <- rG * dt
    
    # update variables
    logitP <- logitP + dlogitP
    logP   <- logP   + dlogP
    logitC <- logitC + dlogitC
    logC   <- logC   + dlogC
    logitU <- logitU + dlogitU
    logU   <- logU   + dlogU
    logitG <- logitG + dlogitG
    
    # implement interventions
    if (!is.null(intervention_vars)) {
      for (j in 1:length(intervention_vars)) {
        if (times[i] >= intervention_times[j]) {
          if (intervention_vars[j] == 1) logitP <- intervention_values[j]
          if (intervention_vars[j] == 2) logP   <- intervention_values[j]
          if (intervention_vars[j] == 3) logitC <- intervention_values[j]
          if (intervention_vars[j] == 4) logC   <- intervention_values[j]
          if (intervention_vars[j] == 5) logitU <- intervention_values[j]
          if (intervention_vars[j] == 6) logU   <- intervention_values[j]
          if (intervention_vars[j] == 7) logitG <- intervention_values[j]
        }
      }
    }
    
    # record time step
    out[i, ] <- c(logitP, logP, logitC, logC, logitU, logU, logitG)
    
  }
  
  # return results
  out
  
}

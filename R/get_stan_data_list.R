#' Get list of data for Stan model
#'
#' @param data Tibble of data
#' @param prior_only Logical. If \code{TRUE}, the model ignores the likelihood
#'   and samples from the prior only.
#'
#' @returns A named list
#'
get_stan_data_list <- function(data, prior_only = FALSE) {
  
  # filter and sort data
  data <- 
    data |>
    # drop one row with missing date
    drop_na(Date) |>
    mutate(
      # normalise dates between 0 and 1
      Date = (Date - min(Date)) / (max(Date) - min(Date)),
      # scale population size
      pop_scaled = PolityPop / mean(PolityPop, na.rm = TRUE)
    ) |>
    # sort by dates
    arrange(Date)
  
  # get times between unique dates
  dates <- data$Date
  unique_dates <- unique(dates)
  ts <- c(0, diff(unique_dates))
  ts_idx <- sapply(dates, function(x) which(x == unique_dates))
  
  # get data list for stan
  list(
    N_sites     = nrow(data),
    N_times     = length(ts),
    gini        = data$Gini,
    pop_scaled  = ifelse(is.na(data$pop_scaled), -9999, data$pop_scaled),
    ts          = ts,
    ts_idx      = ts_idx,
    prior_only  = as.numeric(prior_only)
  )
  
}

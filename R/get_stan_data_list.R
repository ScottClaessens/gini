#' Get list of data for Stan model
#'
#' @param data Tibble of data
#'
#' @returns A named list
#'
get_stan_data_list <- function(data) {
  
  # filter and sort data
  data <- 
    data |>
    # drop one row with missing date
    drop_na(Date) |>
    mutate(
      # if gini lower bound < 0, set to just above zero
      Lower_B = ifelse(Lower_B < 0, 0.01, Lower_B),
      # normalise dates between 0 and 1
      Date = (Date - min(Date)) / (max(Date) - min(Date)),
      # log and standardise population size
      log_pop_std = as.numeric(scale(log(PolityPop)))
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
    gini_lower  = data$Lower_B,
    gini_upper  = data$Upper_B,
    log_pop_std = ifelse(is.na(data$log_pop_std), -9999, data$log_pop_std),
    ts          = ts,
    ts_idx      = ts_idx
  )
  
}

#' Get list of data for Stan model
#'
#' @param data Tibble of data
#'
#' @returns A named list
#'
get_stan_data_list <- function(data) {
  
  # filter data
  data <- 
    data |>
    # remove rows with NAs for start or end dates
    filter(!is.na(BeginDate) & !is.na(EndDate)) |>
    mutate(
      # if start date == end date, add small offset
      EndDate = ifelse(BeginDate == EndDate, EndDate + 0.1, EndDate),
      # if gini lower bound < 0, set to just above zero
      Lower_B = ifelse(Lower_B < 0, 0.01, Lower_B),
      # convert dates to 0-1 range
      BeginDate01 = 
        (BeginDate - min(BeginDate)) / (max(EndDate) - min(BeginDate)),
      EndDate01 = 
        (EndDate - min(BeginDate)) / (max(EndDate) - min(BeginDate))
    )
  
  # get data list for stan
  list(
    N          = nrow(data),
    gini_lower = data$Lower_B,
    gini_upper = data$Upper_B,
    time_lower = data$BeginDate01,
    time_upper = data$EndDate01
  )
  
}

# function to return data list for stan
get_stan_data_list <- function(data) {
  # filter data
  data <- 
    data |>
    # remove rows with NAs or non-positives for pop_size
    filter(!is.na(PolityPop) | PolityPop <= 0) |>
    # remove rows with NAs for start or end dates
    filter(!is.na(BeginDate) & !is.na(EndDate)) |>
    # if start date == end date, add small offset
    mutate(EndDate = ifelse(BeginDate == EndDate, EndDate + 0.1, EndDate))
  # get data list for stan
  list(
    N          = nrow(data),
    gini_obs   = data$Gini,
    gini_se    = (data$Upper_B - data$Lower_B) / (2 * 1.282), # assuming 80% CIs
    pop_size   = data$PolityPop,
    start_year = data$BeginDate,
    end_year   = data$EndDate
  )
}

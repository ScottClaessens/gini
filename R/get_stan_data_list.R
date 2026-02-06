# function to return data list for stan
get_stan_data_list <- function(data) {
  # for now, remove rows with NAs for pop_size
  data <- filter(data, !is.na(PolityPop))
  # get data list for stan
  list(
    N        = nrow(data),
    gini_obs = data$Gini,
    gini_se  = (data$Upper_B - data$Lower_B) / (2 * 1.282), # assuming 80% CIs
    pop_size = data$PolityPop
  )
}

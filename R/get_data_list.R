#' Get data list for Stan model
#'
#' @param data Dataset in tibble format
#'
#' @returns List
#'
get_data_list <- function(data) {
  list(
    N        = nrow(data),
    N_times  = length(unique(data$t)),
    gini     = data$gini,
    pop_size = data$polity_pop,
    t        = unique(data$t),
    t_idx    = sapply(data$t, function(x) which(x == unique(data$t)))
  )
}

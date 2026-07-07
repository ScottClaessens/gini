#' Get data list for Stan model
#'
#' @param data Dataset in tibble format
#'
#' @returns List
#'
get_data_list <- function(data) {
  list(
    N          = nrow(data),
    N_dates    = length(unique(data$date)),
    N_regions  = length(unique(data$subregion)),
    N_obs_pop  = sum(!is.na(data$pop_size)),
    N_obs_crop = sum(!is.na(data$cropland)),
    date       = sort(unique(data$date)) / 100,
    i0         = which(sort(unique(data$date)) == 0),
    i1600      = which(sort(unique(data$date)) == 1600),
    region     = as.numeric(factor(data$subregion)),
    pop_size   = data$pop_size[!is.na(data$pop_size)] + 0.001, # positive real
    cropland   = data$cropland[!is.na(data$cropland)] + 0.001, # positive real
    date_idx   = sapply(data$date, \(x) which(x == sort(unique(data$date)))),
    pop_idx    = which(!is.na(data$pop_size)),
    crop_idx   = which(!is.na(data$cropland))
  )
}

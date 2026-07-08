#' Get data list for Stan model
#'
#' @param data Dataset in tibble format
#'
#' @returns List
#'
get_data_list <- function(data) {
  
  # get longitude latitude coordinates for regions
  lon_lat <- 
    data |>
    mutate(subregion = factor(subregion)) |>
    group_by(subregion) |>
    summarise(
      longitude = unique(longitude),
      latitude = unique(latitude)
    )
  
  # convert degrees to radians
  xlon <- lon_lat$longitude * pi / 180
  xlat <- lon_lat$latitude * pi / 180
  
  # get x,y,z coordinates on a unit sphere
  coords <- matrix(nrow = length(xlon), ncol = 3)
  coords[, 1] <- cos(xlat) * cos(xlon)
  coords[, 2] <- cos(xlat) * sin(xlon)
  coords[, 3] <- sin(xlat)
  
  # normalise x,y,z coordinates so that maximum distance = 1
  coords <- coords / max(stats::dist(coords))
  
  # return list for stan
  list(
    N          = nrow(data),
    N_dates    = length(unique(data$date)),
    N_regions  = length(unique(data$subregion)),
    N_obs_pop  = sum(!is.na(data$pop_size)),
    N_obs_crop = sum(!is.na(data$cropland)),
    N_obs_gini = sum(!is.na(data$gini)),
    date       = sort(unique(data$date)) / 100,
    i0         = which(sort(unique(data$date)) == 0),
    i1600      = which(sort(unique(data$date)) == 1600),
    region     = as.numeric(factor(data$subregion)),
    pop_size   = data$pop_size[!is.na(data$pop_size)],
    cropland   = data$cropland[!is.na(data$cropland)],
    gini       = data$gini[!is.na(data$gini)],
    date_idx   = sapply(data$date, \(x) which(x == sort(unique(data$date)))),
    pop_idx    = which(!is.na(data$pop_size)),
    crop_idx   = which(!is.na(data$cropland)),
    gini_idx   = which(!is.na(data$gini)),
    coords     = coords
  )
  
}

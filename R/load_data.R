#' Load data on archaeological sites
#'
#' @param file_gini File path to Gini data
#' @param file_hyde File path to HYDE data
#'
#' @returns A tibble
#'
load_data <- function(file_gini, file_hyde) {
  
  # load Gini data
  gini <-
    file_gini |>
    read_csv(show_col_types = FALSE) |>
    transmute(
      id         = `...1`,
      region     = Region,
      subregion  = Subregion,
      latitude   = Latitude,
      longitude  = Longitude,
      date       = Date,
      gini       = Gini
    ) |>
    # nudge lon-lat coordinate for Kodiak Island to match hyde data
    mutate(
      longitude = ifelse(subregion == "North Pacific", -153.279, longitude),
      latitude = ifelse(subregion == "North Pacific", 57.176, latitude)
    ) |>
    # drop sites with no information on dates or longitude/latitude
    drop_na(any_of(c("date", "longitude", "latitude"))) |>
    # filter to later than 10,000 BCE
    filter(date >= -10000)
  
  # load HYDE data
  hyde <- read_csv(file_hyde, show_col_types = FALSE)
  
  # join datasets by lon/lat and closest hyde date prior to gini date
  left_join(
    x = gini, 
    y = hyde,
    by = join_by(longitude, latitude, closest(date >= date))
  ) |>
    # rename resulting date columns
    rename(
      date_gini = date.x,
      date_hyde = date.y
    ) |>
    # remove few cases with missing hyde data
    drop_na(any_of(c("pop_size", "urban", "cropland", "irrigated")))
  
}

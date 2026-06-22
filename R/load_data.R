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
      site       = Site,
      bigregion  = Bigregion,
      region     = Region,
      subregion  = Subregion,
      subarea    = Subarea,
      latitude   = Latitude,
      longitude  = Longitude,
      date       = Date,
      gini       = Gini
    ) |>
    # drop sites with no date information
    drop_na(date)
  
  # load HYDE data
  hyde <- read_csv(file_hyde, show_col_types = FALSE)
  
  # join datasets and return
  left_join(gini, hyde, by = c("site", "longitude", "latitude", "date")) |>
    # drop sites with no HYDE data
    drop_na(c(cropland, pop_size)) |>
    # normalise dates into 0-1 range
    mutate(t = (date - min(date)) / (max(date) - min(date))) |>
    arrange(t)
  
}

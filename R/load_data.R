#' Load data on archaeological sites
#'
#' @param file_gini File path to Gini data
#'
#' @returns A tibble
#'
load_data <- function(file_gini) {
  
  # load Gini data
  file_gini |>
    read_csv(show_col_types = FALSE) |>
    transmute(
      id         = `...1`,
      bigregion  = Bigregion,
      region     = Region,
      subregion  = Subregion,
      latitude   = Latitude,
      longitude  = Longitude,
      date       = Date,
      gini       = Gini
    )
  
}

#' Load data on archaeological sites
#'
#' @param data_file File path to dataset
#'
#' @returns A tibble
#'
load_data <- function(data_file) {
  data_file |>
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
      gini       = Gini,
      polity_pop = PolityPop
    ) |>
    drop_na(date) |>
    mutate(t = (date - min(date)) / (max(date) - min(date))) |>
    arrange(t)
}

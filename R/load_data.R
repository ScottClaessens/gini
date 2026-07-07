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
    # drop sites with no information on dates or longitude/latitude
    drop_na(c(date, longitude, latitude)) |>
    # nudge lon-lat coordinate for Kodiak Island to match hyde data
    mutate(
      longitude = ifelse(subregion == "North Pacific", -153.279, longitude),
      latitude = ifelse(subregion == "North Pacific", 57.176, latitude)
    ) |>
    # fix some subregions to ensure one subregion per lon-lat point
    mutate(
      subregion = paste(region, subregion),
      subregion = ifelse(id == 230, "W Europe Switzerland", subregion),
      subregion = ifelse(id == 664, "Great Britain Central Belt", subregion),
      subregion = ifelse(id == 902, "Great Britain East", subregion),
      subregion = ifelse(id == 916, "SE Europe Bulgaria", subregion)
    )
  
  # load HYDE data
  hyde <- read_csv(file_hyde, show_col_types = FALSE)
  
  # get unique lon-lat coordinates and subregions
  coords <- 
    dplyr::select(gini, longitude, latitude, subregion) |>
    unique()
  
  # summarise hyde data by subregion
  d <-
    left_join(
      hyde,
      coords, 
      by = c("longitude", "latitude")
    ) |>
    # summarise by subregion
    group_by(date, subregion) |>
    summarise(
      across(longitude:pop_size, function(x) mean(x, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    arrange(subregion, date)
  
  # add gini data and return
  bind_rows(d, dplyr::select(gini, date, subregion, gini)) |>
    fill(
      c(longitude, latitude),
      .by = "subregion",
      .direction = "updown"
    ) |>
    # retain only records after 10,000 BC
    filter(date >= -10000) |>
    # remove Rapa Nui as there is no HYDE data
    filter(subregion != "Polynesia Rapa Nui") |>
    # !!!!!! temporarily slice dataset
    slice(1:332) |>
    # arrange dataset
    arrange(subregion, date)
  
}

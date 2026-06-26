# Title: R script to generate hyde.csv
#
# Author: Scott Claessens
#
# This code loops over time slices in the HYDE dataset. For each time slice, the
# code downloads raster data from HYDE and links it to longitude-latitude 
# coordinates from the Gini dataset. The resulting data file is written to 
# hyde.csv in the data/hyde/ directory.
#
# Estimated run time on a MacBook Pro is around 20 minutes.
#

library(terra)      # v1.9-27
library(tidyverse)  # v2.0.0

# get time slices in hyde dataset up to 1980 (maximum date in gini dataset)
time_slices <- 
  c(
    seq(from = -10000, to = 0,    by = 1000),
    seq(from = 100,    to = 1700, by = 100),
    seq(from = 1710,   to = 1950, by = 10),
    seq(from = 1951,   to = 1980, by = 1)
  )

# get names of time slices
names(time_slices) <- 
  paste0(
    abs(time_slices), ifelse(time_slices >= 0, "AD", "BC")
  )

# load gini dataset for linking
gini <- read_csv("data/gini/SiteGiniLevel.csv", show_col_types = FALSE)

# drop sites without longitude-latitude coordinates
gini <- drop_na(gini, c(Longitude, Latitude))

# get unique longitude-latitude coordinates
unique_lon_lat <- unique(gini[, c("Longitude", "Latitude")])

# get lon-lat coordinates as points
points <- vect(
  unique_lon_lat,
  geom = c("Longitude", "Latitude"),
  crs = "EPSG:4326"
)

# initialise hyde dataset
out <- tibble()

# loop over time slices
for (i in 1:length(time_slices)) {
  
  # print progress message
  message("Downloading HYDE data for ", names(time_slices)[i])
  
  # initialise time slice data
  d <- bind_cols(date = time_slices[i], unique_lon_lat)
  
  # loop over different variables to extract
  for (variable in c("cropland", "popc")) {
    
    # get url of zip file for data
    url <- paste0(
      "https://geo.public.data.uu.nl/vault-hyde/",
      "hyde35_c9_apr2025%5B1749214444%5D/original/",
      "gbc2025_7apr_base/zip/",
      names(time_slices)[i],
      "_",
      ifelse(variable == "cropland", "lu", "pop"),
      ".zip"
    )
    
    # download zip file into data/hyde/ directory
    zip_file <- tempfile()
    download.file(url, zip_file, quiet = TRUE)
    
    # identify raster file
    raster_file <- paste0(
      variable,
      ifelse(variable == "cropland", "", "_"),
      names(time_slices)[i],
      ".asc"
    )
    
    # unzip into hyde directory
    unzip(
      zipfile = zip_file,
      files = raster_file,
      exdir = "data/hyde"
    )
    
    # load raster file
    r <- rast(paste0("data/hyde/", raster_file))
    
    # extract data from raster for lon-lat points
    d[[variable]] <- extract(r, points)[, 2]
    
    # clean up
    file.remove(paste0("data/hyde/", raster_file))
    
  }
  
  # bind rows
  out <- bind_rows(out, d)
  
}

# write to csv
out |>
  transmute(
    date      = date,
    longitude = Longitude,
    latitude  = Latitude,
    cropland  = round(cropland, 2),
    pop_size  = round(popc, 2)
  ) |>
  write_csv("data/hyde/hyde.csv")

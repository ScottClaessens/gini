# Title: R script to generate hyde.csv
#
# Author: Scott Claessens
#
# This code loops over time slices in the HYDE dataset. For each time slice, the
# code downloads raster data from HYDE and links it to the Gini dataset. If 
# there are no archaeological sites in the time slices, the loop skips the
# download. The resulting data file is written to hyde.csv in the data/hyde/ 
# directory.
#
# Estimated run time on a MacBook Pro is around 10 minutes.
#

library(raster)     # v1.9-27
library(tidyverse)  # v2.0.0

# get time slices in hyde dataset
time_slices <- 
  c(
    seq(from = -10000, to = 0,    by = 1000),
    seq(from = 100,    to = 1700, by = 100),
    seq(from = 1710,   to = 1950, by = 10),
    seq(from = 1951,   to = 2025, by = 1)
  )

# get names of time slices
names(time_slices) <- 
  paste0(
    abs(time_slices), ifelse(time_slices >= 0, "AD", "BC")
  )

# load gini dataset for linking
gini <- read_csv("data/gini/SiteGiniLevel.csv", show_col_types = FALSE)

# initialise hyde dataset
out <- tibble() 

# loop over time slices
for (i in 1:(length(time_slices) - 1)) {
  
  # retain archaeological sites in this time slice
  d <- filter(gini, Date >= time_slices[i] & Date < time_slices[i + 1])
  
  # if there are archaeological sites in this time slice, load hyde data
  if (nrow(d) > 0) {
    
    # print progress message
    message("Downloading HYDE data for ", names(time_slices)[i])
    
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
      
      # get lon-lat points
      points <- vect(
        d[, c("Longitude", "Latitude")],
        geom = c("Longitude", "Latitude"),
        crs = "EPSG:4326"
      )
      
      # load raster file
      r <- rast(paste0("data/hyde/", raster_file))
      
      # project to a metric CRS
      points_m <- project(points, "EPSG:3857")
      r_m <- project(r, "EPSG:3857")
      
      # create 5km buffers around site lon-lat points
      buffers <- buffer(points_m, width = 5000)
      
      # extract data from raster for these sites
      d[[variable]] <- extract(r_m, buffers, fun = mean, na.rm = TRUE)[, 2]
      
      # clean up
      file.remove(paste0("data/hyde/", raster_file))
      
    }
    
    # bind rows
    out <- bind_rows(out, d)
    
  } else {
    
    # print skip message
    message("Skipping HYDE data for ", names(time_slices)[i])
    
  }
  
}

# write to csv
out |>
  transmute(
    site      = Site,
    longitude = Longitude,
    latitude  = Latitude,
    date      = Date,
    cropland  = cropland,
    pop_size  = popc
  ) |>
  write_csv("data/hyde/hyde.csv")

options(tidyverse.quiet = TRUE)
library(targets)
library(tarchetypes)
library(tidyverse)
tar_option_set(packages = c("tidyverse"))
tar_source()

list(
  
  # ────────────────────────────────────────────
  # Load data
  # ────────────────────────────────────────────
  
  tar_target(file_gini, "data/gini/SiteGiniLevel.csv", format = "file"),
  tar_target(data, load_data(file_gini))

)

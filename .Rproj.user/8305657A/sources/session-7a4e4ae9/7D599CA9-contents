library(dplyr)
library(ebirdst)
library(fields)
library(ggplot2)
library(lubridate)
library(rnaturalearth)
library(sf)
library(terra)
library(tidyr)
library(prioritizr)
extract <- terra::extract

# Read the CSV file
bird_data <- read.csv("data/CO_regional_status_2023.csv")

# Extract the top 100 species with lowest state_trend_median values
list_species <- bird_data %>%
  # Sort by state_trend_median in ascending order (lowest first)
  arrange(desc(percent_pop_breeding)) %>%
  # Take the first 100 rows after sorting
  slice_head(n = 25) %>%
  # Select only the species_code column
  pull(common_name)


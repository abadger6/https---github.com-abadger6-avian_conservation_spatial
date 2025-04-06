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

# Display the results as a vector
print(list_species)

# download seasonal relative abundance data
ebirdst_download_status("wesmea",
                        pattern = "abundance_seasonal_mean")

# load seasonal mean relative abundance at 3km resolution
abd_seasonal <- load_raster("wesmea", 
                            product = "abundance", 
                            period = "seasonal",
                            metric = "mean",
                            resolution = "3km")

# extract just the breeding season relative abundance
abd_breeding <- abd_seasonal[["breeding"]]

# region boundary
region_boundary <- ne_states(iso_a2 = "US") |> 
  filter(name == "Colorado") |>
  st_transform(st_crs(abd_breeding)) |> 
  vect()


# Initialize an empty list to store individual species rasters
prop_pop <- list()

# Create a list to store rasters for the multi-layer rast
raster_list <- list()
species_names <- c()

# Loop through each species with error handling
for (sp in list_species) {
  tryCatch({
    # Process species
    cat("Processing species:", sp, "\n")
    
    # download seasonal abundance at 3km
    ebirdst_download_status(sp,
                          pattern = "proportion-population_seasonal_mean_3km")

    # load breeding season proportion of population
    pp <- load_raster(sp,
                    product = "proportion-population",
                    period = "seasonal") |>
      subset("breeding")
    
    # Convert to SpatRaster if it's not already
    # (load_raster might return a raster object rather than a terra object)
    if (!inherits(pp, "SpatRaster")) {
      pp <- terra::rast(pp)
    }
    
    # crop and mask to region using terra functions
    species_raster <- terra::mask(terra::crop(pp, region_boundary), region_boundary)
    
    # Store in the original list indexed by species
    prop_pop[[sp]] <- species_raster
    
    # Add to our list for multi-layer raster creation
    raster_list[[length(raster_list) + 1]] <- species_raster
    species_names <- c(species_names, sp)
    
  }, error = function(e) {
    # Log the error message and continue to the next species
    cat("Skipping species", sp, "- Error:", conditionMessage(e), "\n")
    # Continue to the next iteration of the loop
  })
}

# Create a multi-layer SpatRaster from the list of individual rasters
# Only proceed if we have at least one successful species
if (length(raster_list) > 0) {
  bird_features <- terra::rast(raster_list)
  
  # Name the layers by species
  names(bird_features) <- species_names
} else {
  cat("Warning: No species were successfully processed.\n")
}



# Create a multi-layer SpatRaster from the list of individual rasters
bird_features <- terra::rast(raster_list)

# Name the layers by species
names(bird_features) <- species_names

# Create a visual check of the study area
png("plots/first9_birds.png", width=800, height=800, res=150)
plot(bird_features[[1:9]], nr = 3, axes = FALSE)
dev.off()

#-----------------------------------------
# Now load and transform all other data to match abundance CRS
#-----------------------------------------
# 1. Import the land cost TIF file
cost_raster <- rast("data/places_fmv_pnas_dryad/1 estimates/places_fmv_all.tif")

reference_crs <- crs(bird_features)  # Get the CRS of the bird features raster

# Transform cost raster to match abundance CRS if different
if (crs(cost_raster) != reference_crs) {
  cat("Reprojecting cost raster to match abundance CRS...\n")
  cost_raster <- project(cost_raster, reference_crs)
}

# 3. Crop cost raster to study area
cost_raster_crop <- crop(cost_raster, region_boundary)

# Resample to match multi-layer dimensions
cost_final <- terra::resample(cost_raster_crop, bird_features, method = "bilinear")

# calculate budget
budget <- 10000 #terra::global(cost_raster_crop, "sum", na.rm = TRUE)[[1]] * 0.05

print(bird_features)
print(cost_final)

# create problem
p1 <-
  problem(cost_final, features = bird_features) %>%
  add_min_shortfall_objective(budget) %>%
  add_relative_targets(0.3) %>%
  add_binary_decisions() %>%
  add_default_solver(gap = 0.1, verbose = TRUE)

# print problem
print(p1)

# solve the problem
s1 <- solve(p1, force=TRUE)

# extract the objective
print(attr(s1, "objective"))

# extract time spent solving the problem
print(attr(s1, "runtime"))

# extract state message from the solver
print(attr(s1, "status"))

# plot the solution
plot(s1, main = "Solution", axes = FALSE)
# Add a basemap to the plot
library(rasterVis)

# Create a basemap using the region boundary
basemap <- terra::mask(terra::crop(cost_raster, region_boundary), region_boundary)

# Plot the solution with the basemap
png("plots/solution_with_basemap.png", width = 800, height = 800, res = 150)
levelplot(s1, margin = FALSE, main = "Solution with Basemap") +
    layer(sp.polygons(as_Spatial(region_boundary), col = "black")) +
    layer(sp.raster(basemap, col.regions = terrain.colors))
dev.off()

# calculate target coverage for the solution
p1_target_coverage <- eval_target_coverage_summary(p1, s1)
print(p1_target_coverage)
library(dplyr)
library(ebirdst)
library(rnaturalearth)
library(terra)
library(prioritizr)
library(sf)
extract <- terra::extract

# Read the CSV file
bird_data <- read.csv("../R-data/CO_regional_status_2023.csv")

# Extract the top 100 species with lowest state_trend_median values
list_species <- bird_data %>%
  # Sort by state_trend_median in ascending order (lowest first)
  arrange(desc(percent_pop_breeding)) %>%
  # Take the first 100 rows after sorting
  slice_head(n = 3) %>%
  # Select only the species_code column
  pull(common_name)

# Display the results as a vector
print(list_species)

# download seasonal relative abundance data
ebirdst_download_status("wesmea",
                        pattern = "proportion-population_seasonal_mean_3km")

# load seasonal mean relative abundance at 3km resolution
abd_breeding <- load_raster("wesmea",
                    product = "proportion-population",
                    period = "seasonal") |>
                    subset("breeding")

# region boundary
region_boundary <- ne_states(iso_a2 = "US") |> 
  filter(name == "Colorado")

# project boundary to match raster data
region_boundary_proj <- st_transform(region_boundary, st_crs(abd_breeding))

# project boundary to match raster data
region_boundary_bird_crs <- st_transform(region_boundary, st_crs(abd_breeding))

#-----------------------------------------
# Now load bird data
#-----------------------------------------

# find the centroid of the region
region_centroid <- region_boundary |> 
  st_geometry() |> 
  st_transform(crs = 4326) |> 
  st_centroid() |> 
  st_coordinates() |> 
  round(1)

# define projection
crs_laea <- paste0("+proj=laea +lat_0=", region_centroid[2],
                   " +lon_0=", region_centroid[1])

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

    #pp <- terra::mask(terra::crop(pp, region_boundary_bird_crs), region_boundary_bird_crs)
    #pp <- terra::crop(pp, region_boundary_bird_crs)
    pp <- crop(pp, region_boundary_proj) |> 
      mask(region_boundary_proj)
    # transform to the custom projection using nearest neighbor resampling
    species_raster <- project(pp, crs_laea, method = "near") |> 
      trim() 
      # remove areas of the raster containing no data
      

    # Convert to SpatRaster if it's not already
    # (load_raster might return a raster object rather than a terra object)
    if (!inherits(pp, "SpatRaster")) {
      pp <- terra::rast(pp)
    }
    
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
print(bird_features)

#print(bird_features)
# Create a visual check of the study area
png("plots/first9_birds.png", width=800, height=800, res=150)
plot(bird_features[[1:9]], nr = 3, axes = FALSE)
dev.off()
print("Bird features processed!")

ext_bird <- ext(bird_features)

#-----------------------------------------
# Now load and transform all other data to match abundance CRS
#-----------------------------------------
# 1. Import the land cost TIF file
#cost_raster <- terra::rast("../R-data/places_fmv_pnas_dryad/1 estimates/output.tif")#, ext=ext_bird)
cost_raster <- terra::rast("../R-data/places_fmv_pnas_dryad/1 estimates/test.tif")#, ext=ext_bird)
#https://www.pnas.org/doi/10.1073/pnas.2012865117
#print(ext(cost_raster))
#print(ext_bird)
#print(ext(cost_raster))

png("plots/base_cost.png", width=800, height=800, res=150)
plot(cost_raster, axes = FALSE)
dev.off()

cost_raster <- project(cost_raster, crs_laea, res=terra::res(bird_features)[1], method = "near") |> 
      trim() 

# Resample to match multi-layer dimensions
cost_final <- terra::resample(cost_raster, bird_features, method = "bilinear")

# calculate budget
budget <- 10000 #terra::global(cost_raster_crop, "sum", na.rm = TRUE)[[1]] * 0.05

png("plots/final_cost.png", width=800, height=800, res=150)
plot(cost_final, axes = FALSE)
dev.off()

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

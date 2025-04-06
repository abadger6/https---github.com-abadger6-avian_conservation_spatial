# Load required libraries
library(ebirdst)
library(prioritizr)
library(sf)
library(terra)
library(tidyverse)

#-----------------------------------------
# Define a square study area by providing:
# - center_x: longitude of center point
# - center_y: latitude of center point
# - width_km: width of square in kilometers
#-----------------------------------------
define_study_area <- function(center_x, center_y, width_km) {
  # Convert width from km to degrees (approximate)
  # 1 degree latitude ≈ 111 km
  # 1 degree longitude varies with latitude
  lat_degree_km <- 111
  long_degree_km <- 111 * cos(center_y * pi/180)
  
  width_lat <- width_km / lat_degree_km
  width_lon <- width_km / long_degree_km
  
  # Create a square
  xmin <- center_x - width_lon/2
  xmax <- center_x + width_lon/2
  ymin <- center_y - width_lat/2
  ymax <- center_y + width_lat/2
  
  # Create a polygon
  square <- st_polygon(list(rbind(
    c(xmin, ymin),
    c(xmax, ymin),
    c(xmax, ymax),
    c(xmin, ymax),
    c(xmin, ymin)
  )))
  
  # Convert to sf object with CRS
  square_sf <- st_sfc(square, crs = 4326)  # WGS84
  return(square_sf)
}

#-----------------------------------------
# USER INPUT: Define your study area here
#-----------------------------------------
study_center_x <- -93.2650  # Minneapolis longitude
study_center_y <- 44.9778   # Minneapolis latitude
study_width_km <- 500       # 500km square

# Create the study area boundary (initially in WGS84)
study_area <- define_study_area(study_center_x, study_center_y, study_width_km)

#-----------------------------------------
# FIRST: Get the reference CRS from abundance layer
#-----------------------------------------
# Download species data first (before loading other data)
species_code <- "grbher3" # Great Blue Heron
path <- ebirdst_download_status(species = species_code)

# Load abundance data to get reference CRS
abundance <- load_raster(species_code, resolution = "27km")
reference_crs <- crs(abundance)
cat("Using reference CRS from abundance layer:\n")
print(reference_crs)

# Create a visual check of the study area
png("plots/study_area_map.png", width=800, height=800, res=150)
plot(st_geometry(study_area), main="Study Area", col="transparent", border="red", lwd=2)
dev.off()

#-----------------------------------------
# Now load and transform all other data to match abundance CRS
#-----------------------------------------
# 1. Import the land cost TIF file
cost_raster <- rast("data/places_fmv_pnas_dryad/1 estimates/places_fmv_all.tif")

# Transform cost raster to match abundance CRS if different
if (crs(cost_raster) != reference_crs) {
  cat("Reprojecting cost raster to match abundance CRS...\n")
  cost_raster <- project(cost_raster, reference_crs)
}

# 2. Transform study area to match abundance CRS
study_area_transformed <- st_transform(study_area, reference_crs)

# 3. Crop cost raster to study area
cost_raster_crop <- crop(cost_raster, vect(study_area_transformed))

# Plot the cropped cost raster
png("plots/cropped_land_cost_map.png", width=1200, height=800, res=300)
plot(cost_raster_crop, col=hcl.colors(100, "YlOrRd"), 
     main="Land Costs - Study Area")
# Add study area boundary
plot(vect(study_area_transformed), add=TRUE, border="black", lwd=2)
dev.off()

# 4. Crop abundance data to study area
abundance_crop <- crop(abundance, vect(study_area_transformed))
#print(crs(abundance_crop)) # Verify CRS

# 5. Define planning units
planning_units <- st_read("data/cb_2018_us_county_500k/cb_2018_us_county_500k.shp")

# Transform planning units to match abundance CRS
planning_units <- st_transform(planning_units, reference_crs)

# Find which planning units intersect with our study area
planning_units_filtered <- planning_units[st_intersects(planning_units, study_area_transformed, sparse=FALSE), ]

# Convert filtered planning units to SpatVector
planning_units_vect <- vect(planning_units_filtered)

# Plot the filtered planning units
png("plots/filtered_planning_units.png", width=1200, height=800, res=300)
plot(planning_units_vect, main="Planning Units in Study Area")
plot(vect(study_area_transformed), add=TRUE, border="red", lwd=2)
dev.off()

# Check spatial alignment - all should now be in the same CRS
png("plots/debug_spatial_alignment.png", width=1200, height=800, res=150)
par(mfrow=c(1,3))
plot(abundance_crop, main="Abundance")
plot(cost_raster_crop, main="Cost")
plot(planning_units_vect, main="Planning Units")
dev.off()

# Verify all layers have the same CRS
cat("\nVerifying all layers have the same CRS:\n")
cat("Abundance CRS matches reference CRS:", crs(abundance_crop) == reference_crs, "\n")
cat("Cost CRS matches reference CRS:", crs(cost_raster_crop) == reference_crs, "\n")
cat("Planning units CRS matches reference CRS:", crs(planning_units_vect) == reference_crs, "\n")
#-----------------------------------------
# Extract values for planning units - FIXED VERSION
#-----------------------------------------
# 6. Extract abundance values for each planning unit
abundance_values <- terra::extract(abundance_crop, planning_units_vect)

# Print diagnostic information for debugging
cat("\nDiagnostic Information:\n")
cat("Number of extracted abundance values:", nrow(abundance_values), "\n")
cat("Columns in abundance_values:", paste(colnames(abundance_values), collapse=", "), "\n")

# Add error checking
if(nrow(abundance_values) == 0) {
  stop("No abundance values were extracted. Check if the rasters and planning units overlap.")
}

# 7. Calculate sum for each planning unit with direct column references
# First identify the value column name (excluding the ID column)
value_col_name <- setdiff(colnames(abundance_values), "ID")[1]
cat("Using value column:", value_col_name, "\n")

# Use tapply for aggregation instead of formula-based aggregate
# This handles multiple values per ID properly
abundance_sums <- tapply(abundance_values[[value_col_name]], 
                        abundance_values$ID, 
                        FUN = function(x) sum(x, na.rm = TRUE))

# Convert to data frame for easier handling
abundance_summary <- data.frame(
  ID = as.integer(names(abundance_sums)),
  abundance = as.numeric(abundance_sums)
)

cat("Aggregation result dimensions:", nrow(abundance_summary), "x", ncol(abundance_summary), "\n")

# 8. Add abundance to planning units with explicit checks
# First create a data frame with all planning unit IDs
pu_ids <- data.frame(ID = 1:nrow(planning_units_filtered))

# Join using merge with explicit column matching
planning_units_data <- merge(pu_ids, abundance_summary, by = "ID", all.x = TRUE)

# Now assign to planning units safely
planning_units_filtered$abundance <- planning_units_data$abundance
planning_units_filtered$abundance[is.na(planning_units_filtered$abundance)] <- 0

# 9. Extract cost values for each planning unit - using the same approach
cost_values <- terra::extract(cost_raster_crop, planning_units_vect)

if(nrow(cost_values) == 0) {
  stop("No cost values were extracted. Check if the cost raster and planning units overlap.")
}

# Find cost column name
cost_col_name <- setdiff(colnames(cost_values), "ID")[1]
cat("Using cost column:", cost_col_name, "\n")

# Use tapply for cost aggregation (using mean instead of sum)
cost_means <- tapply(cost_values[[cost_col_name]], 
                    cost_values$ID, 
                    FUN = function(x) mean(x, na.rm = TRUE))

# Convert to data frame
cost_summary <- data.frame(
  ID = as.integer(names(cost_means)),
  cost = as.numeric(cost_means)
)

# Join with planning units
planning_units_data <- merge(pu_ids, cost_summary, by = "ID", all.x = TRUE)
planning_units_filtered$cost <- planning_units_data$cost
planning_units_filtered$cost[is.na(planning_units_filtered$cost)] <- 0

print(planning_units_filtered[, c("ID", "abundance", "cost")])  # Print summary of abundance and cost

# Create a plot to check the spatial overlap between abundance and planning units
png("plots/abundance_overlap_check.png", width=1200, height=800, res=150)
plot(abundance_crop, main="Abundance with Planning Units")
plot(planning_units_vect, add=TRUE, border="red")
dev.off()

# Check if abundance raster has any data in the study area
has_abundance_data <- !all(is.na(values(abundance_crop)))
cat("Abundance raster has data in study area:", has_abundance_data, "\n")

# Visualize the planning units that have non-zero abundance values
png("plots/abundance_by_county.png", width=1200, height=800, res=150)
planning_units_filtered$has_data <- planning_units_filtered$abundance > 0
plot(planning_units_filtered["has_data"], main="Counties with Abundance Data")
dev.off()

# Create a comparison plot of uncropped vs cropped abundance data
png("plots/abundance_comparison.png", width=1600, height=800, res=150)

# Set up a 1x2 plotting grid
par(mfrow=c(1,2))

# Plot 1: Original uncropped abundance raster
plot(abundance, 
     main="Original Abundance Data",
     col=hcl.colors(100, "viridis"),
     legend=FALSE)  # Hide legend from first plot

# Add study area boundary to uncropped plot
plot(vect(study_area_transformed), add=TRUE, border="red", lwd=2)

# Plot 2: Cropped abundance raster
plot(abundance_crop, 
     main="Cropped Abundance Data",
     col=hcl.colors(100, "viridis"))

# Add study area boundary to cropped plot
plot(vect(study_area_transformed), add=TRUE, border="red", lwd=2)

# Add planning units outline to see overlap
plot(planning_units_vect, add=TRUE, border="white", lwd=0.5)

# Add a common legend
par(mfrow=c(1,1))
# Reset to single plot layout

# Create a more detailed diagnostic visualization
dev.off()

# Create a detailed abundance map with planning units
png("plots/abundance_detail.png", width=1200, height=800, res=200)
plot(abundance_crop, 
     main="Great Blue Heron Abundance with Planning Units",
     col=hcl.colors(100, "viridis"))
plot(planning_units_vect, add=TRUE, border="white", lwd=0.8)
plot(vect(study_area_transformed), add=TRUE, border="red", lwd=2.5)

# Add a small map showing non-NA cells
inset_data <- !is.na(abundance_crop)
par(fig=c(0.02, 0.25, 0.7, 0.95), new=TRUE)
plot(inset_data, 
     legend=FALSE, 
     main="Data Coverage",
     col=c("transparent", "darkgreen"))
plot(vect(study_area_transformed), add=TRUE, border="red", lwd=1.5)
dev.off()

#-----------------------------------------
# Create and solve the conservation problem
#-----------------------------------------
# 12. Calculate budget (5% of total cost)
budget <- sum(planning_units_filtered$cost, na.rm=TRUE) * 0.05

# 13. Create the conservation problem
prob <- problem(planning_units_filtered, features = "abundance", cost_column = "cost") %>%
  add_min_shortfall_objective(budget) %>%
  add_relative_targets(0.3) %>%
  add_binary_decisions() %>%
  add_default_solver(verbose = TRUE)

# 14. Solve the problem
solution <- solve(prob)

# 15. Visualize results
png("plots/solution_map.png", width=1200, height=800, res=300)
plot(solution["solution_1"], main="Conservation Solution")
plot(vect(study_area_transformed), add=TRUE, border="red", lwd=2)
dev.off()

# 16. Save the solution for further analysis
st_write(solution, "output/conservation_solution.shp", append=FALSE)

# 17. Print summary of solution
selected_units <- sum(solution$solution_1)
total_cost <- sum(solution$solution_1 * solution$cost)
total_abundance_protected <- sum(solution$solution_1 * solution$abundance)
percent_abundance_protected <- total_abundance_protected / sum(solution$abundance) * 100

cat("\n--- SOLUTION SUMMARY ---\n")
cat("Planning units selected:", selected_units, "out of", nrow(solution), "\n")
cat("Total cost:", total_cost, "\n")
cat("Total abundance protected:", total_abundance_protected, "\n")
cat("Percent of abundance protected:", percent_abundance_protected, "%\n")
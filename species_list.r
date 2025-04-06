# Load required packages
library(ebirdst)
library(rebird)
library(dplyr)

#' Get top 100 most common bird species abbreviations for a state
#' 
#' @param state_code Two-letter state code (e.g., "CA" for California)
#' @return Character vector of species codes (abbreviations) for the most common species
get_top_species_abbrev <- function(state_code) {
  # Ensure state_code is uppercase
  state_code <- toupper(state_code)
  
  # Create the region code for the eBird API
  region_code <- paste0("US-", state_code)
  
  # Get recent observations from the state (past year)
  # This provides a good proxy for common species
  recent_obs <- ebirdregion(region = region_code, back = 365, max = 10000)
  
  # Count occurrences of each species and sort by frequency
  species_counts <- recent_obs %>%
    count(speciesCode) %>%
    arrange(desc(n))
  
  # Get species available in the eBird Status and Trends dataset
  st_species <- ebirdst_runs()
  
  # Find intersection of observed species and those in Status and Trends
  # This ensures we only work with species that have complete data
  common_species <- species_counts %>%
    filter(speciesCode %in% st_species$species_code) %>%
    head(100)
  
  # Return the species codes (abbreviations)
  return(common_species$speciesCode)
}

# Example usage
ca_top_species <- get_top_species_abbrev("CA")
print(ca_top_species)
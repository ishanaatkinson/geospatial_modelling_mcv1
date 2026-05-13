# ================================
# Libraries
# ================================
library(dplyr)
library(purrr)
library(terra)
library(sf)
library(haven)

# ================================
# Load DHS data
# ================================
all_KR_surveys_with_gps_HR_IR <- readRDS(
  "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/all_KR_surveys_with_gps_HR_IR.rds"
)

# ================================
# Load rasters
# ================================
travel_time_to_city_folder <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/MAP distance times to HF/2015_accessibility_to_cities_v1.0"
motorized_travel_time_to_HC_folder <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/MAP distance times to HF/2020_motorized_travel_time_to_healthcare"
walking_travel_time_to_HC_folder <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/MAP distance times to HF/2020_walking_only_travel_time_to_healthcare"

city_raster  <- rast(file.path(travel_time_to_city_folder, "2015_accessibility_to_cities_v1.0.tif"))
motor_raster <- rast(file.path(motorized_travel_time_to_HC_folder, "2020_motorized_travel_time_to_healthcare.geotiff"))
walk_raster  <- rast(file.path(walking_travel_time_to_HC_folder, "2020_walking_only_travel_time_to_healthcare.geotiff"))

# ================================
# Initialize output list
# ================================
all_KR_surveys_with_gps_HR_IR_maps <- list()

# ================================
# Loop through surveys
# ================================
for(survey_name in names(all_KR_surveys_with_gps_HR_IR)) {
  
  message("Processing: ", survey_name)
  
  df <- all_KR_surveys_with_gps_HR_IR[[survey_name]]
  
  # -------------------------------
  # Clean labelled variables
  # -------------------------------
  df <- df %>%
    mutate(across(where(haven::is.labelled), haven::zap_labels))
  
  # -------------------------------
  # Create cluster-level coordinates
  # -------------------------------
  cluster_df <- df %>%
    mutate(cluster = if("v001" %in% names(.)) v001 else hv001) %>%
    select(cluster, LATNUM, LONGNUM) %>%
    rename(lat = LATNUM, lon = LONGNUM) %>%
    filter(!is.na(lat) & !is.na(lon)) %>%
    distinct(cluster, lat, lon)
  
  # -------------------------------
  # Convert to spatial
  # -------------------------------
  points_sf  <- st_as_sf(cluster_df, coords = c("lon", "lat"), crs = 4326)
  points_vect <- vect(points_sf)
  
  # -------------------------------
  # Extract travel times
  # -------------------------------
  city_time  <- terra::extract(city_raster, points_vect)[,2]
  motor_time <- terra::extract(motor_raster, points_vect)[,2]
  walk_time  <- terra::extract(walk_raster, points_vect)[,2]
  
  # -------------------------------
  # Attach travel times to clusters
  # -------------------------------
  cluster_df <- cluster_df %>%
    mutate(
      travel_time_to_city = city_time,
      travel_time_to_HC_motor = motor_time,
      travel_time_to_HC_walk = walk_time
    )
  
  # -------------------------------
  # Merge back to full dataset
  # -------------------------------
  df <- df %>%
    mutate(cluster = if("v001" %in% names(.)) v001 else hv001) %>%
    left_join(cluster_df, by = c("cluster", "LATNUM" = "lat", "LONGNUM" = "lon"))
  
  # -------------------------------
  # Save to list
  # -------------------------------
  all_KR_surveys_with_gps_HR_IR_maps[[survey_name]] <- df
}

# ================================
# Save final dataset
# ================================
saveRDS(
  all_KR_surveys_with_gps_HR_IR_maps,
  "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/all_KR_surveys_with_gps_HR_IR_maps.rds"
)

cat("✅ All surveys processed and saved.\n")


library(dplyr)
library(purrr)
library(tibble)

# Define the variables you want to check
travel_vars <- c(
  "travel_time_to_city",
  "travel_time_to_HC_motor",
  "travel_time_to_HC_walk"
)

# Loop through all surveys
missing_summary <- map_dfr(names(all_KR_surveys_with_gps_HR_IR_maps), function(survey_name) {
  
  df <- all_KR_surveys_with_gps_HR_IR_maps[[survey_name]]
  
  # For each variable, check if it exists
  map_dfr(travel_vars, function(var) {
    
    if(!var %in% names(df)) {
      tibble(
        survey = survey_name,
        variable = var,
        exists_in_data = FALSE,
        missing_count = NA_integer_,
        total_rows = nrow(df),
        missing_pct = NA_real_
      )
    } else {
      missing_count <- sum(is.na(df[[var]]))
      tibble(
        survey = survey_name,
        variable = var,
        exists_in_data = TRUE,
        missing_count = missing_count,
        total_rows = nrow(df),
        missing_pct = 100 * missing_count / nrow(df)
      )
    }
  })
})

# View summary
missing_summary

# container for results
pipeline_checks_list <- list()

############################################
# 1️⃣ LOAD KR ONLY
############################################

all_KR_surveys <- readRDS("C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/all_KR_surveys.rds")

for(survey_name in names(all_KR_surveys)){
  
  df <- haven::zap_labels(all_KR_surveys[[survey_name]])
  
  year_kr <- if("SurveyYear" %in% names(df)) unique(df$SurveyYear)
  else if("v007" %in% names(df)) unique(df$v007)
  else NA
  
  hh_kr <- if(all(c("v001","v002") %in% names(df))) {
    nrow(unique(df[c("v001","v002")]))
  } else NA
  
  pipeline_checks_list[[survey_name]] <- data.frame(
    
    survey_name = survey_name,
    year_KR = paste(year_kr, collapse = ","),
    hh_count_KR = hh_kr,
    
    rows_KR = nrow(df),
    cols_KR = ncol(df),
    
    stringsAsFactors = FALSE
  )
}

rm(all_KR_surveys)
# gc()

############################################
# 2️⃣ LOAD KR + GPS
############################################

all_KR_surveys_with_gps <- readRDS("C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/all_KR_surveys_with_gps.rds")

for(survey_name in names(all_KR_surveys_with_gps)){
  
  df <- haven::zap_labels(all_KR_surveys_with_gps[[survey_name]])
  
  year_gps <- if("SurveyYear" %in% names(df)) unique(df$SurveyYear)
  else if("v007" %in% names(df)) unique(df$v007)
  else NA
  
  hh_gps <- if(all(c("v001","v002") %in% names(df))) {
    nrow(unique(df[c("v001","v002")]))
  } else NA
  
  pipeline_checks_list[[survey_name]]$year_gps <- paste(year_gps, collapse = ",")
  pipeline_checks_list[[survey_name]]$hh_count_gps <- hh_gps
  
  pipeline_checks_list[[survey_name]]$rows_gps <- nrow(df)
  pipeline_checks_list[[survey_name]]$cols_gps <- ncol(df)
}

rm(all_KR_surveys_with_gps)
# gc()

############################################
# 3️⃣ LOAD KR + GPS + HR
############################################

all_KR_surveys_with_gps_HR <- readRDS("C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/all_KR_surveys_with_gps_HR.rds")

for(survey_name in names(all_KR_surveys_with_gps_HR)){
  
  df <- haven::zap_labels(all_KR_surveys_with_gps_HR[[survey_name]])
  
  year_HR <- if("SurveyYear" %in% names(df)) unique(df$SurveyYear)
  else if("v007" %in% names(df)) unique(df$v007)
  else NA
  
  hh_HR <- if(all(c("v001","v002") %in% names(df))) {
    nrow(unique(df[c("v001","v002")]))
  } else NA
  
  pipeline_checks_list[[survey_name]]$year_HR <- paste(year_HR, collapse = ",")
  pipeline_checks_list[[survey_name]]$hh_count_HR <- hh_HR
  
  pipeline_checks_list[[survey_name]]$rows_HR <- nrow(df)
  pipeline_checks_list[[survey_name]]$cols_HR <- ncol(df)
}

rm(all_KR_surveys_with_gps_HR)
# gc()

############################################
# 4️⃣ LOAD KR + GPS + HR + IR
############################################

all_KR_surveys_with_gps_HR_IR <- readRDS("C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/all_KR_surveys_with_gps_HR_IR.rds")

for(survey_name in names(all_KR_surveys_with_gps_HR_IR)){
  
  df <- haven::zap_labels(all_KR_surveys_with_gps_HR_IR[[survey_name]])
  
  year_IR <- if("SurveyYear" %in% names(df)) unique(df$SurveyYear)
  else if("v007" %in% names(df)) unique(df$v007)
  else NA
  
  hh_IR <- if(all(c("v001","v002") %in% names(df))) {
    nrow(unique(df[c("v001","v002")]))
  } else NA
  
  pipeline_checks_list[[survey_name]]$year_IR <- paste(year_IR, collapse = ",")
  pipeline_checks_list[[survey_name]]$hh_count_IR <- hh_IR
  
  pipeline_checks_list[[survey_name]]$rows_IR <- nrow(df)
  pipeline_checks_list[[survey_name]]$cols_IR <- ncol(df)
}

rm(all_KR_surveys_with_gps_HR_IR)
# gc()

############################################
# 5️⃣ LOAD KR + GPS + HR + IR + MAPS
############################################

all_KR_surveys_with_gps_HR_IR_maps <- readRDS("C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/all_KR_surveys_with_gps_HR_IR_maps.rds")

for(survey_name in names(all_KR_surveys_with_gps_HR_IR_maps)){
  
  df <- haven::zap_labels(all_KR_surveys_with_gps_HR_IR_maps[[survey_name]])
  
  year_maps <- if("SurveyYear" %in% names(df)) unique(df$SurveyYear)
  else if("v007" %in% names(df)) unique(df$v007)
  else NA
  
  hh_maps <- if(all(c("v001","v002") %in% names(df))) {
    nrow(unique(df[c("v001","v002")]))
  } else NA
  
  pipeline_checks_list[[survey_name]]$year_maps <- paste(year_maps, collapse = ",")
  pipeline_checks_list[[survey_name]]$hh_count_maps <- hh_maps
  
  pipeline_checks_list[[survey_name]]$rows_maps <- nrow(df)
  pipeline_checks_list[[survey_name]]$cols_maps <- ncol(df)
  
  # Check that the 3 travel time variables were added
  travel_vars <- c("travel_time_to_city", "travel_time_to_HC_motor", "travel_time_to_HC_walk")
  has_travel  <- travel_vars %in% names(df)
  pipeline_checks_list[[survey_name]]$maps_vars_present <- paste(travel_vars[has_travel], collapse = ",")
  pipeline_checks_list[[survey_name]]$maps_vars_missing <- paste(travel_vars[!has_travel], collapse = ",")
  
  # Non-NA coverage for first travel var found
  if(any(has_travel)){
    first_var <- travel_vars[has_travel][1]
    pipeline_checks_list[[survey_name]]$maps_pct_nonNA <- round(100 * sum(!is.na(df[[first_var]])) / nrow(df), 1)
  } else {
    pipeline_checks_list[[survey_name]]$maps_pct_nonNA <- NA_real_
  }
}

rm(all_KR_surveys_with_gps_HR_IR_maps)
# gc()

############################################
# 6️⃣ LOAD KR + GPS + HR + IR + MAPS + GC
############################################

all_KR_surveys_with_GC <- readRDS("C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/all_KR_surveys_with_gps_HR_IR_maps_GC.rds")

for(survey_name in names(all_KR_surveys_with_GC)){
  
  df <- haven::zap_labels(all_KR_surveys_with_GC[[survey_name]])
  
  year_GC <- if("SurveyYear" %in% names(df)) unique(df$SurveyYear)
  else if("v007" %in% names(df)) unique(df$v007)
  else NA
  
  hh_GC <- if(all(c("v001","v002") %in% names(df))) {
    nrow(unique(df[c("v001","v002")]))
  } else NA
  
  pipeline_checks_list[[survey_name]]$year_GC <- paste(year_GC, collapse = ",")
  pipeline_checks_list[[survey_name]]$hh_count_GC <- hh_GC
  
  pipeline_checks_list[[survey_name]]$rows_GC <- nrow(df)
  pipeline_checks_list[[survey_name]]$cols_GC <- ncol(df)
}

rm(all_KR_surveys_with_GC)
# gc()

############################################
# FINAL TABLE
############################################

pipeline_checks_wide <- bind_rows(pipeline_checks_list)

# ── Reorder columns: group by metric for easier cross-stage comparison ──
pipeline_checks_wide <- pipeline_checks_wide %>%
  select(
    survey_name,
    # years together
    year_KR, year_gps, year_HR, year_IR, year_maps, year_GC,
    # households together
    hh_count_KR, hh_count_gps, hh_count_HR, hh_count_IR, hh_count_maps, hh_count_GC,
    # rows together
    rows_KR, rows_gps, rows_HR, rows_IR, rows_maps, rows_GC,
    # cols together
    cols_KR, cols_gps, cols_HR, cols_IR, cols_maps, cols_GC,
    
    # anything else that may exist
    everything()
  ) %>% 
  select(-c(maps_vars_present, maps_vars_missing, maps_pct_nonNA))

# ── Year alignment (all 6 stages) ──
pipeline_checks_wide$year_alignment <-
  ifelse(
    pipeline_checks_wide$year_KR   == pipeline_checks_wide$year_gps &
      pipeline_checks_wide$year_gps  == pipeline_checks_wide$year_HR &
      pipeline_checks_wide$year_HR   == pipeline_checks_wide$year_IR &
      pipeline_checks_wide$year_IR   == pipeline_checks_wide$year_maps &
      pipeline_checks_wide$year_maps == pipeline_checks_wide$year_GC,
    "OK", "Mismatch"
  )

# ── Household alignment (all 6 stages) ──
pipeline_checks_wide$hh_alignment <-
  ifelse(
    pipeline_checks_wide$hh_count_KR   == pipeline_checks_wide$hh_count_gps &
      pipeline_checks_wide$hh_count_gps  == pipeline_checks_wide$hh_count_HR &
      pipeline_checks_wide$hh_count_HR   == pipeline_checks_wide$hh_count_IR &
      pipeline_checks_wide$hh_count_IR   == pipeline_checks_wide$hh_count_maps &
      pipeline_checks_wide$hh_count_maps == pipeline_checks_wide$hh_count_GC,
    "OK", "Mismatch"
  )

# ── Row alignment (all 6 stages) ──
pipeline_checks_wide$rows_alignment <-
  ifelse(
    pipeline_checks_wide$rows_KR   == pipeline_checks_wide$rows_gps &
      pipeline_checks_wide$rows_KR   == pipeline_checks_wide$rows_HR &
      pipeline_checks_wide$rows_KR   == pipeline_checks_wide$rows_IR &
      pipeline_checks_wide$rows_KR   == pipeline_checks_wide$rows_maps &
      pipeline_checks_wide$rows_KR   == pipeline_checks_wide$rows_GC,
    "OK", "Mismatch"
  )

# ── Column alignment (cols should increase monotonically through pipeline) ──


writexl::write_xlsx(
  pipeline_checks_wide,
  "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/pipeline_checks_wide.xlsx"
)

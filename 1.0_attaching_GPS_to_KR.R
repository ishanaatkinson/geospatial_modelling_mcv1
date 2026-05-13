library(dplyr)
library(sf)
library(tools)
library(rdhs)
library(haven)
library(writexl)

# ================================================================
# CONFIGURATION
# ================================================================

zip_folder <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS (download from site, zip only)"

input_rds   <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/all_KR_surveys.rds"
output_rds  <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/all_KR_surveys_with_gps.rds"
output_log  <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/KR_gps_log.xlsx"


# ================================================================
# 1. Load the already-processed KR surveys
# ================================================================

message("Loading KR surveys from RDS...")
all_KR_surveys <- readRDS(input_rds)
cat("KR surveys loaded from RDS:", length(all_KR_surveys), "\n")


# ================================================================
# 2. Get metadata for KR and GE datasets from rdhs
#    Mirror the approach in the KR load script: use all_datasets
#    filtered by file_type, then join on SurveyId
# ================================================================

message("Fetching dataset metadata from rdhs...")

all_datasets <- dhs_datasets(fileFormat = "FL") %>%
  dplyr::select(SurveyId, CountryName, SurveyYear,
                DHS_CountryCode, FileName) %>%
  mutate(
    file_type   = substr(FileName, 3, 4),
    survey_name = tools::file_path_sans_ext(FileName),
    SurveyYear  = as.integer(SurveyYear)
  ) %>% 
  mutate(
    FileName = toupper(FileName),
    survey_name = toupper(survey_name),
    file_type = toupper(file_type)
  )

kr_meta <- all_datasets %>%
  filter(file_type == "KR") %>%
  dplyr::rename(KR_file = FileName) %>%
  dplyr::select(SurveyId, KR_file, CountryName,
                SurveyYear, DHS_CountryCode, survey_name)

ge_meta <- all_datasets %>%
  filter(file_type == "GE") %>%
  dplyr::rename(GE_file = FileName) %>%
  dplyr::select(SurveyId, GE_file)

cat("KR surveys in rdhs metadata:", nrow(kr_meta), "\n")
cat("GE files in rdhs metadata:  ", nrow(ge_meta),  "\n")


# ================================================================
# 3. Find GE zip files on disk
#    Mirror the KR load script disk-matching logic exactly
# ================================================================

zip_files    <- list.files(zip_folder, pattern = "\\.zip$",
                           full.names = TRUE, ignore.case = TRUE)
GE_zip_files <- zip_files[grepl("GE", basename(zip_files),
                                ignore.case = TRUE)]

cat("GE zip files found on disk:", length(GE_zip_files), "\n")

# Build disk lookup: survey_name → full path (same pattern as KR load script)
disk_lookup_ge <- data.frame(
  zip_path    = GE_zip_files,
  survey_name = tools::file_path_sans_ext(basename(GE_zip_files)),
  stringsAsFactors = FALSE
)

# Join GE metadata to its disk path via survey_name
ge_on_disk <- ge_meta %>%
  mutate(survey_name = tools::file_path_sans_ext(GE_file)) %>%
  left_join(disk_lookup_ge, by = "survey_name") %>%
  mutate(
    ge_on_disk = !is.na(zip_path) & file.exists(coalesce(zip_path, ""))
  )

# Join KR metadata to GE info via SurveyId — one row per KR survey
kr_ge_lookup <- kr_meta %>%
  left_join(
    ge_on_disk %>% dplyr::select(SurveyId, GE_file, zip_path, ge_on_disk),
    by = "SurveyId"
  ) %>%
  mutate(
    # KR must also be in the loaded RDS to be processed
    kr_loaded   = survey_name %in% names(all_KR_surveys),
    ge_on_disk  = coalesce(ge_on_disk, FALSE)
  )


# Manual check of GE files that exist not aligned with rdhs that may be outdated

# KR          || GE
# ---------------------------
# KEKR42FL.ZIP = KEGE43FL.ZIP
# GHKR41FL.ZIP = GHGE42FL.ZIP



cat("KR surveys with a matching GE in metadata:  ",
    sum(!is.na(kr_ge_lookup$GE_file)), "\n")
cat("KR surveys with GE file on disk:            ",
    sum(kr_ge_lookup$ge_on_disk), "\n")
cat("KR surveys loaded from RDS:                 ",
    sum(kr_ge_lookup$kr_loaded), "\n")
cat("KR surveys loaded AND GE on disk:           ",
    sum(kr_ge_lookup$kr_loaded & kr_ge_lookup$ge_on_disk), "\n")
cat("KR surveys loaded but NO GE on disk:        ",
    sum(kr_ge_lookup$kr_loaded & !kr_ge_lookup$ge_on_disk), "\n")


# ================================================================
# 4. Initialise tracking log (mirrors KR load script schema)
# ================================================================

gps_log <- data.frame(
  survey_name       = character(),
  SurveyId          = character(),
  CountryName       = character(),
  SurveyYear        = integer(),
  DHS_CountryCode   = character(),
  KR_file           = character(),
  GE_file           = character(),
  status            = character(),   # OK / NO_GE_FILE / NO_SHP / FAILED_SHP / NO_DHSCLUST / NOT_LOADED
  n_children        = integer(),
  n_clusters_kr     = integer(),
  n_clusters_ge     = integer(),
  clusters_match    = logical(),
  kr_interview_year = character(),
  ge_year           = integer(),
  year_match        = character(),
  cols_added        = character(),
  error_message     = character(),
  stringsAsFactors  = FALSE
)

all_KR_surveys_with_gps <- list()


# ================================================================
# 5. Attach GPS to each loaded KR survey
# ================================================================

# Only process surveys that are in the loaded RDS
kr_to_process <- kr_ge_lookup %>% filter(kr_loaded)

message("\nAttaching GPS to ", nrow(kr_to_process), " KR surveys...\n")

for (i in seq_len(nrow(kr_to_process))) {
  
  row         <- kr_to_process[i, ]
  survey_name <- row$survey_name
  ge_zip      <- if (row$ge_on_disk) row$zip_path else NA_character_
  
  message(sprintf("[%d/%d] %s", i, nrow(kr_to_process), survey_name))
  
  # Retrieve the already-loaded KR data frame
  df <- all_KR_surveys[[survey_name]]
  
  # ── Extract interview year from KR data ──────────────────────
  if ("v007" %in% names(df)) {
    raw_years <- unique(df$v007)
    raw_years <- ifelse(raw_years < 100, raw_years + 1900, raw_years)
    kr_interview_year <- paste(sort(raw_years), collapse = ", ")
  } else {
    kr_interview_year <- NA_character_
  }
  
  n_children    <- nrow(df)
  n_clusters_kr <- if ("v001" %in% names(df)) n_distinct(df$v001) else NA_integer_
  
  # ── No GE file available ─────────────────────────────────────
  if (is.na(ge_zip)) {
    
    df$LATNUM   <- NA_real_
    df$LONGNUM  <- NA_real_
    df$gps_file <- NA_character_
    
    gps_log <- bind_rows(gps_log, data.frame(
      survey_name       = survey_name,
      SurveyId          = row$SurveyId,
      CountryName       = row$CountryName,
      SurveyYear        = as.integer(row$SurveyYear),
      DHS_CountryCode   = row$DHS_CountryCode,
      KR_file           = row$KR_file,
      GE_file           = coalesce(row$GE_file, NA_character_),
      status            = "NO_GE_FILE",
      n_children        = n_children,
      n_clusters_kr     = n_clusters_kr,
      n_clusters_ge     = NA_integer_,
      clusters_match    = NA,
      kr_interview_year = kr_interview_year,
      ge_year           = NA_integer_,
      year_match        = NA_character_,
      cols_added        = NA_character_,
      error_message     = NA_character_,
      stringsAsFactors  = FALSE
    ))
    
    all_KR_surveys_with_gps[[survey_name]] <- df
    message("  ⚠️  No GE file on disk")
    next
  }
  
  # ── Unzip and locate shapefile ────────────────────────────────
  temp_dir <- file.path(tempdir(), survey_name)
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  
  unzip_result <- tryCatch(
    { unzip(ge_zip, exdir = temp_dir); "ok" },
    error = function(e) e$message
  )
  
  shp_file <- list.files(temp_dir, pattern = "\\.shp$",
                         full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  
  if (unzip_result != "ok" || length(shp_file) == 0) {
    
    df$LATNUM   <- NA_real_
    df$LONGNUM  <- NA_real_
    df$gps_file <- basename(ge_zip)
    
    gps_log <- bind_rows(gps_log, data.frame(
      survey_name       = survey_name,
      SurveyId          = row$SurveyId,
      CountryName       = row$CountryName,
      SurveyYear        = as.integer(row$SurveyYear),
      DHS_CountryCode   = row$DHS_CountryCode,
      KR_file           = row$KR_file,
      GE_file           = row$GE_file,
      status            = "NO_SHP",
      n_children        = n_children,
      n_clusters_kr     = n_clusters_kr,
      n_clusters_ge     = NA_integer_,
      clusters_match    = NA,
      kr_interview_year = kr_interview_year,
      ge_year           = NA_integer_,
      year_match        = NA_character_,
      cols_added        = NA_character_,
      error_message     = ifelse(unzip_result != "ok",
                                 paste("Unzip failed:", unzip_result),
                                 "No .shp found after unzip"),
      stringsAsFactors  = FALSE
    ))
    
    all_KR_surveys_with_gps[[survey_name]] <- df
    message("  ⚠️  No shapefile found in GE zip")
    next
  }
  
  # ── Read shapefile ────────────────────────────────────────────
  ge_sf <- tryCatch(
    sf::st_read(shp_file[1], quiet = TRUE),
    error = function(e) {
      message("  ERROR reading shapefile: ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(ge_sf)) {
    
    df$LATNUM   <- NA_real_
    df$LONGNUM  <- NA_real_
    df$gps_file <- basename(ge_zip)
    
    gps_log <- bind_rows(gps_log, data.frame(
      survey_name       = survey_name,
      SurveyId          = row$SurveyId,
      CountryName       = row$CountryName,
      SurveyYear        = as.integer(row$SurveyYear),
      DHS_CountryCode   = row$DHS_CountryCode,
      KR_file           = row$KR_file,
      GE_file           = row$GE_file,
      status            = "FAILED_SHP",
      n_children        = n_children,
      n_clusters_kr     = n_clusters_kr,
      n_clusters_ge     = NA_integer_,
      clusters_match    = NA,
      kr_interview_year = kr_interview_year,
      ge_year           = NA_integer_,
      year_match        = NA_character_,
      cols_added        = NA_character_,
      error_message     = "st_read() failed",
      stringsAsFactors  = FALSE
    ))
    
    all_KR_surveys_with_gps[[survey_name]] <- df
    next
  }
  
  # ── Check for DHSCLUST join key ───────────────────────────────
  if (!"DHSCLUST" %in% names(ge_sf)) {
    
    df$LATNUM   <- NA_real_
    df$LONGNUM  <- NA_real_
    df$gps_file <- basename(ge_zip)
    
    gps_log <- bind_rows(gps_log, data.frame(
      survey_name       = survey_name,
      SurveyId          = row$SurveyId,
      CountryName       = row$CountryName,
      SurveyYear        = as.integer(row$SurveyYear),
      DHS_CountryCode   = row$DHS_CountryCode,
      KR_file           = row$KR_file,
      GE_file           = row$GE_file,
      status            = "NO_DHSCLUST",
      n_children        = n_children,
      n_clusters_kr     = n_clusters_kr,
      n_clusters_ge     = NA_integer_,
      clusters_match    = NA,
      kr_interview_year = kr_interview_year,
      ge_year           = NA_integer_,
      year_match        = NA_character_,
      cols_added        = NA_character_,
      error_message     = "DHSCLUST column absent from shapefile",
      stringsAsFactors  = FALSE
    ))
    
    all_KR_surveys_with_gps[[survey_name]] <- df
    message("  ⚠️  DHSCLUST missing from shapefile")
    next
  }
  
  # ── Join GE columns to KR ─────────────────────────────────────
  ge_df       <- sf::st_drop_geometry(ge_sf)
  ge_year     <- if ("DHSYEAR" %in% names(ge_df)) as.integer(ge_df$DHSYEAR[1]) else NA_integer_
  
  # Only add columns not already in KR (keep DHSCLUST as join key)
  cols_to_add <- union("DHSCLUST", setdiff(names(ge_df), names(df)))
  ge_df       <- ge_df %>% dplyr::select(dplyr::all_of(cols_to_add))
  
  df <- df %>%
    dplyr::left_join(ge_df, by = c("v001" = "DHSCLUST"))
  
  df$gps_file <- basename(ge_zip)
  
  # ── Diagnostics ───────────────────────────────────────────────
  n_clusters_ge  <- n_distinct(ge_df$DHSCLUST)
  clusters_match <- n_clusters_kr == n_clusters_ge
  
  if (!clusters_match) {
    message(sprintf("  WARNING: KR has %d clusters, GE has %d clusters",
                    n_clusters_kr, n_clusters_ge))
  }
  
  # Year match
  kr_years   <- suppressWarnings(as.integer(
    unlist(strsplit(kr_interview_year, ",\\s*"))
  ))
  year_match <- if (is.na(ge_year) || all(is.na(kr_years))) {
    "not_matched"
  } else if (ge_year %in% kr_years) {
    "match"
  } else {
    "not_matched"
  }
  
  cols_added_str <- paste(
    setdiff(names(ge_df), "DHSCLUST"), collapse = ", "
  )
  
  gps_log <- bind_rows(gps_log, data.frame(
    survey_name       = survey_name,
    SurveyId          = row$SurveyId,
    CountryName       = row$CountryName,
    SurveyYear        = as.integer(row$SurveyYear),
    DHS_CountryCode   = row$DHS_CountryCode,
    KR_file           = row$KR_file,
    GE_file           = row$GE_file,
    status            = "OK",
    n_children        = n_children,
    n_clusters_kr     = n_clusters_kr,
    n_clusters_ge     = n_clusters_ge,
    clusters_match    = clusters_match,
    kr_interview_year = kr_interview_year,
    ge_year           = ge_year,
    year_match        = year_match,
    cols_added        = cols_added_str,
    error_message     = NA_character_,
    stringsAsFactors  = FALSE
  ))
  
  all_KR_surveys_with_gps[[survey_name]] <- df
  message(sprintf("  ✅ GPS attached: %d clusters | year_match=%s | cols: %s",
                  n_clusters_ge, year_match, cols_added_str))
}


# ================================================================
# 6. Summary
# ================================================================

cat("\n=== GPS ATTACHMENT SUMMARY ===\n")
cat("Surveys processed:            ", length(all_KR_surveys_with_gps), "\n")
cat("GPS attached (OK):            ",
    sum(gps_log$status == "OK",          na.rm = TRUE), "\n")
cat("No GE file on disk:           ",
    sum(gps_log$status == "NO_GE_FILE",  na.rm = TRUE), "\n")
cat("No shapefile in GE zip:       ",
    sum(gps_log$status == "NO_SHP",      na.rm = TRUE), "\n")
cat("Shapefile read failed:        ",
    sum(gps_log$status == "FAILED_SHP",  na.rm = TRUE), "\n")
cat("DHSCLUST absent:              ",
    sum(gps_log$status == "NO_DHSCLUST", na.rm = TRUE), "\n")
cat("Cluster count mismatches:     ",
    sum(!gps_log$clusters_match, na.rm = TRUE), "\n")
cat("Year mismatches:              ",
    sum(gps_log$year_match == "not_matched", na.rm = TRUE), "\n")

cat("\nCountries with GPS attached:\n")
print(sort(unique(
  gps_log$CountryName[gps_log$status == "OK"]
)))


# ================================================================
# 7. Save
# ================================================================

saveRDS(all_KR_surveys_with_gps, file = output_rds)
message("\n✅ Saved: ", output_rds)

writexl::write_xlsx(
  list(
    gps_log     = gps_log,
    kr_ge_lookup = kr_ge_lookup
  ),
  output_log
)
message("✅ Log saved: ", output_log)

message("\nKR + GPS dataset ready for pipeline.")


library(dplyr)
library(purrr)
library(tools)
library(rdhs)
library(haven)
library(writexl)
library(stringr)


# CONVERT ALL NAMES TO UPPERCASE FOR EASY MATCHING 


# ================================================================
# CONFIGURATION
# ================================================================

zip_folder <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS (download from site, zip only)"

output_rds  <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/all_KR_surveys.rds"
output_log  <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/KR_load_log.xlsx"


zip_files <- list.files(zip_folder, pattern = "\\.zip$", ignore.case = TRUE, full.names = TRUE)

for (f in zip_files) {
  new_name <- file.path(dirname(f), toupper(basename(f)))
  
  if (f != new_name) {
    file.rename(f, new_name)
  }
}

zip_files <- list.files(zip_folder, pattern = "\\.zip$", ignore.case = TRUE, full.names = TRUE)

for (zip_path in zip_files) {

  message("Processing: ", basename(zip_path))
  
  # Temporary extraction folder
  temp_dir <- file.path(zip_folder, "temp_unzip")
  dir.create(temp_dir, showWarnings = FALSE)
  
  # Unzip
  unzip(zip_path, exdir = temp_dir)
  
  # List extracted files
  extracted_files <- list.files(temp_dir, recursive = TRUE, full.names = TRUE)
  
  # Rename all files to uppercase
  for (f in extracted_files) {
    new_f <- file.path(dirname(f), toupper(basename(f)))
    
    if (f != new_f) {
      file.rename(f, new_f)
    }
  }
  
  # Remove original zip
  file.remove(zip_path)
  
  # Recreate zip with uppercase contents
  new_zip <- file.path(zip_folder, toupper(basename(zip_path)))
  
  old_wd <- getwd()
  setwd(temp_dir)
  
  files_to_zip <- list.files(".", recursive = TRUE)
  zip::zip(zipfile = new_zip, files = files_to_zip)
  
  setwd(old_wd)
  
  # Clean temp folder
  unlink(temp_dir, recursive = TRUE)
}


# ================================================================
# 1. Get metadata for all DHS datasets from rdhs
#    This gives us SurveyId, CountryName, SurveyYear etc.
#    for every file — we filter to KR only
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

# KR metadata: one row per KR survey
kr_meta <- all_datasets %>%
  filter(file_type == "KR") %>%
  dplyr::rename(KR_file = FileName) %>%
  dplyr::select(SurveyId, KR_file, CountryName,
                SurveyYear, DHS_CountryCode) %>%
  mutate(survey_name = tools::file_path_sans_ext(KR_file))

cat("KR surveys in rdhs metadata:", nrow(kr_meta), "\n")

# ================================================================
# 2. Find KR zip files on disk
#    Match disk files to rdhs metadata via filename
# ================================================================

zip_files    <- list.files(zip_folder, pattern = "\\.zip$",
                           full.names = TRUE, ignore.case = TRUE)
KR_zip_files <- zip_files[grepl("KR", basename(zip_files),
                                ignore.case = TRUE)]

cat("KR zip files found on disk:", length(KR_zip_files), "\n")

# Build disk lookup: survey_name → full path
disk_lookup <- data.frame(
  zip_path    = KR_zip_files,
  survey_name = tools::file_path_sans_ext(basename(KR_zip_files)),
  stringsAsFactors = FALSE
)

# Match disk files to rdhs metadata
kr_lookup <- kr_meta %>%
  left_join(disk_lookup, by = "survey_name") %>%
  mutate(
    on_disk = !is.na(zip_path) & file.exists(coalesce(zip_path, ""))
  )

cat("KR surveys matched to disk file:", sum(kr_lookup$on_disk), "\n")
cat("KR surveys in metadata but not on disk:",
    sum(!kr_lookup$on_disk), "\n")

# ================================================================
# 3. Initialise tracking with consistent data types
# ================================================================

pipeline_log <- data.frame(
  survey_name     = character(),
  SurveyId        = character(),
  CountryName     = character(),
  SurveyYear      = integer(),      # ← Consistent integer type
  DHS_CountryCode = character(),
  KR_file         = character(),
  status          = character(),
  n_rows_raw      = integer(),
  n_cols_raw      = integer(),
  n_children      = integer(),
  n_women         = integer(),
  n_clusters      = integer(),
  has_caseid      = logical(),
  has_bidx        = logical(),
  has_v001        = logical(),
  has_v002        = logical(),
  has_v003        = logical(),
  has_b3          = logical(),
  has_h9          = logical(),
  error_message   = character(),
  stringsAsFactors = FALSE
)




all_KR_surveys <- list()

# ================================================================
# 4. Load each KR survey
# ================================================================

kr_to_load <- kr_lookup %>% filter(on_disk)

message("\nLoading ", nrow(kr_to_load), " KR surveys...\n")

for (i in seq_len(nrow(kr_to_load))) {
  
  row         <- kr_to_load[i, ]
  survey_name <- row$survey_name
  zip_path    <- row$zip_path
  
  message(sprintf("[%d/%d] %s", i, nrow(kr_to_load), survey_name))
  
  # ── Read flat file ───────────────────────────────────────────
  df <- tryCatch(
    rdhs:::read_dhs_flat(zip_path),
    error = function(e) {
      message("  ERROR reading: ", e$message)
      return(NULL)  # Return NULL instead of list
    }
  )
  
  # ── Handle read failure ───────────────────────────────────────
  if (is.null(df)) {
    pipeline_log <- bind_rows(pipeline_log, data.frame(
      survey_name     = survey_name,
      SurveyId        = row$SurveyId,
      CountryName     = row$CountryName,
      SurveyYear      = as.integer(row$SurveyYear),  # ← Consistent type
      DHS_CountryCode = row$DHS_CountryCode,
      KR_file         = row$KR_file,
      status          = "FAILED_READ",
      n_rows_raw      = NA_integer_,
      n_cols_raw      = NA_integer_,
      n_children      = NA_integer_,
      n_women         = NA_integer_,
      n_clusters      = NA_integer_,
      has_caseid      = FALSE,
      has_bidx        = FALSE,
      has_v001        = FALSE,
      has_v002        = FALSE,
      has_v003        = FALSE,
      has_b3          = FALSE,
      has_h9          = FALSE,
      error_message   = paste("File read failed:", survey_name),
      stringsAsFactors = FALSE
    ))
    next
  }
  
  # ── Strip haven labels for cleaner downstream use ─────────────
  df <- df %>% mutate(across(where(is.labelled), zap_labels))
  
  # ── Attach survey metadata from rdhs ─────────────────────────
  df$SurveyId        <- row$SurveyId
  df$CountryName     <- row$CountryName
  df$SurveyYear      <- as.integer(row$SurveyYear)
  df$DHS_CountryCode <- row$DHS_CountryCode
  
  # ── Verify essential identifiers ─────────────────────────────
  required_id_vars <- c("caseid", "bidx", "v001", "v002", "v003")
  missing_ids      <- required_id_vars[!required_id_vars %in% names(df)]
  
  if (length(missing_ids) > 0) {
    message("  WARNING: missing identifier(s): ",
            paste(missing_ids, collapse = ", "))
  }
  
  # ── Count unique children, women, clusters ────────────────────
  n_children <- if (all(c("caseid","bidx") %in% names(df)))
    n_distinct(paste(df$caseid, df$bidx)) else NA_integer_
  
  n_women    <- if ("caseid" %in% names(df))
    n_distinct(df$caseid) else NA_integer_
  
  n_clusters <- if ("v001" %in% names(df))
    n_distinct(df$v001) else NA_integer_
  
  # ── Check for duplicates (should be 0 in raw KR) ─────────────
  n_rows <- nrow(df)
  if (!is.na(n_children) && n_rows != n_children) {
    message("  WARNING: ", n_rows - n_children,
            " duplicate rows in raw KR — deduplicating")
    df <- df %>%
      group_by(caseid, bidx) %>%
      slice(1) %>%
      ungroup()
  }
  
  # ── Log ───────────────────────────────────────────────────────
  pipeline_log <- bind_rows(pipeline_log, data.frame(
    survey_name     = survey_name,
    SurveyId        = row$SurveyId,
    CountryName     = row$CountryName,
    SurveyYear      = as.integer(row$SurveyYear),  # ← Consistent type
    DHS_CountryCode = row$DHS_CountryCode,
    KR_file         = row$KR_file,
    status          = ifelse(length(missing_ids) == 0,
                             "OK", "MISSING_IDS"),
    n_rows_raw      = n_rows,
    n_cols_raw      = ncol(df),
    n_children      = nrow(df),
    n_women         = n_women,
    n_clusters      = n_clusters,
    has_caseid      = "caseid" %in% names(df),
    has_bidx        = "bidx"   %in% names(df),
    has_v001        = "v001"   %in% names(df),
    has_v002        = "v002"   %in% names(df),
    has_v003        = "v003"   %in% names(df),
    has_b3          = "b3"     %in% names(df),
    has_h9          = "h9"     %in% names(df),
    error_message   = ifelse(length(missing_ids) > 0,
                             paste("Missing:", paste(missing_ids,
                                                     collapse = ",")),
                             NA_character_),
    stringsAsFactors = FALSE
  ))
  
  all_KR_surveys[[survey_name]] <- df
  message("  ✅ Loaded: ", nrow(df), " children | ",
          n_women, " women | ", n_clusters, " clusters")
}





# ================================================================
# 5. Summary
# ================================================================

cat("\n=== LOAD SUMMARY ===\n")
cat("Surveys successfully loaded: ", length(all_KR_surveys), "\n")
cat("Surveys failed:              ",
    sum(pipeline_log$status == "FAILED_READ", na.rm = TRUE), "\n")
cat("Surveys with missing IDs:    ",
    sum(pipeline_log$status == "MISSING_IDS", na.rm = TRUE), "\n")
cat("Total children loaded:       ",
    format(sum(pipeline_log$n_children, na.rm = TRUE),
           big.mark = ","), "\n")
cat("Total surveys with h9 (MCV1):",
    sum(pipeline_log$has_h9, na.rm = TRUE), "\n")
cat("Total surveys without h9:    ",
    sum(!pipeline_log$has_h9, na.rm = TRUE), "\n")

cat("\nCountries represented:\n")
print(sort(unique(pipeline_log$CountryName[pipeline_log$status == "OK"])))

cat("\nDHS phases represented:\n")
print(table(substr(names(all_KR_surveys), 5, 5)))

# ================================================================
# 6. Save
# ================================================================

saveRDS(all_KR_surveys, file = output_rds)
message("\n✅ Saved: ", output_rds)

writexl::write_xlsx(
  list(
    load_log       = pipeline_log,
    kr_lookup      = kr_lookup,
    not_on_disk    = kr_lookup %>% filter(!on_disk)
  ),
  output_log
)
message("✅ Log saved: ", output_log)

message("\nBase KR dataset ready for pipeline.")
message("Next step: attach GPS coordinates (all_KR_surveys_with_gps)")


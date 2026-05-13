library(dplyr)
library(sf)
library(tools)
library(rdhs)
library(haven)
library(writexl)
library(readr)
library(foreign)

# ================================================================
# CONFIGURATION
# ================================================================

zip_folder <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS (download from site, zip only)"

input_rds  <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/all_KR_surveys_with_gps_HR_IR_maps.rds"
output_rds <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/all_KR_surveys_with_gps_HR_IR_maps_GC.rds"
output_log <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/KR_GC_log.xlsx"


# ================================================================
# 1. Load the already-processed KR surveys (with GPS, HR, IR, maps)
# ================================================================

message("Loading KR surveys from RDS...")
all_KR_surveys <- readRDS(input_rds)
cat("KR surveys loaded from RDS:", length(all_KR_surveys), "\n")


# ================================================================
# 2. Get metadata for KR and GC datasets from rdhs
#    Mirror the approach in the GPS script: use all_datasets
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
    FileName    = toupper(FileName),
    survey_name = toupper(survey_name),
    file_type   = toupper(file_type)
  )

kr_meta <- all_datasets %>%
  filter(file_type == "KR") %>%
  dplyr::rename(KR_file = FileName) %>%
  dplyr::select(SurveyId, KR_file, CountryName,
                SurveyYear, DHS_CountryCode, survey_name)

gc_meta <- all_datasets %>%
  filter(file_type == "GC") %>%
  dplyr::rename(GC_file = FileName) %>%
  dplyr::select(SurveyId, GC_file)

cat("KR surveys in rdhs metadata:", nrow(kr_meta), "\n")
cat("GC files in rdhs metadata:  ", nrow(gc_meta),  "\n")


# ================================================================
# 3. Find GC zip files on disk
#    Mirror the GPS script disk-matching logic exactly
# ================================================================

zip_files    <- list.files(zip_folder, pattern = "\\.zip$",
                           full.names = TRUE, ignore.case = TRUE)
GC_zip_files <- zip_files[grepl("GC", basename(zip_files),
                                ignore.case = TRUE)]

cat("GC zip files found on disk:", length(GC_zip_files), "\n")

# Build disk lookup: survey_name -> full path
disk_lookup_gc <- data.frame(
  zip_path    = GC_zip_files,
  survey_name = toupper(tools::file_path_sans_ext(basename(GC_zip_files))),
  stringsAsFactors = FALSE
)

# Join GC metadata to its disk path via survey_name
gc_on_disk <- gc_meta %>%
  mutate(survey_name = toupper(tools::file_path_sans_ext(GC_file))) %>%
  left_join(disk_lookup_gc, by = "survey_name") %>%
  mutate(
    gc_on_disk = !is.na(zip_path) & file.exists(coalesce(zip_path, ""))
  )

# Join KR metadata to GC info via SurveyId — one row per KR survey
kr_gc_lookup <- kr_meta %>%
  left_join(
    gc_on_disk %>% dplyr::select(SurveyId, GC_file, zip_path, gc_on_disk),
    by = "SurveyId"
  ) %>%
  mutate(
    # KR must also be in the loaded RDS to be processed
    kr_loaded   = survey_name %in% names(all_KR_surveys),
    gc_on_disk  = coalesce(gc_on_disk, FALSE)
  )

cat("KR surveys with a matching GC in metadata:  ",
    sum(!is.na(kr_gc_lookup$GC_file)), "\n")
cat("KR surveys with GC file on disk:            ",
    sum(kr_gc_lookup$gc_on_disk), "\n")
cat("KR surveys loaded from RDS:                 ",
    sum(kr_gc_lookup$kr_loaded), "\n")
cat("KR surveys loaded AND GC on disk:           ",
    sum(kr_gc_lookup$kr_loaded & kr_gc_lookup$gc_on_disk), "\n")
cat("KR surveys loaded but NO GC on disk:        ",
    sum(kr_gc_lookup$kr_loaded & !kr_gc_lookup$gc_on_disk), "\n")


# ================================================================
# 4. Helper: read a flat file from inside a zip
#    GC datasets are distributed as .dta (Stata), .sav (SPSS),
#    or .dat/.csv flat files. This function tries each format.
# ================================================================

read_gc_flat_file <- function(temp_dir) {
  
  # Priority order: .dta > .sav > .dbf > .csv > .dat
  dta_files <- list.files(temp_dir, pattern = "\\.dta$",
                          full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (length(dta_files) > 0) {
    message("    Reading .dta file: ", basename(dta_files[1]))
    gc_df <- haven::read_dta(dta_files[1])
    return(gc_df)
  }
  
  sav_files <- list.files(temp_dir, pattern = "\\.sav$",
                          full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (length(sav_files) > 0) {
    message("    Reading .sav file: ", basename(sav_files[1]))
    gc_df <- haven::read_sav(sav_files[1])
    return(gc_df)
  }
  
  dbf_files <- list.files(temp_dir, pattern = "\\.dbf$",
                          full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (length(dbf_files) > 0) {
    message("    Reading .dbf file: ", basename(dbf_files[1]))
    gc_df <- foreign::read.dbf(dbf_files[1], as.is = TRUE)
    return(gc_df)
  }
  
  csv_files <- list.files(temp_dir, pattern = "\\.csv$",
                          full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (length(csv_files) > 0) {
    message("    Reading .csv file: ", basename(csv_files[1]))
    gc_df <- readr::read_csv(csv_files[1], show_col_types = FALSE)
    return(gc_df)
  }
  
  dat_files <- list.files(temp_dir, pattern = "\\.dat$",
                          full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (length(dat_files) > 0) {
    message("    Reading .dat file: ", basename(dat_files[1]))
    gc_df <- readr::read_csv(dat_files[1], show_col_types = FALSE)
    return(gc_df)
  }
  
  return(NULL)
}


# ================================================================
# 5. Initialise tracking log
# ================================================================

gc_log <- data.frame(
  survey_name       = character(),
  SurveyId          = character(),
  CountryName       = character(),
  SurveyYear        = integer(),
  DHS_CountryCode   = character(),
  KR_file           = character(),
  GC_file           = character(),
  status            = character(),   # OK / NO_GC_FILE / NO_FLAT_FILE / FAILED_READ / NO_DHSCLUST / NOT_LOADED
  n_children        = integer(),
  n_clusters_kr     = integer(),
  n_clusters_gc     = integer(),
  clusters_match    = logical(),
  n_gc_vars_added   = integer(),
  cols_added        = character(),
  error_message     = character(),
  stringsAsFactors  = FALSE
)

all_KR_surveys_with_GC <- list()


# ================================================================
# 6. Attach GC covariates to each loaded KR survey
# ================================================================

kr_to_process <- kr_gc_lookup %>% filter(kr_loaded)

message("\nAttaching GC covariates to ", nrow(kr_to_process), " KR surveys...\n")

for (i in seq_len(nrow(kr_to_process))) {
  
  row         <- kr_to_process[i, ]
  survey_name <- row$survey_name
  gc_zip      <- if (row$gc_on_disk) row$zip_path else NA_character_
  
  message(sprintf("[%d/%d] %s", i, nrow(kr_to_process), survey_name))
  
  # Retrieve the already-loaded KR data frame
  df <- all_KR_surveys[[survey_name]]
  
  n_children    <- nrow(df)
  n_clusters_kr <- if ("v001" %in% names(df)) n_distinct(df$v001) else NA_integer_
  
  # ── No GC file available ─────────────────────────────────────
  if (is.na(gc_zip)) {
    
    gc_log <- bind_rows(gc_log, data.frame(
      survey_name       = survey_name,
      SurveyId          = row$SurveyId,
      CountryName       = row$CountryName,
      SurveyYear        = as.integer(row$SurveyYear),
      DHS_CountryCode   = row$DHS_CountryCode,
      KR_file           = row$KR_file,
      GC_file           = coalesce(row$GC_file, NA_character_),
      status            = "NO_GC_FILE",
      n_children        = n_children,
      n_clusters_kr     = n_clusters_kr,
      n_clusters_gc     = NA_integer_,
      clusters_match    = NA,
      n_gc_vars_added   = 0L,
      cols_added        = NA_character_,
      error_message     = NA_character_,
      stringsAsFactors  = FALSE
    ))
    
    all_KR_surveys_with_GC[[survey_name]] <- df
    message("  ⚠️  No GC file on disk")
    next
  }
  
  # ── Unzip GC file ────────────────────────────────────────────
  temp_dir <- file.path(tempdir(), paste0("GC_", survey_name))
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  
  unzip_result <- tryCatch(
    { unzip(gc_zip, exdir = temp_dir); "ok" },
    error = function(e) e$message
  )
  
  if (unzip_result != "ok") {
    
    gc_log <- bind_rows(gc_log, data.frame(
      survey_name       = survey_name,
      SurveyId          = row$SurveyId,
      CountryName       = row$CountryName,
      SurveyYear        = as.integer(row$SurveyYear),
      DHS_CountryCode   = row$DHS_CountryCode,
      KR_file           = row$KR_file,
      GC_file           = row$GC_file,
      status            = "UNZIP_FAILED",
      n_children        = n_children,
      n_clusters_kr     = n_clusters_kr,
      n_clusters_gc     = NA_integer_,
      clusters_match    = NA,
      n_gc_vars_added   = 0L,
      cols_added        = NA_character_,
      error_message     = paste("Unzip failed:", unzip_result),
      stringsAsFactors  = FALSE
    ))
    
    all_KR_surveys_with_GC[[survey_name]] <- df
    message("  ⚠️  Unzip failed")
    next
  }
  
  # ── Read the flat file inside the zip ─────────────────────────
  gc_df <- tryCatch(
    read_gc_flat_file(temp_dir),
    error = function(e) {
      message("  ERROR reading GC file: ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(gc_df)) {
    
    gc_log <- bind_rows(gc_log, data.frame(
      survey_name       = survey_name,
      SurveyId          = row$SurveyId,
      CountryName       = row$CountryName,
      SurveyYear        = as.integer(row$SurveyYear),
      DHS_CountryCode   = row$DHS_CountryCode,
      KR_file           = row$KR_file,
      GC_file           = row$GC_file,
      status            = "NO_FLAT_FILE",
      n_children        = n_children,
      n_clusters_kr     = n_clusters_kr,
      n_clusters_gc     = NA_integer_,
      clusters_match    = NA,
      n_gc_vars_added   = 0L,
      cols_added        = NA_character_,
      error_message     = "No readable flat file (.dta/.sav/.dbf/.csv/.dat) found in zip",
      stringsAsFactors  = FALSE
    ))
    
    all_KR_surveys_with_GC[[survey_name]] <- df
    message("  ⚠️  No readable flat file found in GC zip")
    next
  }
  
  # ── Zap labels from GC data for clean joins ──────────────────
  gc_df <- gc_df %>%
    mutate(across(where(haven::is.labelled), haven::zap_labels))
  
  # ── Standardise column names to uppercase (DHS convention) ───
  names(gc_df) <- toupper(names(gc_df))
  
  # ── Check for DHSCLUST join key ──────────────────────────────
  if (!"DHSCLUST" %in% names(gc_df)) {
    
    gc_log <- bind_rows(gc_log, data.frame(
      survey_name       = survey_name,
      SurveyId          = row$SurveyId,
      CountryName       = row$CountryName,
      SurveyYear        = as.integer(row$SurveyYear),
      DHS_CountryCode   = row$DHS_CountryCode,
      KR_file           = row$KR_file,
      GC_file           = row$GC_file,
      status            = "NO_DHSCLUST",
      n_children        = n_children,
      n_clusters_kr     = n_clusters_kr,
      n_clusters_gc     = NA_integer_,
      clusters_match    = NA,
      n_gc_vars_added   = 0L,
      cols_added        = NA_character_,
      error_message     = paste("DHSCLUST column absent. Columns found:",
                                paste(names(gc_df), collapse = ", ")),
      stringsAsFactors  = FALSE
    ))
    
    all_KR_surveys_with_GC[[survey_name]] <- df
    message("  ⚠️  DHSCLUST missing from GC file")
    next
  }
  
  # ── Join GC columns to KR ─────────────────────────────────────
  # Ensure v001 in KR matches DHSCLUST type in GC
  gc_df$DHSCLUST <- as.numeric(gc_df$DHSCLUST)
  
  # Only add columns not already in KR (keep DHSCLUST as join key)
  # Exclude DHSID, DHSCC, DHSYEAR, DHSREGCO, DHSREGNA if already present
  # (these are identifiers, not covariates — but include them if new)
  cols_to_add <- union("DHSCLUST", setdiff(names(gc_df), names(df)))
  gc_join     <- gc_df %>% dplyr::select(dplyr::all_of(cols_to_add))
  
  # Deduplicate GC rows per cluster (should be 1:1 but just in case)
  gc_join <- gc_join %>% distinct(DHSCLUST, .keep_all = TRUE)
  
  # Ensure v001 is numeric for a clean join
  df_joined <- df %>%
    mutate(v001 = as.numeric(v001)) %>%
    dplyr::left_join(gc_join, by = c("v001" = "DHSCLUST"))
  
  # ── Diagnostics ───────────────────────────────────────────────
  n_clusters_gc  <- n_distinct(gc_join$DHSCLUST)
  clusters_match <- n_clusters_kr == n_clusters_gc
  
  if (!clusters_match) {
    message(sprintf("  WARNING: KR has %d clusters, GC has %d clusters",
                    n_clusters_kr, n_clusters_gc))
  }
  
  new_cols        <- setdiff(names(gc_join), "DHSCLUST")
  cols_added_str  <- paste(new_cols, collapse = ", ")
  n_gc_vars_added <- length(new_cols)
  
  gc_log <- bind_rows(gc_log, data.frame(
    survey_name       = survey_name,
    SurveyId          = row$SurveyId,
    CountryName       = row$CountryName,
    SurveyYear        = as.integer(row$SurveyYear),
    DHS_CountryCode   = row$DHS_CountryCode,
    KR_file           = row$KR_file,
    GC_file           = row$GC_file,
    status            = "OK",
    n_children        = n_children,
    n_clusters_kr     = n_clusters_kr,
    n_clusters_gc     = n_clusters_gc,
    clusters_match    = clusters_match,
    n_gc_vars_added   = n_gc_vars_added,
    cols_added        = cols_added_str,
    error_message     = NA_character_,
    stringsAsFactors  = FALSE
  ))
  
  all_KR_surveys_with_GC[[survey_name]] <- df_joined
  message(sprintf("  ✅ GC attached: %d clusters | %d new vars | cols: %s",
                  n_clusters_gc, n_gc_vars_added, cols_added_str))
  
  # Clean up temp directory
  unlink(temp_dir, recursive = TRUE)
}


# ================================================================
# 7. Handle any loaded surveys not in kr_to_process
#    (edge case: loaded but not in rdhs metadata)
# ================================================================

missing_surveys <- setdiff(names(all_KR_surveys), names(all_KR_surveys_with_GC))
if (length(missing_surveys) > 0) {
  message("\nCarrying forward ", length(missing_surveys),
          " surveys not found in rdhs metadata (no GC attempted)...")
  for (sname in missing_surveys) {
    all_KR_surveys_with_GC[[sname]] <- all_KR_surveys[[sname]]
  }
}


# ================================================================
# 8. Summary
# ================================================================

cat("\n=== GC COVARIATE ATTACHMENT SUMMARY ===\n")
cat("Surveys processed:            ", length(all_KR_surveys_with_GC), "\n")
cat("GC attached (OK):             ",
    sum(gc_log$status == "OK",           na.rm = TRUE), "\n")
cat("No GC file on disk:           ",
    sum(gc_log$status == "NO_GC_FILE",   na.rm = TRUE), "\n")
cat("Unzip failed:                 ",
    sum(gc_log$status == "UNZIP_FAILED", na.rm = TRUE), "\n")
cat("No flat file in GC zip:       ",
    sum(gc_log$status == "NO_FLAT_FILE", na.rm = TRUE), "\n")
cat("DHSCLUST absent:              ",
    sum(gc_log$status == "NO_DHSCLUST",  na.rm = TRUE), "\n")
cat("Cluster count mismatches:     ",
    sum(!gc_log$clusters_match, na.rm = TRUE), "\n")

cat("\nCountries with GC attached:\n")
print(sort(unique(
  gc_log$CountryName[gc_log$status == "OK"]
)))

# ── Show typical GC variables added ────────────────────────────
ok_logs <- gc_log %>% filter(status == "OK")
if (nrow(ok_logs) > 0) {
  all_gc_cols <- ok_logs$cols_added %>%
    strsplit(", ") %>%
    unlist() %>%
    table() %>%
    sort(decreasing = TRUE)
  cat("\nMost common GC variables added (across OK surveys):\n")
  print(head(all_gc_cols, 30))
}


# ================================================================
# 9. Save
# ================================================================

saveRDS(all_KR_surveys_with_GC, file = output_rds)
message("\n✅ Saved: ", output_rds)

writexl::write_xlsx(
  list(
    gc_log       = gc_log,
    kr_gc_lookup = kr_gc_lookup
  ),
  output_log
)
message("✅ Log saved: ", output_log)

message("\nKR + GPS + HR + IR + Maps + GC dataset ready for pipeline.")


# ================================================================
# 10. QC: Check GC variable coverage across surveys
# ================================================================

message("\n=== QC: GC VARIABLE COVERAGE ===\n")

gc_coverage <- purrr::map_dfr(names(all_KR_surveys_with_GC), function(sname) {
  
  df <- all_KR_surveys_with_GC[[sname]]
  
  # Typical GC variable prefixes (uppercase after our standardisation)
  # Common GC vars: UN_*, GN_*, EP_*, CH_*, CL_*, etc.
  gc_var_candidates <- grep("^(UN_|GN_|EP_|CH_|CL_|FP_|PO_|WS_|AN_|AH_|ML_|CN_)",
                            names(df), value = TRUE, ignore.case = TRUE)
  
  if (length(gc_var_candidates) == 0) {
    return(tibble::tibble(
      survey       = sname,
      n_gc_vars    = 0L,
      n_non_na     = 0L,
      pct_with_gc  = 0
    ))
  }
  
  n_non_na <- sum(!is.na(df[[gc_var_candidates[1]]]))
  
  tibble::tibble(
    survey       = sname,
    n_gc_vars    = length(gc_var_candidates),
    n_non_na     = n_non_na,
    pct_with_gc  = round(100 * n_non_na / nrow(df), 1)
  )
})

cat("Surveys with GC vars:  ", sum(gc_coverage$n_gc_vars > 0), "/",
    nrow(gc_coverage), "\n")
cat("Surveys with 0 GC vars:", sum(gc_coverage$n_gc_vars == 0), "\n")

# Show surveys with low coverage
low_coverage <- gc_coverage %>% filter(n_gc_vars > 0, pct_with_gc < 80)
if (nrow(low_coverage) > 0) {
  cat("\nSurveys with GC vars but <80% coverage:\n")
  print(low_coverage)
}

cat("\n✅ GC covariate attachment complete.\n")


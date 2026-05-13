library(dplyr)
library(stringr)
library(tools)
library(rdhs)
library(readxl)
library(purrr)
library(haven)

# ================================================================
# 1. Load variable list
# ================================================================
# 
vars <- read_xlsx("C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/variables_to_keep_classified.xlsx")

variables_to_keep <- vars %>%
  mutate(flag = coalesce(keep)) %>%
  filter(flag == 1) %>%
  distinct(variable) %>%
  pull(variable)

# ================================================================
# 2. Load base KR + GPS data (one row per child)
# ================================================================

all_KR_surveys_with_gps <- readRDS(
  "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/all_KR_surveys_with_gps.rds"
)

message("Surveys loaded: ", length(all_KR_surveys_with_gps))

# ================================================================
# 2. Pull HR file metadata from rdhs and match by SurveyId
# ================================================================

all_datasets <- dhs_datasets(fileFormat = "FL") %>%
  dplyr::select(SurveyId, CountryName, SurveyYear,
                DHS_CountryCode, FileName) %>%
  mutate(
    file_type = substr(FileName, 3, 4),  # e.g. HR, KR, IR, BR
    survey_name = tools::file_path_sans_ext(FileName)
  ) %>% 
  mutate(
    FileName = toupper(FileName),
    survey_name = toupper(survey_name),
    file_type = toupper(file_type)
  )

# Keep only HR files
hr_meta <- all_datasets %>%
  filter(file_type == "HR") %>%
  dplyr::rename(HR_file = FileName) %>%
  dplyr::select(SurveyId, CountryName, SurveyYear,
                DHS_CountryCode, HR_file)

# Keep only KR files to get SurveyId for each KR survey name
kr_meta <- all_datasets %>%
  filter(file_type == "KR") %>%
  dplyr::rename(KR_file = FileName) %>%
  dplyr::select(SurveyId, KR_file) %>%
  mutate(survey_name = tools::file_path_sans_ext(KR_file))

# Build lookup: KR survey name → HR file path via shared SurveyId
hr_lookup <- kr_meta %>%
  # Match KR survey names to those in our loaded data
  filter(survey_name %in% names(all_KR_surveys_with_gps)) %>%
  # Join to HR metadata on SurveyId
  left_join(hr_meta, by = "SurveyId") %>%
  mutate(
    local_hr_path = ifelse(
      !is.na(HR_file),
      file.path(
        "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS (download from site, zip only)",
        HR_file
      ),
      NA_character_
    ),
    hr_exists = !is.na(local_hr_path) &
      file.exists(coalesce(local_hr_path, ""))
  )

cat("KR surveys in data:          ", length(names(all_KR_surveys_with_gps)), "\n")
cat("KR surveys matched to HR:    ", sum(!is.na(hr_lookup$HR_file)), "\n")
cat("HR files found on disk:      ", sum(hr_lookup$hr_exists), "\n")
cat("KR surveys with no HR match: ", sum(is.na(hr_lookup$HR_file)), "\n")


# Check for any KR surveys not in hr_lookup at all
missing_from_lookup <- setdiff(
  names(all_KR_surveys_with_gps),
  hr_lookup$survey_name
)
if (length(missing_from_lookup) > 0) {
  cat("KR surveys not found in rdhs metadata:",
      length(missing_from_lookup), "\n")
  print(missing_from_lookup)
}

# ================================================================
# 4. Define join keys
#
# KR uses: v001 (cluster), v002 (household)
# HR uses: hv001 (cluster), hv002 (household)
#
# HR is ONE ROW PER HOUSEHOLD — joining on cluster + household
# gives a MANY-TO-ONE join (many children in same household
# all get the same household-level variables added as columns)
# Row count must NOT increase after join
# ================================================================

kr_id_vars <- c("v001", "v002", "v003", "caseid",
                "SurveyYear", "DHS_CountryCode", "CountryName",
                "LATNUM", "LONGNUM", "gps_file")

hr_join_keys_kr <- c("v001", "v002")   # names in KR
hr_join_keys_hr <- c("hv001", "hv002") # names in HR

# ================================================================
# 5. Main loop
# ================================================================

hr_status   <- list()
output_list <- list()

for (survey_name in names(all_KR_surveys_with_gps)) {
  
  message("Processing: ", survey_name)
  
  kr_data <- all_KR_surveys_with_gps[[survey_name]]
  if (is.null(kr_data) || nrow(kr_data) == 0) next
  
  # ── Step 1: Verify KR is one row per child ────────────────────
  n_before     <- nrow(kr_data)
  n_unique_ch  <- n_distinct(paste(kr_data$caseid, kr_data$bidx))
  
  if (n_before != n_unique_ch) {
    message("  WARNING: KR has ", n_before - n_unique_ch,
            " duplicate rows — deduplicating to one row per child")
    kr_data <- kr_data %>%
      group_by(caseid, bidx) %>%
      slice(1) %>%
      ungroup()
  }
  
  # ── Step 2: Keep only required KR variables ───────────────────
  kr_keep <- intersect(
    c(kr_id_vars, variables_to_keep),
    names(kr_data)
  )

  kr_data <- kr_data %>% dplyr::select(all_of(kr_keep))
  
  # ── Step 3: Get KR interview year ─────────────────────────────
  kr_year_str <- if ("v007" %in% names(kr_data)) {
    yr <- unique(as.numeric(zap_labels(kr_data$v007)))
    yr <- ifelse(yr < 100, yr + 1900, yr)
    paste(sort(yr), collapse = ", ")
  } else NA_character_
  
  # ── Step 4: Find HR file ──────────────────────────────────────
  hr_row <- hr_lookup %>% filter(survey_name == !!survey_name)
  hr_zip <- if (nrow(hr_row) == 1 && isTRUE(hr_row$hr_exists)) {
    hr_row$local_hr_path
  } else NA_character_
  
  if (is.na(hr_zip)) {
    message("  No HR file found — keeping KR only")
    
    # Attach survey metadata from hr_lookup if available
    if (nrow(hr_row) == 1) {
      kr_data$SurveyYear      <- hr_row$SurveyYear
      kr_data$DHS_CountryCode <- hr_row$DHS_CountryCode
      kr_data$CountryName     <- hr_row$CountryName
    }
    
    hr_status[[survey_name]] <- data.frame(
      survey_name     = survey_name,
      kr_year         = kr_year_str,
      hr_year         = NA_integer_,
      has_hr          = FALSE,
      hr_file         = NA_character_,
      n_children      = nrow(kr_data),
      n_hr_rows       = NA_integer_,
      hr_vars_added   = 0L,
      rows_after_join = nrow(kr_data),
      join_clean      = TRUE
    )
    output_list[[survey_name]] <- kr_data
    next
  }
  
  # ── Step 5: Read HR ───────────────────────────────────────────
  hr_data <- tryCatch(
    rdhs:::read_dhs_flat(hr_zip),
    error = function(e) {
      message("  Failed to read HR: ", e$message)
      NULL
    }
  )
  
  if (is.null(hr_data)) {
    hr_status[[survey_name]] <- data.frame(
      survey_name     = survey_name,
      kr_year         = kr_year_str,
      hr_year         = NA_integer_,
      has_hr          = FALSE,
      hr_file         = basename(hr_zip),
      n_children      = nrow(kr_data),
      n_hr_rows       = NA_integer_,
      hr_vars_added   = 0L,
      rows_after_join = nrow(kr_data),
      join_clean      = TRUE
    )
    output_list[[survey_name]] <- kr_data
    next
  }
  
  # ── Step 6: Extract HR year ───────────────────────────────────
  hr_year <- if ("hv007" %in% names(hr_data)) {
    yr <- as.numeric(zap_labels(hr_data$hv007))[1]
    as.integer(ifelse(yr < 100, yr + 1900, yr))
  } else NA_integer_
  
  # ── Step 7: Deduplicate HR to ONE ROW PER HOUSEHOLD ──────────
  # HR should be one row per household but enforce it
  # Join key: hv001 + hv002 uniquely identifies a household
  n_hr_before <- nrow(hr_data)
  hr_data <- hr_data %>%
    mutate(across(all_of(hr_join_keys_hr), as.numeric)) %>%
    group_by(hv001, hv002) %>%
    slice(1) %>%
    ungroup()
  
  if (nrow(hr_data) != n_hr_before) {
    message("  HR had ", n_hr_before - nrow(hr_data),
            " duplicate household rows — deduplicated")
  }
  
  # ── Step 8: Select only variables to add from HR ─────────────
  # Keep join keys + variables in variables_to_keep
  # Drop columns already present in KR to avoid conflicts
  kr_existing <- setdiff(names(kr_data), hr_join_keys_kr)
  
  hr_vars_to_add <- hr_data %>%
    dplyr::select(
      all_of(hr_join_keys_hr),
      any_of(variables_to_keep)
    ) %>%
    # Remove columns whose KR-equivalent already exists
    dplyr::select(-any_of(kr_existing)) %>%
    names() %>%
    setdiff(hr_join_keys_hr)
  
  if (length(hr_vars_to_add) == 0) {
    message("  No new HR variables to add")
    kr_data$SurveyYear      <- hr_row$SurveyYear
    kr_data$DHS_CountryCode <- hr_row$DHS_CountryCode
    kr_data$CountryName     <- hr_row$CountryName
    
    hr_status[[survey_name]] <- data.frame(
      survey_name     = survey_name,
      kr_year         = kr_year_str,
      hr_year         = hr_year,
      has_hr          = TRUE,
      hr_file         = basename(hr_zip),
      n_children      = nrow(kr_data),
      n_hr_rows       = nrow(hr_data),
      hr_vars_added   = 0L,
      rows_after_join = nrow(kr_data),
      join_clean      = TRUE
    )
    output_list[[survey_name]] <- kr_data
    next
  }
  
  hr_slim <- hr_data %>%
    dplyr::select(all_of(c(hr_join_keys_hr, hr_vars_to_add)))
  
  # ── Step 9: Join HR to KR ─────────────────────────────────────
  # MANY-TO-ONE: many children in same household → one household row
  # Row count must stay the same as KR
  n_kr_before <- nrow(kr_data)
  
  kr_data <- kr_data %>%
    mutate(across(all_of(hr_join_keys_kr), as.numeric)) %>%
    left_join(
      hr_slim,
      by = setNames(hr_join_keys_hr, hr_join_keys_kr),
      relationship = "many-to-one"  # enforces no fan-out
    )
  
  n_kr_after <- nrow(kr_data)
  join_clean  <- n_kr_after == n_kr_before
  
  if (!join_clean) {
    message("  ❌ JOIN CREATED DUPLICATES: ",
            n_kr_before, " → ", n_kr_after, " rows — forcing dedup")
    kr_data <- kr_data %>%
      group_by(caseid, bidx) %>%
      slice(1) %>%
      ungroup()
  } else {
    message("  ✅ Join clean: ", n_kr_after, " rows (",
            length(hr_vars_to_add), " HR variables added)")
  }
  
  # ── Step 10: Attach survey metadata ──────────────────────────
  kr_data$SurveyYear      <- hr_row$SurveyYear
  kr_data$DHS_CountryCode <- hr_row$DHS_CountryCode
  kr_data$CountryName     <- hr_row$CountryName
  
  # ── Record status ─────────────────────────────────────────────
  hr_status[[survey_name]] <- data.frame(
    survey_name     = survey_name,
    kr_year         = kr_year_str,
    hr_year         = hr_year,
    has_hr          = TRUE,
    hr_file         = basename(hr_zip),
    n_children      = nrow(kr_data),
    n_hr_rows       = nrow(hr_data),
    hr_vars_added   = length(hr_vars_to_add),
    rows_after_join = n_kr_after,
    join_clean      = join_clean
  )
  
  output_list[[survey_name]] <- kr_data
  
  rm(hr_data, hr_slim)
  gc()
}

# ================================================================
# 6. Compile status and verify
# ================================================================

hr_status_df <- dplyr::bind_rows(hr_status)

# Year match
hr_status_df$year_match <- sapply(seq_len(nrow(hr_status_df)), function(i) {
  kr_yrs <- as.numeric(unlist(strsplit(hr_status_df$kr_year[i], ",")))
  hr_yr  <- hr_status_df$hr_year[i]
  if (is.na(hr_yr) || all(is.na(kr_yrs))) "not matched"
  else if (hr_yr %in% kr_yrs) "match"
  else "not matched"
})

cat("\n=== HR MERGE STATUS SUMMARY ===\n")
cat("Surveys processed:           ", nrow(hr_status_df), "\n")
cat("Surveys with HR attached:    ", sum(hr_status_df$has_hr, na.rm = TRUE), "\n")
cat("Surveys without HR:          ", sum(!hr_status_df$has_hr, na.rm = TRUE), "\n")
cat("Clean joins:                 ", sum(hr_status_df$join_clean, na.rm = TRUE), "\n")
cat("Joins that created duplicates:", sum(!hr_status_df$join_clean, na.rm = TRUE), "\n")
cat("Year matches:                ", sum(hr_status_df$year_match == "match", na.rm = TRUE), "\n")
cat("Total children:              ",
    format(sum(hr_status_df$n_children, na.rm = TRUE), big.mark = ","), "\n")

# Final verification
cat("\n=== FINAL VERIFICATION ===\n")
final_check <- purrr::map_dfr(names(output_list), function(sn) {
  df <- output_list[[sn]]
  data.frame(
    survey   = sn,
    n_rows   = nrow(df),
    n_unique = n_distinct(paste(df$caseid, df$bidx)),
    clean    = nrow(df) == n_distinct(paste(df$caseid, df$bidx))
  )
})

cat("Surveys at one row per child:", sum(final_check$clean), "/",
    nrow(final_check), "\n")

dirty <- final_check %>% filter(!clean)
if (nrow(dirty) > 0) {
  cat("❌ Surveys still with duplicates:\n")
  print(dirty)
} else {
  cat("✅ All surveys confirmed at one row per child\n")
}

# ================================================================
# 7. Save
# ================================================================

saveRDS(
  output_list,
  "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/all_KR_surveys_with_gps_HR.rds"
)

writexl::write_xlsx(
  list(
    hr_merge_status = hr_status_df,
    final_check     = final_check
  ),
  "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/hr_merge_status.xlsx"
)

message("\n✅ Done. Saved to all_KR_surveys_with_gps_HR_rds")
message("Total surveys: ", length(output_list))
message("Total children: ", format(sum(final_check$n_rows), big.mark = ","))


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
 
vars <- read_xlsx("C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/variables_to_keep_classified.xlsx")
 
variables_to_keep <- vars %>%
  mutate(flag = coalesce(keep)) %>%
  filter(flag == 1) %>%
  distinct(variable) %>%
  pull(variable)
 
# ================================================================
# 2. Load base KR+GPS+HR data (one row per child — already clean)
# ================================================================
 
all_KR_surveys_with_gps_HR <- readRDS(
  "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/all_KR_surveys_with_gps_HR.rds"
)
 
message("Surveys loaded: ", length(all_KR_surveys_with_gps_HR))
 
# ================================================================
# 3. Pull IR file metadata from rdhs — match via SurveyId
#    (same approach as the HR merge script)
# ================================================================
 
all_datasets <- dhs_datasets(fileFormat = "FL") %>%
  dplyr::select(SurveyId, CountryName, SurveyYear,
                DHS_CountryCode, FileName) %>%
  mutate(
    file_type   = substr(FileName, 3, 4),
    survey_name = tools::file_path_sans_ext(FileName)
  ) %>%
  mutate(
    FileName    = toupper(FileName),
    survey_name = toupper(survey_name),
    file_type   = toupper(file_type)
  )
 
# Keep only IR files
ir_meta <- all_datasets %>%
  filter(file_type == "IR") %>%
  dplyr::rename(IR_file = FileName) %>%
  dplyr::select(SurveyId, CountryName, SurveyYear,
                DHS_CountryCode, IR_file)
 
# Keep only KR files to get SurveyId for each KR survey name
kr_meta <- all_datasets %>%
  filter(file_type == "KR") %>%
  dplyr::rename(KR_file = FileName) %>%
  dplyr::select(SurveyId, KR_file) %>%
  mutate(survey_name = tools::file_path_sans_ext(KR_file))
 
# Build lookup: KR survey name → IR file path via shared SurveyId
ir_lookup <- kr_meta %>%
  filter(survey_name %in% names(all_KR_surveys_with_gps_HR)) %>%
  left_join(ir_meta, by = "SurveyId") %>%
  mutate(
    local_ir_path = ifelse(
      !is.na(IR_file),
      file.path(
        "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS (download from site, zip only)",
        IR_file
      ),
      NA_character_
    ),
    ir_exists = !is.na(local_ir_path) &
      file.exists(coalesce(local_ir_path, ""))
  )
 
cat("KR surveys in data:          ", length(names(all_KR_surveys_with_gps_HR)), "\n")
cat("KR surveys matched to IR:    ", sum(!is.na(ir_lookup$IR_file)), "\n")
cat("IR files found on disk:      ", sum(ir_lookup$ir_exists), "\n")
cat("KR surveys with no IR match: ", sum(is.na(ir_lookup$IR_file)), "\n")
 
# Check for any KR surveys not in ir_lookup at all
missing_from_lookup <- setdiff(
  names(all_KR_surveys_with_gps_HR),
  ir_lookup$survey_name
)
if (length(missing_from_lookup) > 0) {
  cat("KR surveys not found in rdhs metadata:",
      length(missing_from_lookup), "\n")
  print(missing_from_lookup)
}
 
# ================================================================
# 4. Define join keys
#
# IR is ONE ROW PER WOMAN (mother)
# KR is ONE ROW PER CHILD
# Join key: v001 (cluster) + v002 (household) + v003 (woman line)
# This uniquely identifies the mother
# MANY-TO-ONE join: many children → one mother
# Row count must NOT increase after join
#
# Variables added from IR are mother-level characteristics
# (maternal education detail, maternal health, partner/paternal
#  variables, autonomy, media exposure, etc.)
# that are relevant to each child but stored only in IR
# ================================================================
 
ir_join_keys <- c("v001", "v002", "v003")
 
# ================================================================
# 5. Main loop: attach IR variables to each survey
# ================================================================
 
ir_status   <- list()
output_list <- list()
 
for (survey_name in names(all_KR_surveys_with_gps_HR)) {
 
  message("Processing: ", survey_name)
 
  kr_data <- all_KR_surveys_with_gps_HR[[survey_name]]
  
  if (is.null(kr_data) || nrow(kr_data) == 0) next
 
  # ── Step 1: Verify KR is already one row per child ────────────
  n_before    <- nrow(kr_data)
  n_unique_ch <- n_distinct(paste(kr_data$caseid, kr_data$bidx))
 
  if (n_before != n_unique_ch) {
    message("  WARNING: KR has ", n_before - n_unique_ch,
            " duplicate rows — deduplicating to one row per child")
    kr_data <- kr_data %>%
      group_by(caseid, bidx) %>%
      slice(1) %>%
      ungroup()
  }
 
  # ── Step 2: Get KR interview year ─────────────────────────────
  kr_year_str <- if ("v007" %in% names(kr_data)) {
    yr <- unique(as.numeric(zap_labels(kr_data$v007)))
    yr <- ifelse(yr < 100, yr + 1900, yr)
    paste(sort(yr), collapse = ", ")
  } else NA_character_
 
  # ── Step 3: Find IR file ──────────────────────────────────────
  ir_row <- ir_lookup %>% filter(survey_name == !!survey_name)
  ir_zip <- if (nrow(ir_row) == 1 && isTRUE(ir_row$ir_exists)) {
    ir_row$local_ir_path
  } else NA_character_
 
  if (is.na(ir_zip)) {
    message("  No IR file found — keeping KR only")
    ir_status[[survey_name]] <- data.frame(
      survey_name     = survey_name,
      kr_year         = kr_year_str,
      ir_year         = NA_integer_,
      has_ir          = FALSE,
      ir_file         = NA_character_,
      n_children      = nrow(kr_data),
      n_ir_rows       = NA_integer_,
      ir_vars_added   = 0L,
      rows_after_join = nrow(kr_data),
      join_clean      = TRUE
    )
    output_list[[survey_name]] <- kr_data
    next
  }
 
  # ── Step 4: Read IR ───────────────────────────────────────────
  ir_data <- tryCatch(
    rdhs:::read_dhs_flat(ir_zip),
    error = function(e) {
      message("  Failed to read IR: ", e$message)
      NULL
    }
  )
 
  if (is.null(ir_data)) {
    ir_status[[survey_name]] <- data.frame(
      survey_name     = survey_name,
      kr_year         = kr_year_str,
      ir_year         = NA_integer_,
      has_ir          = FALSE,
      ir_file         = basename(ir_zip),
      n_children      = nrow(kr_data),
      n_ir_rows       = NA_integer_,
      ir_vars_added   = 0L,
      rows_after_join = nrow(kr_data),
      join_clean      = TRUE
    )
    output_list[[survey_name]] <- kr_data
    next
  }
 
  # ── Step 5: Extract IR year ───────────────────────────────────
  ir_year <- if ("v007" %in% names(ir_data)) {
    yr <- as.numeric(zap_labels(ir_data$v007))[1]
    as.integer(ifelse(yr < 100, yr + 1900, yr))
  } else NA_integer_
 
  # ── Step 6: Deduplicate IR to ONE ROW PER WOMAN ──────────────
  # IR should already be one row per woman but enforce it
  # Join key: v001 + v002 + v003 uniquely identifies a woman
  ir_data <- ir_data %>%
    mutate(across(all_of(ir_join_keys), as.numeric))
 
  kr_data <- kr_data %>%
    mutate(across(all_of(ir_join_keys), as.numeric))
 
  n_ir_before <- nrow(ir_data)
  ir_data <- ir_data %>%
    group_by(v001, v002, v003) %>%
    slice(1) %>%
    ungroup()
 
  if (nrow(ir_data) != n_ir_before) {
    message("  IR had ", n_ir_before - nrow(ir_data),
            " duplicate woman rows — deduplicated")
  }
 
  # ── Step 7: Select only variables to add from IR ─────────────
  # Keep join keys + variables in variables_to_keep
  # Drop columns already present in KR to avoid conflicts
  kr_existing <- setdiff(names(kr_data), ir_join_keys)
 
  ir_vars_to_add <- ir_data %>%
    dplyr::select(
      all_of(ir_join_keys),
      any_of(variables_to_keep)
    ) %>%
    dplyr::select(-any_of(kr_existing)) %>%
    names() %>%
    setdiff(ir_join_keys)
 
  if (length(ir_vars_to_add) == 0) {
    message("  No new IR variables to add (all already in KR)")
    ir_status[[survey_name]] <- data.frame(
      survey_name     = survey_name,
      kr_year         = kr_year_str,
      ir_year         = ir_year,
      has_ir          = TRUE,
      ir_file         = basename(ir_zip),
      n_children      = nrow(kr_data),
      n_ir_rows       = nrow(ir_data),
      ir_vars_added   = 0L,
      rows_after_join = nrow(kr_data),
      join_clean      = TRUE
    )
    output_list[[survey_name]] <- kr_data
    next
  }
 
  ir_slim <- ir_data %>%
    dplyr::select(all_of(c(ir_join_keys, ir_vars_to_add)))
 
  # ── Step 8: Join IR to KR ─────────────────────────────────────
  # MANY-TO-ONE: many children → one mother
  # Each child row gets its mother's variables added as columns
# Row count MUST stay the same as KR after join
n_kr_before <- nrow(kr_data)

kr_data <- kr_data %>%
  left_join(
    ir_slim,
    by           = ir_join_keys,
    relationship = "many-to-one"
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
          length(ir_vars_to_add), " IR variables added)")
}

# ── Step 9: Final verification ────────────────────────────────
final_n      <- nrow(kr_data)
final_unique <- n_distinct(paste(kr_data$caseid, kr_data$bidx))

if (final_n != final_unique) {
  message("  ⚠️ Still has duplicates after cleanup — investigate")
}

# ── Record status ─────────────────────────────────────────────
ir_status[[survey_name]] <- data.frame(
  survey_name     = survey_name,
  kr_year         = kr_year_str,
  ir_year         = ir_year,
  has_ir          = TRUE,
  ir_file         = basename(ir_zip),
  n_children      = final_n,
  n_ir_rows       = nrow(ir_data),
  ir_vars_added   = length(ir_vars_to_add),
  rows_after_join = n_kr_after,
  join_clean      = join_clean
)

output_list[[survey_name]] <- kr_data

rm(ir_data, ir_slim)
#gc()
}

# ================================================================
# 6. Compile status report
# ================================================================

ir_status_df <- dplyr::bind_rows(ir_status)

# Year match check (same as HR script)
ir_status_df$year_match <- sapply(seq_len(nrow(ir_status_df)), function(i) {
  kr_yrs <- as.numeric(unlist(strsplit(ir_status_df$kr_year[i], ",")))
  ir_yr  <- ir_status_df$ir_year[i]
  if (is.na(ir_yr) || all(is.na(kr_yrs))) "not matched"
  else if (ir_yr %in% kr_yrs) "match"
  else "not matched"
})

cat("\n=== IR MERGE STATUS SUMMARY ===\n")
cat("Surveys processed:           ", nrow(ir_status_df), "\n")
cat("Surveys with IR attached:    ", sum(ir_status_df$has_ir, na.rm = TRUE), "\n")
cat("Surveys without IR:          ", sum(!ir_status_df$has_ir, na.rm = TRUE), "\n")
cat("Clean joins:                 ", sum(ir_status_df$join_clean, na.rm = TRUE), "\n")
cat("Joins that created duplicates:", sum(!ir_status_df$join_clean, na.rm = TRUE), "\n")
cat("Year matches:                ", sum(ir_status_df$year_match == "match", na.rm = TRUE), "\n")
cat("Total IR variables added:    ", sum(ir_status_df$ir_vars_added, na.rm = TRUE), "\n")
cat("Total children across all:   ",
    format(sum(ir_status_df$n_children, na.rm = TRUE), big.mark = ","), "\n")

# ── Final verification: all surveys are one row per child ───────
cat("\n=== FINAL VERIFICATION ===\n")
final_check <- purrr::map_dfr(names(output_list), function(sn) {
  df <- output_list[[sn]]
  data.frame(
    survey    = sn,
    n_rows    = nrow(df),
    n_unique  = n_distinct(paste(df$caseid, df$bidx)),
    clean     = nrow(df) == n_distinct(paste(df$caseid, df$bidx))
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
# 7. Save outputs
# ================================================================

saveRDS(
  output_list,
  "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/all_KR_surveys_with_gps_HR_IR.rds"
)

writexl::write_xlsx(
  list(
    ir_merge_status = ir_status_df,
    final_check     = final_check
  ),
  "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/ir_merge_status.xlsx"
)

message("\n✅ Done. Clean data saved to all_KR_surveys_with_gps_HR_IR.rds")
message("Total surveys: ", length(output_list))
message("Total children: ",
        format(sum(final_check$n_rows), big.mark = ","))


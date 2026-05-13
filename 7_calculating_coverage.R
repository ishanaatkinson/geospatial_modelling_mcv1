# ================================
# Libraries
# ================================
library(dplyr)
library(purrr)
library(writexl)
library(haven)

# ================================
# Load final DHS data
# ================================
final_DHS_data <- readRDS(
  "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/final_DHS_data.rds"
)

output_folder_table <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/TABLES/mcv1_coverage_dhs"
output_folder_figures <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/FIGURES/mcv1_coverage_dhs"

options(survey.lonely.psu = "adjust")

calculate_mcv1_coverage <- function(
    surveys,
    geo_level     = "country",
    disagg_var    = NULL,
    disagg_groups = NULL
) {
  
  # Force lonely PSU handling inside function scope
  options(survey.lonely.psu = "adjust")
  
  geo_vars <- switch(geo_level,
                     "country" = c("CountryName", "SurveyYear", "country_phase"),
                     "admin1"  = c("CountryName", "SurveyYear", "country_phase",
                                   "REGION_CLEAN", "ADM1NAME_CLEAN"),
                     "cluster" = c("CountryName", "SurveyYear", "country_phase",
                                   "REGION_CLEAN", "ADM1NAME_CLEAN", "cluster_psu"),
                     stop("geo_level must be one of: 'country', 'admin1', 'cluster'")
  )
  
  # All possible mcv1 indicator names (used for output consistency)
  all_mcv1_vars <- c("mcv1_any", "mcv1_card_dated", "mcv1_mother_recall",
                     "mcv1_card_marked", "mcv1_on_time", "mcv1_late")
  
  results <- lapply(names(surveys), function(survey_name) {
    
    df <- surveys[[survey_name]]
    if (nrow(df) == 0) return(NULL)
    
    # ── Zap haven labels early ────────────────────────────────────
    df <- df %>%
      mutate(across(where(haven::is.labelled), haven::zap_labels))
    
    # ── Ensure core ID columns exist and are consistent types ─────
    if (!"CountryName"    %in% names(df)) df$CountryName    <- NA_character_
    if (!"SurveyYear"     %in% names(df)) df$SurveyYear     <- NA_character_
    if (!"ADM1NAME_CLEAN" %in% names(df)) df$ADM1NAME_CLEAN <- NA_character_
    if (!"REGION_CLEAN"   %in% names(df)) df$REGION_CLEAN   <- NA
    
    # Force geo/grouping columns to character for bind_rows safety
    df <- df %>%
      mutate(
        CountryName    = as.character(CountryName),
        SurveyYear     = as.character(SurveyYear),
        country_phase  = if ("country_phase" %in% names(.)) as.character(country_phase) else NA_character_,
        REGION_CLEAN   = as.character(REGION_CLEAN),
        ADM1NAME_CLEAN = as.character(ADM1NAME_CLEAN)
      )
    
    # ── Filter out rows with no clean region assignment ───────────
    if (geo_level %in% c("admin1", "cluster")) {
      n_before <- nrow(df)
      df <- df %>% filter(!is.na(ADM1NAME_CLEAN) & !is.na(REGION_CLEAN))
      n_dropped <- n_before - nrow(df)
      if (n_dropped > 0) {
        message(survey_name, ": dropped ", n_dropped, " rows with missing ADM1NAME_CLEAN/REGION_CLEAN")
      }
    }
    
    # ── Filter to 12–23 month window ─────────────────────────────
    if (!"child_age_months" %in% names(df)) {
      message(survey_name, ": child_age_months missing, skipping")
      return(NULL)
    }
    
    df <- df %>%
      mutate(child_age_months = as.numeric(child_age_months)) %>%
      filter(!is.na(child_age_months),
             child_age_months >= 12,
             child_age_months <= 23)
    
    if (nrow(df) == 0) return(NULL)
    
    # ── MCV1 vaccination source variable ─────────────────────────
    if ("child_vacc_measles" %in% names(df)) {
      df <- df %>% mutate(vacc = as.numeric(child_vacc_measles))
    } else if ("h9" %in% names(df)) {
      df <- df %>%
        mutate(vacc = as.numeric(h9),
               vacc = ifelse(vacc %in% c(8, 9), NA, vacc))
    } else {
      message(survey_name, ": no vaccination variable found, skipping")
      return(NULL)
    }
    
    # 🔴 CRITICAL FIX: restrict to valid DHS codes ONLY
    df <- df %>%
      mutate(vacc = ifelse(vacc %in% c(0,1,2,3), vacc, NA_real_))
    
    # ── MCV1 binary indicators ───────────────────────────────────
    df <- df %>%
      mutate(
        mcv1_any = case_when(
          vacc %in% c(1,2,3) ~ 1L,
          vacc == 0          ~ 0L,
          TRUE               ~ NA_integer_
        ),
        
        mcv1_card_dated = case_when(
          vacc == 1          ~ 1L,
          vacc %in% c(0,2,3) ~ 0L,
          TRUE               ~ NA_integer_
        ),
        
        mcv1_mother_recall = case_when(
          vacc == 2          ~ 1L,
          vacc %in% c(0,1,3) ~ 0L,
          TRUE               ~ NA_integer_
        ),
        
        mcv1_card_marked = case_when(
          vacc == 3          ~ 1L,
          vacc %in% c(0,1,2) ~ 0L,
          TRUE               ~ NA_integer_
        )
      )
    
    # ── Timing variables ─────────────────────────────────────────
    has_timing <- "child_vacc_measles_months" %in% names(df)
    
    if (has_timing) {
      
      df <- df %>%
        mutate(child_vacc_measles_months = as.numeric(child_vacc_measles_months)) %>%
        mutate(
          mcv1_on_time = case_when(
            vacc %in% c(1,3) & !is.na(child_vacc_measles_months) &
              child_vacc_measles_months <= 12 ~ 1L,
            vacc %in% c(1,3) & !is.na(child_vacc_measles_months) &
              child_vacc_measles_months > 12  ~ 0L,
            vacc == 0                         ~ 0L,
            vacc == 2                         ~ NA_integer_,
            TRUE                              ~ NA_integer_
          ),
          
          mcv1_late = case_when(
            vacc %in% c(1,3) & !is.na(child_vacc_measles_months) &
              child_vacc_measles_months > 12  ~ 1L,
            vacc %in% c(1,3) & !is.na(child_vacc_measles_months) &
              child_vacc_measles_months <= 12 ~ 0L,
            vacc == 0                         ~ 0L,
            vacc == 2                         ~ NA_integer_,
            TRUE                              ~ NA_integer_
          )
        )
      
    } else {
      df$mcv1_on_time <- NA_integer_
      df$mcv1_late    <- NA_integer_
    }
    
    # ── Disaggregation ────────────────────────────────────────────
    if (!is.null(disagg_var) && disagg_var %in% names(df)) {
      df <- df %>% mutate(disagg_raw = as.numeric(.data[[disagg_var]]))
      if (!is.null(disagg_groups)) {
        df <- df %>%
          mutate(disagg = case_when(
            !!!purrr::imap(disagg_groups,
                           ~ rlang::expr(disagg_raw %in% !!.x ~ !!.y)),
            TRUE ~ NA_character_
          ))
      } else {
        df <- df %>% mutate(disagg = as.character(disagg_raw))
      }
    } else {
      df <- df %>% mutate(disagg = "all")
    }
    
    # Drop rows with missing disagg
    df <- df %>% filter(!is.na(disagg))
    if (nrow(df) == 0) return(NULL)
    
    # ── Ensure geo vars exist ─────────────────────────────────────
    for (gv in geo_vars) {
      if (!gv %in% names(df)) df[[gv]] <- NA_character_
    }
    
    # ── Force all geo vars to character ────────────────────────────
    df <- df %>%
      mutate(across(all_of(geo_vars), as.character))
    
    # ── Build survey design variables ─────────────────────────────
    df$weight <- as.numeric(df$v005) / 1e6
    
    use_nest <- TRUE
    
    if ("v023" %in% names(df) && !all(is.na(df$v023))) {
      df$strata_var <- as.numeric(df$v023)
    } else if (all(c("v024", "v025") %in% names(df))) {
      df$strata_var <- as.integer(interaction(df$v024, df$v025, drop = TRUE))
    } else {
      df$strata_var <- rep(1L, nrow(df))
      use_nest <- FALSE
    }
    
    df$psu_var <- if ("v021" %in% names(df) && !all(is.na(df$v021))) {
      as.numeric(df$v021)
    } else if ("cluster_psu" %in% names(df) && !all(is.na(df$cluster_psu))) {
      as.numeric(df$cluster_psu)
    } else if ("v001" %in% names(df)) {
      as.numeric(df$v001)
    } else {
      seq_len(nrow(df))
    }
    
    # ── Select ONLY the columns needed for survey analysis ────────
    keep_cols <- unique(c(
      geo_vars, "disagg",
      "mcv1_any", "mcv1_card_dated", "mcv1_mother_recall",
      "mcv1_card_marked", "mcv1_on_time", "mcv1_late",
      "weight", "strata_var", "psu_var"
    ))
    
    df_design <- df %>% select(all_of(keep_cols))
    
    # ── Unweighted counts ─────────────────────────────────────────
    denom <- df %>%
      group_by(across(all_of(c(geo_vars, "disagg")))) %>%
      summarise(
        n_total_12_23mo   = n(),
        n_vaccinated_any  = sum(mcv1_any       == 1, na.rm = TRUE),
        n_card_dated      = sum(mcv1_card_dated == 1, na.rm = TRUE),
        n_mother_recall   = sum(mcv1_mother_recall == 1, na.rm = TRUE),
        n_card_marked     = sum(mcv1_card_marked == 1, na.rm = TRUE),
        n_on_time         = sum(mcv1_on_time    == 1, na.rm = TRUE),
        n_late            = sum(mcv1_late        == 1, na.rm = TRUE),
        n_unknown_timing  = sum(is.na(mcv1_on_time) & mcv1_any == 1,
                                na.rm = TRUE),
        .groups = "drop"
      )
    
    # ── Survey design on clean subset ─────────────────────────────
    design <- tryCatch(
      survey::svydesign(
        ids     = ~psu_var,
        strata  = ~strata_var,
        weights = ~weight,
        data    = df_design,
        nest    = use_nest
      ),
      error = function(e) {
        message(survey_name, ": survey design error - ", e$message)
        NULL
      }
    )
    
    if (is.null(design)) return(NULL)
    
    # ── Build dynamic formula (exclude all-NA variables) ──────────
    mcv1_valid <- all_mcv1_vars[
      sapply(all_mcv1_vars, function(v) any(!is.na(df_design[[v]])))
    ]
    
    if (length(mcv1_valid) == 0) {
      message(survey_name, ": no non-NA mcv1 variables, skipping")
      return(NULL)
    }
    
    # ── Weighted coverage proportions ─────────────────────────────
    group_vars    <- c(geo_vars, "disagg")
    group_formula <- as.formula(paste("~", paste(group_vars, collapse = " + ")))
    
    coverage <- tryCatch({
      
      coverage_list <- lapply(mcv1_valid, function(var) {
        
        f <- as.formula(paste0("~", var))
        
        res <- survey::svyby(
          formula  = f,
          by       = group_formula,
          design   = design,
          FUN      = survey::svymean,
          na.rm    = TRUE,
          keep.var = TRUE
        )
        
        names(res)[names(res) == var] <- var
        names(res)[names(res) == paste0("se.", var)] <- paste0("se_", var)
        
        res
      })
      
      Reduce(function(x, y) {
        dplyr::left_join(x, y, by = group_vars)
      }, coverage_list)
      
    }, error = function(e) {
      message(survey_name, ": svyby error - ", e$message)
      NULL
    })
    
    if (is.null(coverage)) return(NULL)
    
    # ── Add back any excluded variables as NA columns ─────────────
    mcv1_missing <- setdiff(all_mcv1_vars, mcv1_valid)
    for (mv in mcv1_missing) {
      coverage[[mv]]              <- NA_real_
      coverage[[paste0("se_", mv)]] <- NA_real_
    }
    
    # ── Combine weighted estimates with unweighted counts ─────────
    coverage %>%
      left_join(denom, by = c(geo_vars, "disagg")) %>%
      rename_with(~ gsub("^se\\.", "se_", .x)) %>%
      mutate(
        survey_name = survey_name,
        geo_level   = geo_level,
        disagg_var  = ifelse(is.null(disagg_var), "none", disagg_var),
        across(any_of(c("mcv1_any", "mcv1_card_dated", "mcv1_mother_recall",
                        "mcv1_card_marked", "mcv1_on_time", "mcv1_late")),
               ~ round(.x * 100, 1)),
        across(starts_with("se_"),
               ~ round(.x * 100, 1))
      )
    
  }) %>% dplyr::bind_rows()
  
  # ── Save ──────────────────────────────────────────────────────
  disagg_label <- ifelse(is.null(disagg_var), "all", disagg_var)
  out_path <- file.path(
    output_folder_table,
    paste0("mcv1_coverage_", geo_level, "_by_", disagg_label, ".xlsx")
  )
  writexl::write_xlsx(results, out_path)
  message("Saved: ", out_path)
  
  return(results)
}



# ── Example calls ────────────────────────────────────────────────

national_all <- calculate_mcv1_coverage(
  final_DHS_data)

# Overall adm1 coverage
admin1_all <- calculate_mcv1_coverage(
  final_DHS_data,
  geo_level = "admin1"
)

cluster_all <- calculate_mcv1_coverage(
  final_DHS_data, 
  geo_level = "cluster"
)

# National by wealth quintile
admin1_wealth <- calculate_mcv1_coverage(
  final_DHS_data,
  geo_level  = "admin1",
  disagg_var = "wealth_index_quintile",
  disagg_groups = list(
    "Poorest"  = 1,
    "Second"   = 2,
    "Middle"   = 3,
    "Fourth"   = 4,
    "Richest"  = 5
  )
)

# Admin1 (region) by maternal education level
admin1_edu <- calculate_mcv1_coverage(
  final_DHS_data,
  geo_level  = "admin1",
  disagg_var = "respondent_edu_level",
  disagg_groups = list(
    "No education" = 0,
    "Primary"      = 1,
    "Secondary"    = 2,
    "Higher"       = 3
  )
)

# Admin1 by urban/rural
admin1_urban <- calculate_mcv1_coverage(
  final_DHS_data,
  geo_level  = "admin1",
  disagg_var = "urban_rural",
  disagg_groups = list(
    "Urban" = 1,
    "Rural" = 2
  )
)

# By mother's age group
admin1_age_group <- calculate_mcv1_coverage(
  final_DHS_data,
  geo_level  = "admin1",
  disagg_var = "respondent_age_group",
  disagg_groups = list(
    "15-19" = 1,
    "20-24" = 2,
    "25-29" = 3,
    "30-34" = 4,
    "35-39" = 5,
    "40-44" = 6,
    "45-49" = 7
  )
)


# By birth order
admin1_birth_order <- calculate_mcv1_coverage(
  final_DHS_data,
  geo_level  = "admin1",
  disagg_var = "birth_order",
  disagg_groups = list(
    "1st"      = 1,
    "2nd-3rd"  = 2,
    "4th-6th"  = 3,
    "7th+"     = c(4, 5, 6)
  )
)

# admin1 level by healthcare cost barrier
admin1_cost_barrier <- calculate_mcv1_coverage(
  final_DHS_data,
  geo_level  = "admin1",
  disagg_var = "healthcare_barriers_cost",
  disagg_groups = list(
    "No cost barrier" = 0,
    "Cost barrier"    = 1
  )
)


# By child sex
admin1_child_sex <- calculate_mcv1_coverage(
  final_DHS_data,
  geo_level  = "admin1",
  disagg_var = "child_sex",
  disagg_groups = list(
    "Male"   = 1,
    "Female" = 2
  )
)

# By health card
admin1_health_card <- calculate_mcv1_coverage(
  final_DHS_data,
  geo_level  = "admin1",
  disagg_var = "child_health_card",
  disagg_groups = list(
    "No card"  = 0,
    "Has card" = 1
  )
)

# By ANC visits
admin1_anc <- calculate_mcv1_coverage(
  final_DHS_data,
  geo_level  = "admin1",
  disagg_var = "anc_visits_number",
  disagg_groups = list(
    "0 visits"   = 0,
    "1-3 visits" = 1,
    "4+ visits"  = c(2, 3, 4, 5, 6)
  )
)

# By place of delivery
admin1_delivery <- calculate_mcv1_coverage(
  final_DHS_data,
  geo_level  = "admin1",
  disagg_var = "place_of_delivery",
  disagg_groups = list(
    "Home"     = 1,
    "Facility" = c(2, 3)
  )
)

# By healthcare transport barrier
admin1_transport_barrier <- calculate_mcv1_coverage(
  final_DHS_data,
  geo_level  = "admin1",
  disagg_var = "healthcare_barriers_transport",
  disagg_groups = list(
    "No transport barrier" = 0,
    "Transport barrier"    = 1
  )
)

# By any healthcare barrier
admin1_barriers <- calculate_mcv1_coverage(
  final_DHS_data,
  geo_level  = "admin1",
  disagg_var = "healthcare_barriers",
  disagg_groups = list(
    "No barriers"        = 0,
    "At least 1 barrier" = 1
  )
)


# By postnatal care
admin1_pnc <- calculate_mcv1_coverage(
  final_DHS_data,
  geo_level  = "admin1",
  disagg_var = "postnatal_care_newborn",
  disagg_groups = list(
    "No PNC"  = 0,
    "Had PNC" = 1
  )
)



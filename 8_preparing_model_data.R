
# ================================================================
# LIBRARIES
# ================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(spdep)
library(gstat)
library(car)
library(corrplot)
library(patchwork)
library(purrr)
library(writexl)
library(pROC)
library(grid)
library(gridExtra)

inla_available <- requireNamespace("INLA", quietly = TRUE)
if (inla_available) library(INLA)
if (!inla_available) stop("R-INLA is required for this analysis. Install with:\n",
                          '  install.packages("INLA", repos=c(INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)')


# ================================================================
# 0. LOAD DATA
# ================================================================

base_path <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT"

final_DHS_data <- readRDS(
  file.path(base_path, "Data/DHS_final/final_DHS_data.rds")
)

output_folder <- file.path(base_path, "Output/Tables")
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)


# ================================================================
# 1. DEFINE VARIABLES — EXACTLY MATCHING UTAZI TABLE A
# ================================================================
#
# Every variable is treated as a FACTOR in the model, because the
# coding is categorical (e.g. birth_order: 0=1-2, 1=3-5, 2=6+).
# A "1 unit increase in birth_order" is meaningless — what matters
# is the odds ratio for category 1 vs category 0 (reference).
#
# We define three lists corresponding to the three levels of
# Utazi Equation (1): individual, household, community.
# The reference category for each variable is the code listed first
# (usually 0), chosen per Utazi: "the category that had the least
# likelihood of being vaccinated based on prior knowledge."
# ================================================================

# ────────────────────────────────────────────────────────────────
# VARIABLE TYPE KEY:
#   [F] = Factor (categorical) — treated with factor() in R so the
#         model estimates separate ORs per category vs reference.
#   [C] = Continuous — kept numeric; the model estimates a single
#         slope (OR per unit increase). Used only where a unit
#         increase IS meaningful (e.g. years of education, minutes
#         of travel time, wealth z-score, autonomy score 0-6).
#   [B] = Binary (special case of factor with 2 levels)
#
# REDUNDANCY NOTE: Some variables measure the same construct at
# different granularities (e.g. respondent_edu_level vs
# respondent_edu_attainment vs respondent_edu_years). All are
# included as CANDIDATES; the multicollinearity check (GVIF) will
# drop the redundant ones automatically per Utazi's method. This
# is exactly what Utazi did: "We chose between similar covariates
# where both were available for a given country, mostly based on
# completeness of data."
# ────────────────────────────────────────────────────────────────

# ── INDIVIDUAL-LEVEL COVARIATES (x_p,ijkl) ───────────────────────
#    These vary at the child level within households.
individual_vars_factor <- c(
  "child_sex",                     # [F] 0=Male(ref), 1=Female
  "birth_order",                   # [F] 0=1-2(ref), 1=3-5, 2=6+
  "skilled_birth_attendance",      # [B] 0=No(ref), 1=Yes
  "birth_quarter",                 # [F] 0=Jan-Mar(ref), 1=Apr-Jun, 2=Jul-Sep, 3=Oct-Dec
  "anc_visits_number",             # [F] 0=None(ref), 1=1-3, 2=4+
  "anc_first_visit_timing",        # [F] 0=0mo(ref), 1=1-2, 2=3-4, 3=5-6, 4=7-8, 5=9+
  "anc_provider",                  # [F] 1=Skilled, 2=TBA, 3=Relative, 4=None(ref)
  "postnatal_care_newborn",        # [B] 0=No(ref), 1=Yes
  "postnatal_care_newborn_timing", # [B] 0=No/DK(ref), 1=Yes
  "postnatal_care_newborn_person", # [F] 1=Skilled, 2=Traditional/CHW
  "child_health_card",             # [B] 0=No(ref), 1=Yes
  "child_vitamin_a",               # [B] 0=No(ref), 1=Yes
  "place_of_delivery",             # [F] 1=Home(ref), 2=Government, 3=Private
  "birth_interval_preceding",      # [F] 1=9-15mo(ref), 2=16-22, 3=23-29, 4=30-36, 5=37-43, 6=44-50, 7=51+
  "vacc_campaign"                  # [B] 0=No(ref), 1=Yes
)

# ── HOUSEHOLD-LEVEL COVARIATES (x_p,jkl) ─────────────────────────
#    These are constant within a household but vary between households
#    in the same cluster. Includes maternal demographics, SES, assets,
#    media, barriers, and autonomy.

household_vars_factor <- c(
  # Maternal demographics
  "respondent_age",                # [F] 0=15-19(ref), 1=20-29, 2=30-39, 3=40-49
  "marital_status",                # [F] 0=Never(ref), 1=Married, 2=Div/Wid/Sep
  "household_size",                # [F] 0=Large≥9(ref), 1=Medium5-8, 2=Small≤4
  "hh_head_sex",                   # [B] 0=Female(ref), 1=Male
  "hh_head_age",                   # [F] 1=<20(ref), 2=20-29, 3=30-39, ... 8=80+
  "number_living_children",        # [F] 0=0(ref), 1=1-3, 2=4-6, 3=7-9, 4=10-12, 5=13+
  "total_children_ever_born",      # [F] 0=0(ref), 1=1-3, 2=4-6, 3=7-9, 4=10-12, 5=13+
  "length_of_stay",                # [F] 0=<1yr(ref), 1=1-3, 2=4-5, 3=>5/always
  
  # Education (multiple granularities — GVIF will select)
  "respondent_edu_level",          # [F] 0=None(ref), 1=Primary, 2=Secondary, 3=Higher
  "respondent_edu_attainment",     # [F] 0=None(ref), 1=Inc.primary, 2=Comp.primary, 3=Inc.sec, 4=Comp.sec, 5=Higher
  "respondent_literacy",           # [F] 0=Cannot read(ref), 1=Parts, 2=Whole, 3=No card, 4=Blind
  "partner_edu_level",             # [F] 0=None(ref), 1=Primary, 2=Secondary, 3=Higher
  "partner_edu_attainment",        # [F] 0=None(ref)...5=Higher
  
  # Employment / Occupation
  "respondent_working",            # [B] 0=No(ref), 1=Yes
  "respondent_working_12mo",       # [B] 0=No(ref), 1=Yes
  "respondent_occupation",         # [F] 0=Agric(ref), 1=Clerical/Sales, 2=Prof, 3=Manual
  #"partner_occupation",            # [F] 0=Agric(ref), 1=Clerical/Sales, 2=Prof, 3=Manual
  
  # Wealth
  "wealth_index_quintile",         # [F] 0=Poor(ref), 1=Middle, 2=Rich
  
  # Assets
  # "hh_electricity",                # [B] 0=No(ref), 1=Yes # REMOVED DUE TO BEING IN WEALTH
  # "hh_radio",                      # [B] 0=No(ref), 1=Yes # REMOVED DUE TO BEING IN WEALTH
  # "hh_television",                 # [B] 0=No(ref), 1=Yes # REMOVED DUE TO BEING IN WEALTH
  # "hh_mobile_phone",               # [B] 0=No(ref), 1=Yes # REMOVED DUE TO BEING IN WEALTH
  # "hh_bicycle",                    # [B] 0=No(ref), 1=Yes # REMOVED DUE TO BEING IN WEALTH
  # "hh_motorcycle",                 # [B] 0=No(ref), 1=Yes # REMOVED DUE TO BEING IN WEALTH
  # "hh_car",                        # [B] 0=No(ref), 1=Yes # REMOVED DUE TO BEING IN WEALTH
  # "hh_computer",                   # [B] 0=No(ref), 1=Yes # REMOVED DUE TO BEING IN WEALTH
  "hh_slum_dwelling",              # [B] 0=Non-slum(ref), 1=Slum (UN-Habitat)
  
  # Healthcare barriers
  "healthcare_barriers",           # [B] 0=Had problem(ref), 1=No problem
  "healthcare_barriers_cost",      # [B] 0=Had problem(ref), 1=No problem
  "healthcare_barriers_transport", # [B] 0=Had problem(ref), 1=No problem
  "healthcare_barriers_safety",    # [B] 0=Had problem(ref), 1=No problem
  "healthcare_barriers_stockout",  # [B] 0=Had problem(ref), 1=No problem
  
  # Media / communication
  "media_exposure",                # [B] 0=No(ref), 1=Yes (any weekly)
  "media_radio",                   # [F] 0=Never(ref), 1=<weekly, 2=weekly, 3=daily
  "media_television",              # [F] 0=Never(ref), 1=<weekly, 2=weekly, 3=daily
  "media_newspaper",               # [F] 0=Never(ref), 1=<weekly, 2=weekly, 3=daily
  "mobile_internet_use",           # [B] 0=No(ref), 1=Yes
  "internet_use",                  # [B] 0=No(ref), 1=Yes
  
  # Other household
  "health_insurance",              # [B] 0=No(ref), 1=Yes
  "hh_mosquito_net",               # [B] 0=No(ref), 1=Yes
  # "mother_owns_land_or_house",     # [B] 0=No(ref), 1=Yes   # REMOVED DUE TO BEING IN WEALTH
  # "mother_bank_account",           # [B] 0=No(ref), 1=Yes   # REMOVED DUE TO BEING IN WEALTH
  "mother_knows_malaria"           # [B] 0=No(ref), 1=Yes
)

# Variables treated as CONTINUOUS (a unit increase IS meaningful)
household_vars_continuous <- c(
  "respondent_edu_years",          # [C] raw years of education (0-20+)
  "partner_edu_years",             # [C] raw years partner education
  "respondent_edu_grade",          # [C] highest grade at level
  "wealth_index_score",            # [C] continuous DHS wealth z-score
  "womens_autonomy_has_a_say",     # [C] 0-6 count of decision domains
  "womens_autonomy_no_say"         # [C] 0-6 count of decision domains
)

# ── COMMUNITY/CLUSTER-LEVEL COVARIATES (x_p,kl) ─────────────────
# These are constant within a cluster.
# The "key community variables" (travel time, conflict, slum) are
# evaluated separately via PCV (Utazi's method).
community_vars_factor <- c(
  "urban_rural",                   # [B] 0=Rural(ref), 1=Urban
  "conflict_area",                 # [B] 0=No conflict(ref), 1=Conflict
  "urban_slum"                     # [B] 0=Not slum(ref), 1=Slum
)

# Travel time variables — will be TERTILED per country (Utazi Table B)
# so that "Higher" (most remote) is the reference category.
travel_time_vars_to_tertile <- c(
  "travel_time_to_city",           # [C→F] minutes to city ≥50k → tertiled
  "travel_time_to_HC_motor",       # [C→F] minutes to HC by motorbike → tertiled
  "travel_time_to_HC_walk"         # [C→F] minutes to HC by foot → tertiled
)

gc_vars_base <- c(
  "NIGHTLIGHTS_COMPOSITE",         # Nighttime lights (proxy urbanicity/development; Utazi 2018)
  "ELEVATION",                      # Elevation in metres (access barrier)
  "ARIDITY",                        # Aridity index (climate stress)
  "PRECIPITATION",                  # Annual precipitation (agricultural productivity)
  "MEAN_TEMPERATURE",               # Mean annual temperature
  "ENHANCED_VEGETATION_INDEX",      # EVI (vegetation/agricultural activity)
  "UN_POPULATION_DENSITY",          # Population density (WorldPop/UN)
  "TRAVEL_TIMES",                   # Travel time to nearest city (Weiss et al. 2018)
  "GLOBAL_HUMAN_FOOTPRINT",         # Human footprint index
  "MALARIA_PREVALENCE",             # P. falciparum prevalence (MAP)
  "DROUGHT_EPISODES"               # Drought episodes
  
)

# Religion and ethnicity — country-specific factors
# The largest category is used as reference (Utazi: "the largest
# category was used when the smallest yielded imprecise estimates")
#factor_vars_country_specific <- c("religion", "ethnicity")

# ── Combine all candidate lists ──────────────────────────────────
all_factor_vars    <- c(individual_vars_factor, household_vars_factor, community_vars_factor)
gc_continuous_vars <- paste0("gc_", tolower(gc_vars_base))
all_continuous_vars <- c(household_vars_continuous, gc_continuous_vars)
all_model_vars     <- c(all_factor_vars, all_continuous_vars)

cat("Full candidate variable set:\n")
cat("  Individual-level (factor):", length(individual_vars_factor), "\n")
cat("  Household-level (factor):", length(household_vars_factor), "\n")
cat("  Household-level (continuous):", length(all_continuous_vars), "\n")
cat("  Community-level (factor):", length(community_vars_factor), "\n")
cat("  Travel time (to tertile):", length(travel_time_vars_to_tertile), "\n")
cat("  Geospatial covariates:", length(gc_continuous_vars), "\n")
cat("  TOTAL CANDIDATES:", length(all_model_vars) + length(travel_time_vars_to_tertile) +
      length(gc_continuous_vars), "\n")


# ================================================================
# 2. PREPARE INDIVIDUAL-LEVEL ANALYSIS DATASET
# ================================================================

cat("\n--- Preparing individual-level data ---\n")

all_surveys <- bind_rows(final_DHS_data, .id = "survey_idx") 

data.table::fwrite(all_surveys, paste0(base_path, "/Data/DHS_final/all_surveys.csv"))

eligible <- bind_rows(final_DHS_data, .id = "survey_idx") %>%
  filter(
    !is.na(child_age_months), child_age_months >= 12, child_age_months <= 23,
    !is.na(child_vacc_measles),
    !is.na(cluster_psu),
    !is.na(LATNUM), LATNUM != 0,
    !is.na(LONGNUM), LONGNUM != 0
  )

data.table::fwrite(eligible, paste0(base_path, "/Data/DHS_final/model_data.csv"))


all_surveys <- vroom::vroom(paste0(base_path, "/Data/DHS_final/model_data.csv"))

eligible <- all_surveys %>%
  filter(
    !is.na(child_age_months), child_age_months >= 12, child_age_months <= 23,
    !is.na(child_vacc_measles),
    !is.na(cluster_psu),
    !is.na(LATNUM), LATNUM != 0,
    !is.na(LONGNUM), LONGNUM != 0
  ) %>%
  mutate(
    # Binary outcome: MCV1 received (any source)
    mcv1 = as.integer(child_vacc_measles %in% c(1, 2, 3)),
    
    # Stratum identifier (for the stratum-level random effect)
    # DHS strata are region × urban/rural (v023 or v022)
    stratum_id = ifelse("v023" %in% names(.) & !is.na(v023),
                        as.character(v023),
                        paste(CountryName, region, urban_rural, sep = "_")),
    
    # Household identifier (for household-level random effect)
    hh_id = paste(CountryName, cluster_psu,
                  ifelse("v002" %in% names(.), v002,
                         ifelse("hhid" %in% names(.), hhid, row_number())),
                  sep = "_"),
    
    # Cluster identifier (for community-level random effect)
    cluster_id = paste(CountryName, cluster_psu, sep = "_")
  )

data.table::fwrite(eligible, paste0(base_path, "/Data/DHS_final/model_data_v2.csv"))

eligible <- vroom::vroom(paste0(base_path, "/Data/DHS_final/model_data_v2.csv"))


# ── Pick year-matched GC columns ─────────────────────────────────
pick_gc_year <- function(df, base_name) {
  gc_cols <- grep(paste0("^", base_name, "_\\d{4}$"), names(df), value = TRUE)
  if (length(gc_cols) == 0) {
    if (base_name %in% names(df)) return(as.numeric(df[[base_name]]))
    return(rep(NA_real_, nrow(df)))
  }
  gc_years <- as.numeric(gsub(paste0(base_name, "_"), "", gc_cols))
  sy <- if ("v007" %in% names(df)) as.numeric(df$v007) else rep(2015, nrow(df))
  result <- rep(NA_real_, nrow(df))
  best_d <- rep(Inf, nrow(df))
  for (i in seq_along(gc_years)) {
    vals <- as.numeric(df[[gc_cols[i]]])
    d <- abs(sy - gc_years[i])
    closer <- !is.na(d) & d < best_d
    result[closer] <- vals[closer]
    best_d[closer] <- d[closer]
  }
  result
}

gc_col_names <- paste0("gc_", tolower(gc_vars_base))
for (gc in gc_vars_base) {
  eligible[[paste0("gc_", tolower(gc))]] <- pick_gc_year(eligible, gc)
}

cat("GC variables added:", sum(gc_col_names %in% names(eligible)), "\n")

cat("Eligible children:", nrow(eligible), "\n")
cat("Countries:", n_distinct(eligible$CountryName), "\n")
cat("Clusters:", n_distinct(eligible$cluster_id), "\n")
cat("Households:", n_distinct(eligible$hh_id), "\n")
cat("Strata:", n_distinct(eligible$stratum_id), "\n")
cat("MCV1 coverage:", round(mean(eligible$mcv1) * 100, 1), "%\n")


data.table::fwrite(eligible, paste0(base_path, "/Data/DHS_final/model_data_v3.csv"))

eligible <- vroom::vroom(paste0(base_path, "/Data/DHS_final/model_data_v3.csv"))


# ================================================================
# 3. CONVERT VARIABLES TO FACTORS
# ================================================================
#
# THIS IS CRITICAL. All ordinal/categorical variables must be factors
# so that the model estimates odds ratios PER CATEGORY, not a single
# linear slope. For example:
#
#   wealth_index_quintile as numeric: β = 0.5 means "a 1-unit increase
#     in wealth is associated with..." — but moving from 0 (poor) to
#     1 (middle) may have a very different effect than 1 to 2 (rich).
#
#   wealth_index_quintile as factor: the model estimates separate
#     β_middle and β_rich, giving OR_middle vs poor and OR_rich vs poor.
#
# The reference category is always level "0" (the first level of the
# factor), matching Utazi's Table A reference categories.
# ================================================================

cat("\n--- Converting variables to factors ---\n")

# Convert FACTOR variables with "0" as reference
for (v in all_factor_vars) {
  if (v %in% names(eligible)) {
    eligible[[v]] <- factor(eligible[[v]])
    if ("0" %in% levels(eligible[[v]])) {
      eligible[[v]] <- relevel(eligible[[v]], ref = "0")
    }
    # Special cases: anc_provider ref=4 (None), place_of_delivery ref=1 (Home)
    if (v == "anc_provider" && "4" %in% levels(eligible[[v]]))
      eligible[[v]] <- relevel(eligible[[v]], ref = "4")
    if (v == "place_of_delivery" && "1" %in% levels(eligible[[v]]))
      eligible[[v]] <- relevel(eligible[[v]], ref = "1")
    if (v == "birth_interval_preceding" && "1" %in% levels(eligible[[v]]))
      eligible[[v]] <- relevel(eligible[[v]], ref = "1")
    if (v == "hh_head_age" && "1" %in% levels(eligible[[v]]))
      eligible[[v]] <- relevel(eligible[[v]], ref = "1")
  }
}

# Continuous variables: ensure numeric (no factor conversion)
for (v in all_continuous_vars) {
  if (v %in% names(eligible)) {
    eligible[[v]] <- as.numeric(eligible[[v]])
  }
}

cat("Factor variables converted:", sum(all_factor_vars %in% names(eligible)), "\n")
cat("Continuous variables kept numeric:", sum(all_continuous_vars %in% names(eligible)), "\n")



# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  CLUSTER-LEVEL UTAZI MODEL                                               ║
# ║  Binomial outcome: y_c / n_c  (vaccinated / eligible per cluster)        ║
# ║  Covariates aggregated to cluster means / proportions                    ║
# ║  Methods kept as close to Utazi (2020) as possible                       ║
# ╚════════════════════════════════════════════════════════════════════════════╝


# ════════════════════════════════════════════════════════════════════════════
# PRE-LOOP: AGGREGATE INDIVIDUAL DATA TO CLUSTER LEVEL
# ════════════════════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════════════════════\n")
cat("  AGGREGATING INDIVIDUAL DATA TO CLUSTER LEVEL\n")
cat("══════════════════════════════════════════════════════════\n\n")

# ── Identify variable types in the eligible dataset ──────────────
# Factor / categorical variables → will become proportions per level
# Numeric / continuous variables → will become cluster means

all_candidate_vars <- intersect(all_model_vars, names(eligible))

factor_vars_all <- all_candidate_vars[sapply(eligible[all_candidate_vars], function(x) {
  is.factor(x) || is.character(x) || is.logical(x)
})]
numeric_vars_all <- all_candidate_vars[sapply(eligible[all_candidate_vars], is.numeric)]

cat("Candidate variables:", length(all_candidate_vars), "\n")
cat("  Factor/categorical:", length(factor_vars_all), "\n")
cat("  Numeric/continuous:", length(numeric_vars_all), "\n")

# ── Build cluster-level dataset ──────────────────────────────────
# For each cluster we compute:
#   - y: number vaccinated (sum of mcv1)
#   - n_children: number of eligible children
#   - coverage: y / n_children
#   - For each numeric covariate: cluster mean
#   - For each factor covariate: proportion in EACH level
#     (one column per non-reference level, named var_prop_levelname)

# Step A: numeric means per cluster (including mcv1 outcome)
cluster_numeric <- eligible %>%
  group_by(cluster_id, CountryName, LATNUM, LONGNUM, stratum_id) %>%
  summarise(
    y          = sum(mcv1, na.rm = TRUE),
    n_children = n(),
    coverage   = mean(mcv1, na.rm = TRUE),
    across(all_of(numeric_vars_all), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
    .groups = "drop"
  )

cat("Clusters after numeric aggregation:", nrow(cluster_numeric), "\n")

# Step B: factor proportions per cluster
#   For each factor variable, compute the proportion of children
#   in each level within each cluster.
#   This creates columns like:  mother_edu_prop_primary, mother_edu_prop_secondary, etc.
#   The reference level is excluded (its proportion = 1 - sum of others).

prop_cols_map <- list()   # track which proportion columns came from which factor

cluster_props <- cluster_numeric  # start with numeric aggregation

for (fv in factor_vars_all) {
  x <- eligible[[fv]]
  if (is.character(x) || is.logical(x)) x <- factor(x)
  
  levs <- levels(droplevels(x[!is.na(x)]))
  if (length(levs) < 2) next   # skip constant factors
  
  # Reference level = the most common category (following Utazi)
  tab <- sort(table(x, useNA = "no"), decreasing = TRUE)
  ref_lev <- names(tab)[1]
  
  # Compute proportion for ALL levels (including reference)
  # At cluster level, all proportions carry information — unlike
  # individual-level dummies where the reference is implicit.
  # Note: the proportions sum to 1 within each cluster, so one will
  # be dropped during GVIF/collinearity screening (Step 3).
  prop_names <- c()
  
  for (lv in levs) {
    new_col <- paste0(fv, "_prop_", make.names(lv))
    prop_names <- c(prop_names, new_col)
    
    lv_props <- eligible %>%
      group_by(cluster_id) %>%
      summarise(!!new_col := mean(.data[[fv]] == lv, na.rm = TRUE),
                .groups = "drop")
    
    cluster_props <- left_join(cluster_props, lv_props, by = "cluster_id")
  }
  
  prop_cols_map[[fv]] <- list(
    ref_level  = ref_lev,
    prop_cols  = prop_names,
    all_levels = levs
  )
  
  cat("  Factor:", fv, "→ ref =", ref_lev,
      "| created", length(prop_names), "proportion columns (incl. reference)\n")
}

# ── Also handle travel-time tertile variables ────────────────────
# These need to be tertiled WITHIN each country, so we do the
# tertiling inside the loop. But we record which raw travel-time
# variables exist so the loop knows to process them.

cat("\nTravel-time variables to tertile inside loop:\n")
cat(paste("  ", travel_time_vars_to_tertile), sep = "\n")

# ── Final cluster-level dataset ──────────────────────────────────
cluster_data <- cluster_props

writexl::write_xlsx(cluster_data, paste0(base_path, "/Data/DHS_final/cluster_data.xlsx"))


# Identify all covariate columns we created (excluding outcome/ID cols)
id_cols <- c("cluster_id", "CountryName", "LATNUM", "LONGNUM",
             "stratum_id", "y", "n_children", "coverage")
all_cluster_covariates <- setdiff(names(cluster_data), id_cols)

cat("\nTotal cluster-level covariate columns:", length(all_cluster_covariates), "\n")
cat("Total clusters:", nrow(cluster_data), "\n")
cat("Countries:", n_distinct(cluster_data$CountryName), "\n\n")


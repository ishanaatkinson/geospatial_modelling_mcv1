# ════════════════════════════════════════════════════════════════════════════════
# SECTION 0: LIBRARIES
# ════════════════════════════════════════════════════════════════════════════════

library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(spdep)          # Moran's I, spatial weights
library(gstat)          # empirical variograms
library(car)            # VIF / GVIF
library(corrplot)       # correlation heatmaps
library(patchwork)      # combining ggplot panels
library(purrr)          # functional mapping
library(writexl)        # Excel output
library(pROC)           # AUC computation

# INLA (install if needed  not on CRAN, uses its own repository):
# install.packages("INLA",
#   repos = c(INLA = "https://inla.r-inla-download.org/R/stable"),
#   dep = TRUE)
inla_available <- requireNamespace("INLA", quietly = TRUE)
if (inla_available) library(INLA)

if (!inla_available) {
  cat("\n ï¸  R-INLA is not installed.\n")
  cat("To install, run:\n")
  cat('  install.packages("INLA",\n')
  cat('    repos = c(INLA = "https://inla.r-inla-download.org/R/stable"),\n')
  cat('    dep = TRUE)\n')
  cat("Models 2-4 (spatial) will be skipped. Model 1 (GLM) will still run.\n\n")
}

cat("
# ════════════════════════════════════════════════════════════════════════════════
═══                                                                  ═══
═══   MCV1 COVERAGE ANALYSIS PIPELINE                                ═══
═══   Following Utazi et al. (2022) methodology                     ═══
═══                                                                  ═══
# ════════════════════════════════════════════════════════════════════════════════
")


# ═══════════════════════════════════════════════════════════
# SECTION 1: LOAD DATA
# ═════════════════════════════════════════════════════════════════
#
# We load the final DHS data produced by the earlier pipeline scripts:
#   - Individual-level vaccination and covariate data
#   - GPS coordinates (LATNUM, LONGNUM)
#   - GC geospatial covariates (nightlights, elevation, aridity, etc.)
#   - ACLED conflict indicators
#   - Slum classification (Utazi method)
#   - Administrative region names (GADM spatial joins)
# ═════════════════════════════════════════════════════════════════

cat("\n----- SECTION 1: Loading data -----\n")

base_path <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT"

# final_DHS_data <- readRDS(
#   file.path(base_path, "Data/DHS_final/final_DHS_data.rds")
# )


# ═════════════════════════════════════════════════════════════════
# SECTION 2: DEFINE CANDIDATE COVARIATES
# ═════════════════════════════════════════════════════════════════
#
# These variables come from Table A in Utazi et al. (2022) Supplementary,
# augmented with geospatial covariates from Tessema et al. (2024) and
# Mosser et al. (2019).
#
# The key principle: assemble ALL available covariates as candidates, then
# let the DATA decide which enter each country's model through a systematic
# screening pipeline (missingness †’ bivariate OR †’ multicollinearity).
# Different countries may end up with different covariate sets because
# DHS questionnaires vary between countries.
#
# VARIABLE TYPES:
#   - CONTINUOUS variables (e.g. maternal age, wealth score) †’ at cluster
#     level these become the mean across children in that cluster
#   - CATEGORICAL/FACTOR variables (e.g. education level, employment status)
#     †’ at cluster level these become PROPORTIONS (e.g. "what fraction of
#     mothers in this cluster have primary education?")
# ═════════════════════════════════════════════════════════════════

cat("\n----- SECTION 2: Defining candidate covariates -----\n")

# ═══════════════════════════════════════════════════════════
# VARIABLE TYPE CLASSIFICATION
# ═══════════════════════════════════════════════════════════
#
# Each variable falls into one of three types, and this matters for how
# we aggregate from individual children to the cluster level:
#
#   BINARY (0/1):
#     These are yes/no variables like "did the mother have skilled birth
#     attendance?" (0 = No, 1 = Yes).
#     At the cluster level, the MEAN of a binary variable IS the PROPORTION.
#     So if 3 out of 5 mothers in a cluster had SBA, the cluster mean is
#     0.6 = 60%.  We name these columns cl_<varname>_prop to be clear.
#
#   CATEGORICAL (3+ ordered or unordered levels):
#     These are variables like education level (0=None, 1=Primary,
#     2=Secondary, 3=Higher) or birth order (0=1st-2nd, 1=3rd-5th, 2=>5th).
#     The mean of the numeric code is NOT meaningful (e.g. "mean education
#     = 1.3" is not interpretable).
#     Instead, we create ONE PROPORTION COLUMN PER LEVEL at the cluster level:
#       cl_respondent_edu_level_prop_0 = proportion with no education
#       cl_respondent_edu_level_prop_1 = proportion with primary education
#       cl_respondent_edu_level_prop_2 = proportion with secondary education
#       cl_respondent_edu_level_prop_3 = proportion with higher education
#     These proportions sum to 1 within each cluster.
#     Following Utazi et al. (2022), we DROP ONE reference level (the most
#     common category) to avoid perfect collinearity.  The dropped level's
#     proportion is implicit (1 minus the sum of the others).
#
#   CONTINUOUS (truly numeric):
#     These are variables like wealth_index_score, travel time in minutes,
#     nighttime lights intensity, temperature, etc.
#     At the cluster level we compute the MEAN, and also a STANDARDISED
#     z-score (mean=0, SD=1) so that odds ratios represent "per 1-SD
#     increase" for comparability across variables with different units.
#
# WHY THIS MATTERS:
#   If we just take the mean of education level codes (0, 1, 2, 3), a
#   cluster with 50% "None" and 50% "Higher" would get the same mean (1.5)
#   as a cluster with 100% "Primary"  obviously very different populations.
#   By computing separate proportions for each level, we preserve this
#   information and can ask: "does % with primary education in a cluster
#   predict coverage differently from % with secondary education?"
# ═══════════════════════════════════════════════════════════

if (run_type == "full") {

# ----- 2a. INDIVIDUAL-LEVEL (child characteristics) ----------------------------------------------------------------------
#
# From Table A "Individual level": sex, birth order, SBA, birth quarter,
# ANC visits, TT vaccination, PNC, health card, vitamin A
individual_vars <- c(
  "child_sex",                     # Male=0, Female=1
  "birth_order",                   # 1-2=0, 3-5=1, >5=2
  "skilled_birth_attendance",      # No=0, Yes=1
  "birth_quarter",                 # Jan-Mar=0, Apr-Jun=1, Jul-Sep=2, Oct-Dec=3
  "anc_visits_number",             # No ANC=0, 1-3=1, 4+=2
  "anc_provider",                  # Skilled=1, TBA=2, Relative=3, None=4
  "postnatal_care_newborn",        # No=0, Yes=1
  "postnatal_care_newborn_timing", # PNC within 2 months: No=0, Yes=1
  # "child_health_card",             # Has health/vaccination card: No=0, Yes=1
  "child_vitamin_a"                # Received vitamin A supplement: No=0, Yes=1
)




# ----- 2b. HOUSEHOLD-LEVEL (maternal & household characteristics) -----------------------------------”€
#
# From Table A "Household level": maternal age, education, employment,
# media exposure, wealth, healthcare barriers, household composition, etc.
household_vars <- c(
  # Maternal demographics
  "respondent_age",                # 15-19=0, 20-29=1, 30-39=2, 40-49=3
  "marital_status",                # Never=0, Married=1, Div/Wid/Sep=2
  "respondent_working_12mo",       # Worked in past 12mo: No=0, Yes=1
  "respondent_occupation",         # Agric=0, Clerical/Sales=1, Prof=2, Manual=3
  
  # Healthcare barriers (Table A: "Mother had problem seeking medical advice")
  # Polarity: Had problem=0, No problem=1
  "healthcare_barriers",           # ANY barrier composite
  # "healthcare_barriers_cost",      # Cost barrier
  # "healthcare_barriers_transport", # Transport/distance barrier
  # "healthcare_barriers_safety",    # Safety barrier
  # "healthcare_barriers_stockout",  # No provider/drugs barrier
  
  # Education
  "respondent_edu_level",          # None=0, Primary=1, Secondary=2, Higher=3
  "edu_grouped",
  # "partner_edu_level",             # None=0, Primary=1, Secondary=2, Higher=3
  # "respondent_literacy",           # Cannot read=0, Parts=1, Whole=2, No card=3
  
  # Media exposure
  "media_exposure",                # Any media weekly: No=0, Yes=1
  # "media_radio",                   # 0=never, 1=<weekly, 2=weekly, 3=daily
  # "media_television",              # 0=never, 1=<weekly, 2=weekly, 3=daily
  # "media_newspaper",               # 0=never, 1=<weekly, 2=weekly, 3=daily
  
  # Mobile phone / internet
  "mobile_internet_use",           # No=0, Yes=1
  # "hh_mobile_phone",               # Household owns mobile: No=0, Yes=1
  # "internet_use",                  # Uses internet: No=0, Yes=1
  
  # Health insurance & health knowledge
  "health_insurance",              # No=0, Yes=1
  #"hh_mosquito_net",               # Household owns mosquito net: No=0, Yes=1
  #"mother_knows_malaria",          # Mother knows about malaria: No=0, Yes=1
  
  # Wealth
  "wealth_index_quintile",         # Poor=0, Middle=1, Rich=2
  # "wealth_index_score",            # Continuous factor score (z-score)
  
  # Household composition
  "household_size",                # Large(>=9)=0, Medium(5-8)=1, Small(<=4)=2
  "hh_head_sex",                   # Female=0, Male=1
  "place_of_delivery",             # Home=1, Govt=2, Private=3
  "length_of_stay",                # <1yr=0, 1-3=1, 4-5=2, >5/always=3
  
  # Women's autonomy (from Acharya et al. 2018)
  "womens_autonomy_has_a_say",     # Count of decision-making domains (0-6)
  
  # Slum dwelling
  "hh_slum_dwelling"               # UN-Habitat criteria slum indicator
)

# ----- 2c. COMMUNITY/CLUSTER-LEVEL --------------------------------------------------------------------------------------------------------------”€
#
# From Table A "Community/cluster level": urban/rural, conflict, urban slum
community_vars <- c(
  "urban_rural",                   # Rural=0, Urban=1
  "conflict_area",                 # From ACLED merge
  "urban_slum"                     # From Utazi slum classification
)

# ═══════════════════════════════════════════════════════════
# EXPLICIT VARIABLE TYPE CLASSIFICATION
# ═══════════════════════════════════════════════════════════
#
# We explicitly classify every variable as binary, categorical, or continuous.
# This determines how it is aggregated from individual to cluster level.

# BINARY variables (0/1 only  cluster mean = proportion with "1")
# At cluster level: mean = proportion. E.g. cl_skilled_birth_attendance_prop
# = 0.60 means 60% of mothers in that cluster had SBA.
binary_vars <- c(
  "child_sex",                      # 0=Male, 1=Female †’ prop female
  "skilled_birth_attendance",       # 0=No, 1=Yes †’ prop with SBA
  "postnatal_care_newborn",         # 0=No, 1=Yes †’ prop with PNC
  "postnatal_care_newborn_timing",  # 0=No, 1=Yes †’ prop PNC within 2mo
  # "child_health_card",              # 0=No, 1=Yes †’ prop with health card
  # "child_vitamin_a",                # 0=No, 1=Yes †’ prop received vitamin A
  "respondent_working_12mo",        # 0=No, 1=Yes †’ prop mothers working
  "healthcare_barriers",            # 0=Has barrier, 1=No barrier †’ prop no barrier
  # "healthcare_barriers_cost",       # 0=Has cost barrier, 1=No †’ prop no cost barrier
  # "healthcare_barriers_transport",  # 0=Has transport barrier, 1=No
  # "healthcare_barriers_safety",     # 0=Has safety barrier, 1=No
  # "healthcare_barriers_stockout",   # 0=Has stockout barrier, 1=No
  "media_exposure",                 # 0=No, 1=Yes †’ prop with weekly media
  "mobile_internet_use",            # 0=No, 1=Yes †’ prop using mobile/internet
  # "hh_mobile_phone",                # 0=No, 1=Yes †’ prop with mobile phone
  # "internet_use",                   # 0=No, 1=Yes †’ prop using internet
  "health_insurance",               # 0=No, 1=Yes †’ prop insured
  #"hh_mosquito_net",                # 0=No, 1=Yes †’ prop with mosquito net
  #"mother_knows_malaria",           # 0=No, 1=Yes †’ prop with malaria knowledge
  "hh_head_sex",                    # 0=Female, 1=Male †’ prop male-headed HH
  "hh_slum_dwelling",               # 0=No, 1=Yes †’ prop in slum
  "urban_rural",                    # 0=Rural, 1=Urban †’ prop urban
  "conflict_area",                  # 0=Conflict, 1=Non-conflict †’ prop non-conflict
  "urban_slum",                     # 0=Slum, 1=Non-slum †’ prop non-slum
  "literacy_binary",                # 1=can read parts or whole sentence, 0=Can't read
  "anc_skilled",                    # 1=Skilled, 0=Unskilled
  "facility_delivery"               # 1=Home, 2=Private or Govt 
)

# CATEGORICAL variables (3+ discrete levels  need per-level proportions)
# At cluster level: one column per level giving the proportion of children
# in that level.  E.g. for respondent_edu_level with levels {0,1,2,3}:
#   cl_respondent_edu_level_prop_0 = proportion with no education
#   cl_respondent_edu_level_prop_1 = proportion with primary education
#   ... etc.
# The REFERENCE LEVEL (most common category, following Utazi) is dropped
# to avoid perfect collinearity (proportions sum to 1).
categorical_vars <- c(
  "birth_order",                    # 0=1st-2nd, 1=3rd-5th, 2=>5th
  "birth_quarter",                  # 0=Jan-Mar, 1=Apr-Jun, 2=Jul-Sep, 3=Oct-Dec
  "anc_visits_number",              # 0=No ANC, 1=1-3, 2=4+
  #"anc_provider",                   # 1=Skilled, 2=TBA, 3=Relative, 4=None
  "respondent_age",                 # 0=15-19, 1=20-29, 2=30-39, 3=40-49
  "marital_status",                 # 0=Never, 1=Married, 2=Div/Wid/Sep
  "respondent_occupation",          # 0=Agric, 1=Clerical/Sales, 2=Prof, 3=Manual
  #"respondent_edu_level",           # 0=None, 1=Primary, 2=Secondary, 3=Higher
  # "partner_edu_level",              # 0=None, 1=Primary, 2=Secondary, 3=Higher
  # "respondent_literacy",            # 0=Cannot, 1=Parts, 2=Whole, 3=No card
  # "media_radio",                    # 0=never, 1=<weekly, 2=weekly, 3=daily
  # "media_television",               # 0=never, 1=<weekly, 2=weekly, 3=daily
  # "media_newspaper",                # 0=never, 1=<weekly, 2=weekly, 3=daily
  "wealth_index_quintile",          # 0=Poor, 1=Middle, 2=Rich
  "household_size",                 # 0=Large(>=9), 1=Medium(5-8), 2=Small(<=4)
  #"place_of_delivery",              # 1=Home, 2=Govt, 3=Private
  "length_of_stay",                  # 0=<1yr, 1=1-3, 2=4-5, 3=>5/always
  "edu_grouped"
)

# CONTINUOUS variables (truly numeric  cluster mean + z-score)
# At cluster level: mean value across children.
# E.g. cl_wealth_index_score_mean = mean wealth z-score in that cluster.
continuous_vars <- c(
  # "wealth_index_score",             # DHS wealth factor score (continuous)
  "womens_autonomy_has_a_say"       # Count of decision-making domains (0-6)
)

# Geospatial covariates (all continuous  measured at cluster GPS point)
# These don't need aggregation since they're already cluster-level, but
# we include them in the continuous classification for z-scoring.
# gc_vars_base are handled separately via pick_gc_year().

# Travel time variables are continuous
# travel_time_vars are all continuous (minutes)

# Human-readable labels for categorical variable levels
# (used for interpretable column names and summary tables)
categorical_level_labels <- list(
  birth_order          = c("0" = "1st_2nd", "1" = "3rd_5th", "2" = "gt5th"),
  birth_quarter        = c("0" = "jan_mar", "1" = "apr_jun", "2" = "jul_sep", "3" = "oct_dec"),
  anc_visits_number    = c("0" = "no_anc", "1" = "1to3", "2" = "4plus"),
  #anc_provider         = c("1" = "skilled", "2" = "tba/relative/none"),
  respondent_age       = c("0" = "15_19", "1" = "20_29", "2" = "30_39", "3" = "40_49"),
  marital_status       = c("0" = "never_married", "1" = "married", "2" = "div_wid_sep"),
  respondent_occupation = c("0" = "agriculture", "1" = "clerical_sales", "2" = "professional", "3" = "manual"),
  #respondent_edu_level = c("0" = "no_education", "1" = "primary_or_above"),
  # partner_edu_level    = c("0" = "no_education", "1" = "primary", "2" = "secondary", "3" = "higher"),
  #respondent_literacy  = c("0" = "cannot_read", "1" = "parts_or_whole_sentence"),
  # media_radio          = c("0" = "never", "1" = "lt_weekly", "2" = "weekly", "3" = "daily"),
  # media_television     = c("0" = "never", "1" = "lt_weekly", "2" = "weekly", "3" = "daily"),
  # media_newspaper      = c("0" = "never", "1" = "lt_weekly", "2" = "weekly", "3" = "daily"),
  wealth_index_quintile = c("0" = "poor", "1" = "middle", "2" = "rich"),
  household_size       = c("0" = "large_9plus", "1" = "medium_5to8", "2" = "small_lto4"),
  place_of_delivery    = c("1" = "home", "2" = "govt_facility/private_facility"),
  length_of_stay       = c("0" = "lt_1yr", "1" = "1to3yr", "2" = "4to5yr", "3" = "gt5yr_always"),
  edu_grouped          = c("0" = "No_education", "1" = "Primary", "2" ="Seconary_or_higher")
)

cat("Variable type classification:\n")
cat("  BINARY (0/1):     ", length(binary_vars), "variables\n")
cat("  CATEGORICAL (3+): ", length(categorical_vars), "variables\n")
cat("  CONTINUOUS:        ", length(continuous_vars), "+ GC + travel time\n")

# ----- 2d. GEOSPATIAL COVARIATES (DHS GC dataset) ---------------------------------------------------------------------------
#
# Environmental/remoteness variables extracted at cluster GPS locations.
# These have year suffixes (_2000, _2005, _2010, _2015, _2020); we pick
# the year closest to each survey (see pick_gc_year function below).
gc_vars_base <- c(
  "NIGHTLIGHTS_COMPOSITE",         # Nighttime lights (proxy for development)
  "ELEVATION",                     # Elevation in metres (access barrier)
  # "ARIDITY",                       # Aridity index (climate stress)
  # "PRECIPITATION",                 # Annual precipitation
  # "MEAN_TEMPERATURE",              # Mean annual temperature
  # "ENHANCED_VEGETATION_INDEX",     # EVI (vegetation/agriculture)
  # "UN_POPULATION_DENSITY",         # Population density (WorldPop/UN)
  "TRAVEL_TIMES"                  # Travel time to nearest city (Weiss et al.)
  # "GLOBAL_HUMAN_FOOTPRINT",        # Human footprint index
  # "MALARIA_PREVALENCE",            # P. falciparum prevalence (MAP)
  # "DROUGHT_EPISODES"               # Drought episodes
)

# ----- 2e. TRAVEL TIME VARIABLES (from GPS merge) ---------------------------------------------------------------------------
#
# Utazi et al. (2022) Table A: travel time to nearest city ‰¥ 50,000 people,
# classified into country-specific tertiles.
travel_time_vars <- c(
  "travel_time_to_city",            # Minutes to nearest city >= 50k pop
  "travel_time_to_HC_motor",        # Minutes to nearest health centre (motor)
  "travel_time_to_HC_walk"          # Minutes to nearest health centre (walk)
)

# ----- Print summary ------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Candidate variable set:\n")
cat("  Individual-level:", length(individual_vars), "vars\n")
cat("  Household-level: ", length(household_vars), "vars\n")
cat("  Community-level: ", length(community_vars), "vars\n")
cat("  Geospatial (GC): ", length(gc_vars_base), "base vars\n")
cat("  Travel time:     ", length(travel_time_vars), "vars\n")
cat("  TOTAL candidates: ~",
    length(individual_vars) + length(household_vars) +
      length(community_vars) + length(gc_vars_base) +
      length(travel_time_vars), "\n")


# ═════════════════════════════════════════════════════════════════
# SECTION 3: PICK GC YEAR & PREPARE INDIVIDUAL DATA
# ═════════════════════════════════════════════════════════════════
#
# GC variables come in multiple time-stamped versions (2000, 2005, 2010, 2015,
# 2020).  For each child, we pick the version closest to their survey year.
# For example, a child surveyed in 2018 gets the 2020 nightlights value,
# while a child surveyed in 2012 gets the 2010 value.
# ═════════════════════════════════════════════════════════════════

cat("\n----- SECTION 3: Picking GC year & preparing data -----\n")

# Combine all surveys into a single data frame
# all_surveys <- bind_rows(final_DHS_data, .id = "survey_idx")
# rm(final_DHS_data)  # free memory

# Helper: for a given GC base name (e.g. "NIGHTLIGHTS_COMPOSITE"),
# find all year-suffixed columns, and for each row pick the one whose
# year is closest to that row's SurveyYear.
pick_gc_year <- function(df, base_name) {
  gc_cols <- grep(paste0("^", base_name, "_\\d{4}$"), names(df), value = TRUE)
  if (length(gc_cols) == 0) {
    if (base_name %in% names(df)) return(as.numeric(df[[base_name]]))
    return(rep(NA_real_, nrow(df)))
  }
  gc_years <- as.numeric(gsub(paste0(base_name, "_"), "", gc_cols))
  sy <- as.numeric(df$SurveyYear)
  dist_mat <- abs(outer(sy, gc_years, `-`))
  best_idx <- max.col(-dist_mat, ties.method = "first")
  val_mat  <- as.matrix(df[, gc_cols, drop = FALSE])
  storage.mode(val_mat) <- "double"
  val_mat[cbind(seq_len(nrow(df)), best_idx)]
}


# ═════════════════════════════════════════════════════════════════
# SECTION 4: FILTER FOR ELIGIBILITY
# ═════════════════════════════════════════════════════════════════
#
# MCV1 is recommended at age 9€“12 months.  Following standard DHS practice,
# we analyse children aged 12€“23 months  old enough to have had the
# opportunity to be vaccinated, but not so old that recall bias dominates.
#
# We also require:
#   - Non-missing vaccination data (child_vacc_measles)
#   - Valid GPS coordinates (LATNUM/LONGNUM, excluding 0/0 which is a
#     DHS missing value code for "coordinates suppressed")
#   - A valid cluster PSU identifier
# ═════════════════════════════════════════════════════════════════
# 
# cat("\n----- SECTION 4: Filtering for eligible children -----\n")
# 
# eligible <- all_surveys %>%
#   filter(
#     !is.na(child_age_months), child_age_months >= 12, child_age_months <= 23,
#     !is.na(child_vacc_measles),
#     !is.na(cluster_psu),
#     !is.na(LATNUM), LATNUM != 0,
#     !is.na(LONGNUM), LONGNUM != 0
#   ) %>%
#   mutate(
#     # Binary MCV1 indicator: vaccinated from any source
#     # (card dated = 1, mother's report = 2, card marked = 3 in DHS coding)
#     mcv1_received = as.integer(child_vacc_measles %in% c(1, 2, 3))
#   )
# 
# cat("Eligible children (12-23 months with MCV1 data & GPS):", nrow(eligible), "\n")
# cat("Countries:", n_distinct(eligible$CountryName), "\n")
# 
# # Save intermediate file for reproducibility
# data.table::fwrite(eligible, paste0(base_path, "/Data/DHS_final/model_data_v2.csv"))
model_data <- vroom::vroom(paste0(base_path, "/Data/DHS_final/model_data_v2.csv"))




#rm(eligible, all_surveys)  # free memory

# Add year-matched GC columns
gc_col_names <- paste0("gc_", tolower(gc_vars_base))
gc_results   <- lapply(gc_vars_base, function(gc) pick_gc_year(model_data, gc))
names(gc_results) <- gc_col_names
model_data <- cbind(model_data, as.data.frame(gc_results))
gc_available <- gc_col_names[gc_col_names %in% names(model_data)]

# Save with GC columns added
# data.table::fwrite(model_data, paste0(base_path, "/Data/DHS_final/model_data_v3.csv"))
# model_data <- vroom::vroom(paste0(base_path, "/Data/DHS_final/model_data_v3.csv"))
gc_col_names <- paste0("gc_", tolower(gc_vars_base))
gc_available <- gc_col_names[gc_col_names %in% names(model_data)]

# ----- Before Section 5 aggregation, recode categorical variables -----------------------------------
# 
# We collapse levels that are substantively similar into groups.
# This is done on the INDIVIDUAL-level data so that when we compute
# cluster-level proportions, each group is a single clean column.

model_data <- model_data %>%
  filter(SurveyYear >=2000) %>% 
  mutate(
    # Literacy: binary (can read vs cannot)
    literacy_binary = case_when(
      respondent_literacy %in% c(1, 2) ~ 1L,  # can read (parts or whole)
      respondent_literacy %in% c(0, 3) ~ 0L,  # cannot read / no card
      TRUE ~ NA_integer_
    ),
    
    # Education: collapse to 3 groups  
    edu_grouped = case_when(
      respondent_edu_level == 0 ~ 0L,          # no education
      respondent_edu_level == 1 ~ 1L,          # primary
      respondent_edu_level %in% c(2, 3) ~ 2L,  # primary, secondary or higher
      TRUE ~ NA_integer_
    ),
    
    # ANC provider: binary (skilled vs not)
    anc_skilled = case_when(
      anc_provider == 1 ~ 1L,                  # skilled
      anc_provider %in% c(2, 3, 4) ~ 0L,       # TBA/relative/none
      TRUE ~ NA_integer_
    ),
    
    # Place of delivery: binary (facility vs home)
    facility_delivery = case_when(
      place_of_delivery == 1 ~ 1L,                  # home
      place_of_delivery %in% c(2, 3, 4) ~ 0L,       # govt or private
      TRUE ~ NA_integer_
    )
  )




# ═════════════════════════════════════════════════════════════════
# SECTION 5: AGGREGATE TO CLUSTER LEVEL
# ═════════════════════════════════════════════════════════════════
#
# THE KEY IDEA: individual-level DHS data must be aggregated to the cluster
# level for the geostatistical model.  But HOW we aggregate depends on what
# TYPE of variable it is:
#
#   BINARY (0/1) variables like "had skilled birth attendance?"
#     †’ Cluster MEAN of 0/1 = PROPORTION.
#     Example: cl_skilled_birth_attendance = 0.60 means "60% of mothers in
#     this cluster had SBA".  This is already interpretable and correct.
#     Column naming: cl_<varname>   (the mean IS the proportion)
#
#   CATEGORICAL (3+ levels) variables like "education level" (0,1,2,3)
#     †’ Cluster mean of the code is MEANINGLESS.
#     A mean of 1.5 could mean "everyone primary-to-secondary" OR
#     "half no education, half higher education"  totally different.
#     INSTEAD: create separate proportion columns for each level.
#     Example: for respondent_edu_level with levels {0,1,2,3}:
#       cl_edu_no_education = 0.40  (40% of mothers have no education)
#       cl_edu_primary      = 0.30  (30% have primary)
#       cl_edu_secondary    = 0.20  (20% have secondary)
#       cl_edu_higher       = 0.10  (10% have higher)
#     These sum to 1 within each cluster, so we DROP ONE (the reference
#     level  the most common category) to avoid perfect collinearity.
#     Following Utazi et al. (2022): "Wherever possible, we considered
#     the reference category to be the category that had the least
#     likelihood of being vaccinated based on prior knowledge."
#
#   CONTINUOUS (truly numeric) variables like wealth_index_score, travel time
#     †’ Cluster MEAN is appropriate.
#     We also compute a z-score (standardised) version so that odds ratios
#     represent "per 1-SD increase" for comparability.
#     Column naming: cl_<varname> (mean) and cl_<varname>_z (z-score)
#
# WHY THIS MATTERS:
#   If we just average education codes (0,1,2,3), a cluster with 50% "None"
#   and 50% "Higher" gets the same mean (1.5) as a cluster with 100%
#   "Primary".  These are obviously very different populations.  By creating
#   separate proportion columns for each level, we preserve this information
#   and can ask: "does a 10% increase in % with primary education predict
#   coverage differently from a 10% increase in % with higher education?"
# ═════════════════════════════════════════════════════════════════

cat("\n----- SECTION 5: Aggregating to cluster level -----\n")

# ----- 5a. Identify which variables are present in the data --------------------------------------------------”€
all_cov <- c(individual_vars, household_vars, community_vars, travel_time_vars)
all_cov_in_data <- intersect(c(all_cov, gc_available), names(model_data))
cat("Covariates found in data:", length(all_cov_in_data),
    "of", length(all_cov) + length(gc_available), "candidates\n")

# Classify the variables that are actually present in the data
binary_in_data     <- intersect(binary_vars, names(model_data))
categorical_in_data <- intersect(categorical_vars, names(model_data))
continuous_in_data  <- intersect(c(continuous_vars, travel_time_vars), names(model_data))
gc_in_data          <- gc_available  # already cluster-level, just need z-scores

cat("  Binary in data:     ", length(binary_in_data), "\n")
cat("  Categorical in data:", length(categorical_in_data), "\n")
cat("  Continuous in data: ", length(continuous_in_data), "\n")
cat("  GC in data:         ", length(gc_in_data), "\n")


# ----- 5b. Aggregate: outcome + binary variables + continuous variables --------------------”€
#
# Binary: mean of 0/1 = proportion (this is correct and interpretable)
# Continuous: mean across children in the cluster
cat("\n----- 5b: Aggregating binary and continuous covariates -----\n")

# Combine binary + continuous for the simple mean aggregation
simple_mean_vars <- intersect(c(binary_in_data, continuous_in_data, gc_in_data),
                              names(model_data))


cluster_data <- model_data %>%
  group_by(CountryName, country_phase, cluster_psu, SurveyYear, LATNUM, LONGNUM) %>%
  summarise(
    # Outcome
    n_children   = n(),
    n_vaccinated = sum(mcv1, na.rm = TRUE),
    coverage     = n_vaccinated / n_children,
    
    # Binary variables †’ proportion (mean of 0/1)
    # Continuous variables †’ cluster mean
    across(all_of(simple_mean_vars), ~ mean(.x, na.rm = TRUE), .names = "cl_{.col}"),
    
    .groups = "drop"
  ) %>%
  # Empirical logit (Haldane-Anscombe correction):
  # Handles logit(0) = -Inf and logit(1) = +Inf by adding 0.5
  mutate(
    coverage_emp_logit = log((n_vaccinated + 0.5) / (n_children - n_vaccinated + 0.5))
  )

cat("Clusters after binary/continuous aggregation:", nrow(cluster_data), "\n")


# ----- 5c. Aggregate: categorical variables †’ per-level proportions ------------------------------”€
#
# For each categorical variable (e.g. respondent_edu_level with codes 0,1,2,3),
# we compute the PROPORTION of children in each level within each cluster.
# This creates columns like:
#   cl_edu_primary       = fraction of mothers with primary education
#   cl_edu_secondary     = fraction of mothers with secondary education
#   cl_edu_higher        = fraction of mothers with higher education
#   (reference level "no_education" is dropped  its proportion = 1 minus
#    the sum of the others)
#
# WHY DROP THE REFERENCE LEVEL?
# Because the proportions sum to 1 within each cluster, including ALL levels
# would create perfect multicollinearity (the last column is exactly
# determined by the others).  This is the same reason you drop one dummy
# variable in standard regression.  Following Utazi et al. (2022), the
# reference level is typically the most common category.

cat("\n----- 5c: Aggregating categorical covariates †’ per-level proportions -----\n")

prop_cols_created <- list()  # track which columns we created for each variable

for (cv in categorical_in_data) {
  x <- model_data[[cv]]
  if (all(is.na(x))) {
    cat("  ", cv, ": all NA  skipping\n")
    next
  }
  
  # Find the levels present in the data
  levs <- sort(unique(x[!is.na(x)]))
  if (length(levs) < 2) {
    cat("  ", cv, ": only 1 level  skipping\n")
    next
  }
  
  # Identify the reference level (most common = drop)
  #tab <- sort(table(x, useNA = "no"), decreasing = TRUE)
  ref_level <- min(levs)
  
  # Get human-readable labels for this variable
  labels <- categorical_level_labels[[cv]]
  
  # Create proportion columns for each NON-REFERENCE level
  non_ref_levs <- setdiff(as.character(levs), ref_level)
  created_cols <- c()
  
  for (lv in non_ref_levs) {
    # Create a readable column name
    if (!is.null(labels) && lv %in% names(labels)) {
      col_name <- paste0("cl_", cv, "_", labels[lv])
    } else {
      col_name <- paste0("cl_", cv, "_lev", lv)
    }
    # Make syntactically valid R name
    col_name <- make.names(col_name)
    
    # Compute per-cluster proportion for this level
    lv_props <- model_data %>%
      group_by(cluster_psu) %>%
      summarise(!!col_name := mean(.data[[cv]] == as.numeric(lv), na.rm = TRUE),
                .groups = "drop")
    
    cluster_data <- left_join(cluster_data, lv_props, by = "cluster_psu")
    created_cols <- c(created_cols, col_name)
  }
  
  # Also create the reference level proportion for descriptive purposes
  # (not used in the model, but useful for summary tables)
  if (!is.null(labels) && ref_level %in% names(labels)) {
    ref_col_name <- paste0("cl_", cv, "_", labels[1], "_REF")
  } else {
    ref_col_name <- paste0("cl_", cv, "_lev", ref_level, "_REF")
  }
  ref_col_name <- make.names(ref_col_name)
  ref_props <- model_data %>%
    group_by(cluster_psu) %>%
    summarise(!!ref_col_name := mean(.data[[cv]] == as.numeric(ref_level), na.rm = TRUE),
              .groups = "drop")
  cluster_data <- left_join(cluster_data, ref_props, by = "cluster_psu")
  
  prop_cols_created[[cv]] <- list(
    ref_level      = ref_level,
    ref_label      = if (!is.null(labels) && ref_level %in% names(labels)) labels[ref_level] else ref_level,
    ref_col        = ref_col_name,
    non_ref_cols   = created_cols,
    all_levels     = as.character(levs)
  )
  
  cat("  ", cv, ": ref='", ref_col_name, "' (dropped) | ",
      length(created_cols), " proportion cols created\n", sep = "")
}

cat("\nTotal categorical proportion columns created:",
    sum(sapply(prop_cols_created, function(x) length(x$non_ref_cols))), "\n")


# ----- 5d. Remove the old misleading mean-of-code columns for categorical vars ”€
#
# Now that we have proper per-level proportions, REMOVE the old cl_<varname>
# columns for categorical variables  these were the misleading mean-of-codes.
# For example, cl_respondent_edu_level (mean of 0,1,2,3 codes) is now
# replaced by cl_respondent_edu_level_primary, cl_..._secondary, etc.

cat("\n----- 5d: Removing misleading mean-of-code columns for categorical vars -----\n")

old_cat_cols <- paste0("cl_", categorical_in_data)
old_cat_present <- intersect(old_cat_cols, names(cluster_data))

if (length(old_cat_present) > 0) {
  cluster_data <- cluster_data %>% dplyr::select(-all_of(old_cat_present))
  cat("  Removed", length(old_cat_present), "mean-of-code columns:",
      paste(head(old_cat_present, 5), collapse = ", "),
      if (length(old_cat_present) > 5) "..." else "", "\n")
}


# ----- 5e. Standardise continuous covariates (z-scores) ------------------------------------------------------------”€
#
# Z-scores put all continuous variables on the same scale (mean=0, SD=1).
# This is useful because:
#   - Odds ratios become "per 1-SD increase" which is more interpretable
#     than "per 1-unit increase" when units differ (e.g. years vs km)
#   - It helps with model convergence in INLA
#   - It allows direct comparison of effect sizes across covariates
#
# We standardise ACROSS ALL COUNTRIES here (global z-scores).
# NOTE: we only z-score the truly continuous variables and GC variables,
# NOT the binary proportions or categorical proportions (which are already
# on a 0-1 scale and directly interpretable as "per 10% increase").

cat("\n----- 5e: Standardising continuous covariates (z-scores) -----\n")

# Identify continuous cl_ columns (NOT binary proportions, NOT categorical proportions)
continuous_cl_cols <- c(
  paste0("cl_", intersect(continuous_in_data, names(model_data))),
  paste0("cl_", gc_in_data)  # gc_ columns already have cl_ prefix from aggregation
)
# Fix: gc columns are named cl_gc_... from the across() call
continuous_cl_cols <- intersect(
  grep("^cl_", names(cluster_data), value = TRUE),
  c(paste0("cl_", continuous_in_data),
    grep("^cl_gc_", names(cluster_data), value = TRUE),
    grep("^cl_travel_", names(cluster_data), value = TRUE))
)

for (v in continuous_cl_cols) {
  v_z <- paste0(v, "_z")
  x <- cluster_data[[v]]
  s <- sd(x, na.rm = TRUE)
  if (!is.na(s) && s > 1e-10) {
    cluster_data[[v_z]] <- (x - mean(x, na.rm = TRUE)) / s
  }
}

cat("Added", sum(grepl("_z$", names(cluster_data))), "z-score columns for continuous vars\n")


# ----- 5f. Final cluster-level dataset summary -------------------------------------------------------------------------------------

cat("\n----- 5f: Cluster-level dataset summary -----\n")

cl_vars <- grep("^cl_", names(cluster_data), value = TRUE)
cl_vars_model <- cl_vars[!grepl("_REF$", cl_vars)]  # exclude reference-level cols

cat("Cluster-level dataset:", nrow(cluster_data), "clusters across",
    n_distinct(cluster_data$CountryName), "countries\n")
cat("Total cl_ columns:", length(cl_vars), "\n")
cat("  Model covariates (excl. reference):", length(cl_vars_model), "\n")
cat("  Reference-level columns (descriptive only):", sum(grepl("_REF$", cl_vars)), "\n")

# Print a summary of column types
cat("\n  Binary proportion cols (mean of 0/1):\n")
bin_cls <- paste0("cl_", binary_in_data)
bin_cls_present <- intersect(bin_cls, names(cluster_data))
cat("    ", paste(head(bin_cls_present, 6), collapse = ", "),
    if (length(bin_cls_present) > 6) paste0(" ... (", length(bin_cls_present), " total)") else "", "\n")

cat("  Categorical proportion cols (per-level %):\n")
cat_cls <- unlist(lapply(prop_cols_created, `[[`, "non_ref_cols"))
cat("    ", paste(head(cat_cls, 6), collapse = ", "),
    if (length(cat_cls) > 6) paste0(" ... (", length(cat_cls), " total)") else "", "\n")

cat("  Continuous cols (cluster mean):\n")
cat("    ", paste(head(continuous_cl_cols, 6), collapse = ", "),
    if (length(continuous_cl_cols) > 6) paste0(" ... (", length(continuous_cl_cols), " total)") else "", "\n")


# ----- 5g. Save cluster-level dataset ---------------------------------------------------------------------------------------------------------
# data.table::fwrite(cluster_data, paste0(base_path, "/Data/DHS_final/cluster_data_v2.csv"))
# cluster_data <- vroom::vroom(paste0(base_path, "/Data/DHS_final/cluster_data_v2.csv"))


# ═════════════════════════════════════════════════════════════════
# SECTION 6: DESCRIPTIVE SUMMARIES OF CLUSTER-LEVEL COVARIATE DISTRIBUTIONS
# ═════════════════════════════════════════════════════════════════
#
# Now that we have proper cluster-level proportions, we print descriptive
# summaries showing what a "typical cluster" looks like in the data.
# This helps contextualise the model results: if 80% of clusters have
# <10% with higher education, the model can't tell us much about the
# effect of higher education.
#
# These summaries are at the CLUSTER LEVEL (not individual level), because
# the cluster is the unit of analysis in the model.
# ═════════════════════════════════════════════════════════════════

cat("\n----- SECTION 6: Cluster-level covariate distribution summaries -----\n")

# Print summary of binary proportions across clusters
cat("\n----- Binary covariates (mean proportion across all clusters) -----\n")
cat("  (Read as: 'on average, X% of individuals within a cluster have this trait')\n\n")

for (bv in bin_cls_present) {
  x <- cluster_data[[bv]]
  if (!all(is.na(x))) {
    cat(sprintf("  %-45s mean=%.1f%%  median=%.1f%%  SD=%.1f%%\n",
                bv,
                mean(x, na.rm = TRUE) * 100,
                median(x, na.rm = TRUE) * 100,
                sd(x, na.rm = TRUE) * 100))
  }
}

# Print summary of categorical proportions across clusters
cat("\n----- Categorical covariates (mean proportion across all clusters) -----\n")
cat("  (Read as: 'on average, X% of individuals within a cluster are in this category')\n\n")

for (cv in names(prop_cols_created)) {
  info <- prop_cols_created[[cv]]
  cat("  ", cv, " (reference = '", info$ref_label, "', dropped):\n", sep = "")
  
  # Print reference level mean
  if (info$ref_col %in% names(cluster_data)) {
    ref_x <- cluster_data[[info$ref_col]]
    if (!all(is.na(ref_x)))
      cat(sprintf("    %-42s mean=%.1f%% [REF  not in model]\n",
                  info$ref_col, mean(ref_x, na.rm = TRUE) * 100))
  }
  # Print non-reference levels
  for (nc in info$non_ref_cols) {
    if (nc %in% names(cluster_data)) {
      x <- cluster_data[[nc]]
      if (!all(is.na(x)))
        cat(sprintf("    %-42s mean=%.1f%%  SD=%.1f%%\n",
                    nc, mean(x, na.rm = TRUE) * 100, sd(x, na.rm = TRUE) * 100))
    }
  }
}

# Print summary of continuous variables across clusters
cat("\n----- Continuous covariates (distribution across clusters) -----\n\n")
for (cv in continuous_cl_cols) {
  if (cv %in% names(cluster_data)) {
    x <- cluster_data[[cv]]
    if (!all(is.na(x)))
      cat(sprintf("  %-45s mean=%.2f  SD=%.2f  range=[%.2f, %.2f]\n",
                  cv, mean(x, na.rm = TRUE), sd(x, na.rm = TRUE),
                  min(x, na.rm = TRUE), max(x, na.rm = TRUE)))
  }
}


# ═════════════════════════════════════════════════════════════════
# SECTION 7: COVERAGE-BY-CLUSTER-SIZE TABLE (0/1 BOUNDARY DIAGNOSTICS)
# ═════════════════════════════════════════════════════════════════
#
# This table shows why so many clusters have 0% or 100% coverage.
# INTERPRETATION: clusters with few children (1-3) often show extreme
# coverage simply because the denominator is too small.  For example,
# if a cluster has only 2 children, coverage can only be 0%, 50%, or 100%.
# This is a SAMPLING ARTEFACT, not a real signal.  The binomial likelihood
# handles this naturally because 0/2 and 2/2 are perfectly valid binomial
# outcomes  the model knows that a cluster with 2/2 vaccinated is NOT the
# same as a cluster with 20/20.
# ═════════════════════════════════════════════════════════════════

cat("\n----- SECTION 7: Coverage by cluster size (0/1 boundary diagnostic) -----\n")

cat("\nOverall coverage distribution:\n")
cat("  Exactly 0%:       ", round(mean(cluster_data$coverage == 0) * 100, 1), "%\n")
cat("  Exactly 100%:     ", round(mean(cluster_data$coverage == 1) * 100, 1), "%\n")
cat("  Between 0% & 100%:", round(mean(cluster_data$coverage > 0 &
                                         cluster_data$coverage < 1) * 100, 1), "%\n")

cov_by_size <- cluster_data %>%
  mutate(sz = cut(n_children, c(0, 3, 5, 10, 20, Inf),
                  labels = c("1-3", "4-5", "6-10", "11-20", "21+"))) %>%
  group_by(sz) %>%
  summarise(
    n_clusters    = n(),
    pct_zero      = round(mean(coverage == 0) * 100, 1),
    pct_one       = round(mean(coverage == 1) * 100, 1),
    pct_extreme   = round(mean(coverage == 0 | coverage == 1) * 100, 1),
    mean_coverage = round(mean(coverage) * 100, 1),
    .groups       = "drop"
  )

cat("\nCoverage extremes by cluster size:\n")
print(cov_by_size, n = Inf)

cat("\nINTERPRETATION:\n")
cat("  Small clusters (1-3 children) will often show 0% or 100%\n")
cat("  simply because the denominator is too small for intermediate values.\n")
cat("  The binomial likelihood handles this: 0/3 and 3/3 are valid.\n")
cat("  The beta-binomial (Model 4) tests if overdispersion remains.\n")


# ═════════════════════════════════════════════════════════════════
# SECTION 8: COUNTRY-BY-COUNTRY MODELLING LOOP
# ═════════════════════════════════════════════════════════════════
#
# For each country with ‰¥ 50 geo-located clusters, we run the full
# Utazi et al. (2022) pipeline:
#
#   Step 1: Missingness screening (<5% rule)
#   Step 2: Bivariate screening (crude ORs, keep p < 0.2)
#   Step 3: Multicollinearity check (|r| > 0.8 then iterative GVIF)
#   Step 4: Model 1  covariates-only binomial GLM
#   Step 5: Model 2  spatial-only INLA-SPDE
#   Step 6: Model 3  full INLA-SPDE (covariates + spatial + nugget)
#   Step 7: Model 4  beta-binomial INLA-SPDE
#   Step 8: Model comparison + diagnostic plots + Excel output
# ═════════════════════════════════════════════════════════════════

cat("\n----- SECTION 8: Beginning country-by-country modelling loop -----\n")

# Identify countries with enough data
countries <- cluster_data %>%
  count(CountryName) %>%
  filter(n >= 50) %>%
  arrange(desc(n)) %>%
  pull(CountryName)

cat("\n", length(countries), "countries with >= 50 geo-located clusters:\n")
cat(paste("  ", countries), sep = "\n")

} else {
  
  
  # ----- 2a. INDIVIDUAL-LEVEL (child characteristics) ----------------------------------------------------------------------
  #
  # From Table A "Individual level": sex, birth order, SBA, birth quarter,
  # ANC visits, TT vaccination, PNC, health card, vitamin A
  individual_vars <- c(
    "skilled_birth_attendance"      # No=0, Yes=1
  )
  
  
  
  
  # ----- 2b. HOUSEHOLD-LEVEL (maternal & household characteristics) -----------------------------------”€
  #
  # From Table A "Household level": maternal age, education, employment,
  # media exposure, wealth, healthcare barriers, household composition, etc.
  household_vars <- c(
    "respondent_edu_level",          # None=0, Primary=1, Secondary=2, Higher=3
    "edu_grouped",
    "wealth_index_quintile"         # Poor=0, Middle=1, Rich=2

   
  )
  
  # ----- 2c. COMMUNITY/CLUSTER-LEVEL --------------------------------------------------------------------------------------------------------------”€
  #
  # From Table A "Community/cluster level": urban/rural, conflict, urban slum
  community_vars <- c(
    "urban_rural"                   # Rural=0, Urban=1

  )
  
  # ═══════════════════════════════════════════════════════════
  # EXPLICIT VARIABLE TYPE CLASSIFICATION
  # ═══════════════════════════════════════════════════════════
  #
  # We explicitly classify every variable as binary, categorical, or continuous.
  # This determines how it is aggregated from individual to cluster level.
  
  # BINARY variables (0/1 only  cluster mean = proportion with "1")
  # At cluster level: mean = proportion. E.g. cl_skilled_birth_attendance_prop
  # = 0.60 means 60% of mothers in that cluster had SBA.
  binary_vars <- c(
    "skilled_birth_attendance",       # 0=No, 1=Yes †’ prop with SBA
    "urban_rural"                    # 0=Rural, 1=Urban †’ prop urban
  )
  
  # CATEGORICAL variables (3+ discrete levels  need per-level proportions)
  # At cluster level: one column per level giving the proportion of children
  # in that level.  E.g. for respondent_edu_level with levels {0,1,2,3}:
  #   cl_respondent_edu_level_prop_0 = proportion with no education
  #   cl_respondent_edu_level_prop_1 = proportion with primary education
  #   ... etc.
  # The REFERENCE LEVEL (most common category, following Utazi) is dropped
  # to avoid perfect collinearity (proportions sum to 1).
  categorical_vars <- c(
    "wealth_index_quintile",          # 0=Poor, 1=Middle, 2=Rich
    "edu_grouped"
  )
  
  # CONTINUOUS variables (truly numeric  cluster mean + z-score)
  # At cluster level: mean value across children.
  # E.g. cl_wealth_index_score_mean = mean wealth z-score in that cluster.
  continuous_vars <- c()
    # "wealth_index_score",             # DHS wealth factor score (continuous)
  
  
  # Geospatial covariates (all continuous  measured at cluster GPS point)
  # These don't need aggregation since they're already cluster-level, but
  # we include them in the continuous classification for z-scoring.
  # gc_vars_base are handled separately via pick_gc_year().
  
  # Travel time variables are continuous
  # travel_time_vars are all continuous (minutes)
  
  # Human-readable labels for categorical variable levels
  # (used for interpretable column names and summary tables)
  categorical_level_labels <- list(
    wealth_index_quintile = c("0" = "poor", "1" = "middle", "2" = "rich"),
    edu_grouped          = c("0" = "No_education", "1" = "Primary", "2" ="Seconary_or_higher")
  )
  
  cat("Variable type classification:\n")
  cat("  BINARY (0/1):     ", length(binary_vars), "variables\n")
  cat("  CATEGORICAL (3+): ", length(categorical_vars), "variables\n")
  cat("  CONTINUOUS:        ", length(continuous_vars), "+ GC + travel time\n")
  
  # ----- 2d. GEOSPATIAL COVARIATES (DHS GC dataset) ---------------------------------------------------------------------------
  #
  # Environmental/remoteness variables extracted at cluster GPS locations.
  # These have year suffixes (_2000, _2005, _2010, _2015, _2020); we pick
  # the year closest to each survey (see pick_gc_year function below).
  gc_vars_base <- c()
  
  # ----- 2e. TRAVEL TIME VARIABLES (from GPS merge) ---------------------------------------------------------------------------
  #
  # Utazi et al. (2022) Table A: travel time to nearest city ‰¥ 50,000 people,
  # classified into country-specific tertiles.
  travel_time_vars <- c(
    "travel_time_to_HC_walk"          # Minutes to nearest health centre (walk)
  )
  
  # ----- Print summary ------------------------------------------------------------------------------------------------------------------------------------------------------
  cat("Candidate variable set:\n")
  cat("  Individual-level:", length(individual_vars), "vars\n")
  cat("  Household-level: ", length(household_vars), "vars\n")
  cat("  Community-level: ", length(community_vars), "vars\n")
  cat("  Geospatial (GC): ", length(gc_vars_base), "base vars\n")
  cat("  Travel time:     ", length(travel_time_vars), "vars\n")
  cat("  TOTAL candidates: ~",
      length(individual_vars) + length(household_vars) +
        length(community_vars) + length(gc_vars_base) +
        length(travel_time_vars), "\n")
  
  
  # ═════════════════════════════════════════════════════════════════
  # SECTION 3: PICK GC YEAR & PREPARE INDIVIDUAL DATA
  # ═════════════════════════════════════════════════════════════════
  #
  # GC variables come in multiple time-stamped versions (2000, 2005, 2010, 2015,
  # 2020).  For each child, we pick the version closest to their survey year.
  # For example, a child surveyed in 2018 gets the 2020 nightlights value,
  # while a child surveyed in 2012 gets the 2010 value.
  # ═════════════════════════════════════════════════════════════════
  
  model_data <- vroom::vroom(paste0(base_path, "/Data/DHS_final/model_data_v2.csv"))
  
  
  # # Add year-matched GC columns
  # gc_col_names <- paste0("gc_", tolower(gc_vars_base))
  # gc_results   <- lapply(gc_vars_base, function(gc) pick_gc_year(model_data, gc))
  # names(gc_results) <- gc_col_names
  # model_data <- cbind(model_data, as.data.frame(gc_results))
  # gc_available <- gc_col_names[gc_col_names %in% names(model_data)]
  # 
  # # Save with GC columns added
  # # data.table::fwrite(model_data, paste0(base_path, "/Data/DHS_final/model_data_v3.csv"))
  # # model_data <- vroom::vroom(paste0(base_path, "/Data/DHS_final/model_data_v3.csv"))
  # gc_col_names <- paste0("gc_", tolower(gc_vars_base))
  # gc_available <- gc_col_names[gc_col_names %in% names(model_data)]
  
  # ----- Before Section 5 aggregation, recode categorical variables -----------------------------------
  # 
  # We collapse levels that are substantively similar into groups.
  # This is done on the INDIVIDUAL-level data so that when we compute
  # cluster-level proportions, each group is a single clean column.
  
  model_data <- model_data %>%
    filter(SurveyYear >=2000) %>% 
    mutate(
      
      # Education: collapse to 3 groups  
      edu_grouped = case_when(
        respondent_edu_level == 0 ~ 0L,          # no education
        respondent_edu_level == 1 ~ 1L,          # primary
        respondent_edu_level %in% c(2, 3) ~ 2L,  # primary, secondary or higher
        TRUE ~ NA_integer_
      )
    )
  
  
  
  
  # ═════════════════════════════════════════════════════════════════
  # SECTION 5: AGGREGATE TO CLUSTER LEVEL
  # ═════════════════════════════════════════════════════════════════
  #
  # THE KEY IDEA: individual-level DHS data must be aggregated to the cluster
  # level for the geostatistical model.  But HOW we aggregate depends on what
  # TYPE of variable it is:
  #
  #   BINARY (0/1) variables like "had skilled birth attendance?"
  #     †’ Cluster MEAN of 0/1 = PROPORTION.
  #     Example: cl_skilled_birth_attendance = 0.60 means "60% of mothers in
  #     this cluster had SBA".  This is already interpretable and correct.
  #     Column naming: cl_<varname>   (the mean IS the proportion)
  #
  #   CATEGORICAL (3+ levels) variables like "education level" (0,1,2,3)
  #     †’ Cluster mean of the code is MEANINGLESS.
  #     A mean of 1.5 could mean "everyone primary-to-secondary" OR
  #     "half no education, half higher education"  totally different.
  #     INSTEAD: create separate proportion columns for each level.
  #     Example: for respondent_edu_level with levels {0,1,2,3}:
  #       cl_edu_no_education = 0.40  (40% of mothers have no education)
  #       cl_edu_primary      = 0.30  (30% have primary)
  #       cl_edu_secondary    = 0.20  (20% have secondary)
  #       cl_edu_higher       = 0.10  (10% have higher)
  #     These sum to 1 within each cluster, so we DROP ONE (the reference
  #     level  the most common category) to avoid perfect collinearity.
  #     Following Utazi et al. (2022): "Wherever possible, we considered
  #     the reference category to be the category that had the least
  #     likelihood of being vaccinated based on prior knowledge."
  #
  #   CONTINUOUS (truly numeric) variables like wealth_index_score, travel time
  #     †’ Cluster MEAN is appropriate.
  #     We also compute a z-score (standardised) version so that odds ratios
  #     represent "per 1-SD increase" for comparability.
  #     Column naming: cl_<varname> (mean) and cl_<varname>_z (z-score)
  #
  # WHY THIS MATTERS:
  #   If we just average education codes (0,1,2,3), a cluster with 50% "None"
  #   and 50% "Higher" gets the same mean (1.5) as a cluster with 100%
  #   "Primary".  These are obviously very different populations.  By creating
  #   separate proportion columns for each level, we preserve this information
  #   and can ask: "does a 10% increase in % with primary education predict
  #   coverage differently from a 10% increase in % with higher education?"
  # ═════════════════════════════════════════════════════════════════
  
  cat("\n----- SECTION 5: Aggregating to cluster level -----\n")
  
  # ----- 5a. Identify which variables are present in the data --------------------------------------------------”€
  all_cov <- c(individual_vars, household_vars, community_vars, travel_time_vars)
  all_cov_in_data <- intersect(c(all_cov), names(model_data))
  cat("Covariates found in data:", length(all_cov_in_data),
      "of", length(all_cov), "candidates\n")
  
  # Classify the variables that are actually present in the data
  binary_in_data     <- intersect(binary_vars, names(model_data))
  categorical_in_data <- intersect(categorical_vars, names(model_data))
  continuous_in_data  <- intersect(c(continuous_vars, travel_time_vars), names(model_data))
  # gc_in_data          <- gc_available  # already cluster-level, just need z-scores
  
  cat("  Binary in data:     ", length(binary_in_data), "\n")
  cat("  Categorical in data:", length(categorical_in_data), "\n")
  cat("  Continuous in data: ", length(continuous_in_data), "\n")
 #  cat("  GC in data:         ", length(gc_in_data), "\n")
  
  
  # ----- 5b. Aggregate: outcome + binary variables + continuous variables --------------------”€
  #
  # Binary: mean of 0/1 = proportion (this is correct and interpretable)
  # Continuous: mean across children in the cluster
  cat("\n----- 5b: Aggregating binary and continuous covariates -----\n")
  
  # Combine binary + continuous for the simple mean aggregation
  simple_mean_vars <- intersect(c(binary_in_data, continuous_in_data),
                                names(model_data))
  
  
  cluster_data <- model_data %>%
    group_by(CountryName, country_phase, cluster_psu, SurveyYear, LATNUM, LONGNUM) %>%
    summarise(
      # Outcome
      n_children   = n(),
      n_vaccinated = sum(mcv1, na.rm = TRUE),
      coverage     = n_vaccinated / n_children,
      
      # Binary variables †’ proportion (mean of 0/1)
      # Continuous variables †’ cluster mean
      across(all_of(simple_mean_vars), ~ mean(.x, na.rm = TRUE), .names = "cl_{.col}"),
      
      .groups = "drop"
    ) %>%
    # Empirical logit (Haldane-Anscombe correction):
    # Handles logit(0) = -Inf and logit(1) = +Inf by adding 0.5
    mutate(
      coverage_emp_logit = log((n_vaccinated + 0.5) / (n_children - n_vaccinated + 0.5))
    )
  
  cat("Clusters after binary/continuous aggregation:", nrow(cluster_data), "\n")
  
  
  # ----- 5c. Aggregate: categorical variables †’ per-level proportions ------------------------------”€
  #
  # For each categorical variable (e.g. respondent_edu_level with codes 0,1,2,3),
  # we compute the PROPORTION of children in each level within each cluster.
  # This creates columns like:
  #   cl_edu_primary       = fraction of mothers with primary education
  #   cl_edu_secondary     = fraction of mothers with secondary education
  #   cl_edu_higher        = fraction of mothers with higher education
  #   (reference level "no_education" is dropped  its proportion = 1 minus
  #    the sum of the others)
  #
  # WHY DROP THE REFERENCE LEVEL?
  # Because the proportions sum to 1 within each cluster, including ALL levels
  # would create perfect multicollinearity (the last column is exactly
  # determined by the others).  This is the same reason you drop one dummy
  # variable in standard regression.  Following Utazi et al. (2022), the
  # reference level is typically the most common category.
  
  cat("\n----- 5c: Aggregating categorical covariates †’ per-level proportions -----\n")
  
  prop_cols_created <- list()  # track which columns we created for each variable
  
  for (cv in categorical_in_data) {
    x <- model_data[[cv]]
    if (all(is.na(x))) {
      cat("  ", cv, ": all NA  skipping\n")
      next
    }
    
    # Find the levels present in the data
    levs <- sort(unique(x[!is.na(x)]))
    if (length(levs) < 2) {
      cat("  ", cv, ": only 1 level  skipping\n")
      next
    }
    
    # Identify the reference level (most common = drop)
    #tab <- sort(table(x, useNA = "no"), decreasing = TRUE)
    ref_level <- min(levs)
    
    # Get human-readable labels for this variable
    labels <- categorical_level_labels[[cv]]
    
    # Create proportion columns for each NON-REFERENCE level
    non_ref_levs <- setdiff(as.character(levs), ref_level)
    created_cols <- c()
    
    for (lv in non_ref_levs) {
      # Create a readable column name
      if (!is.null(labels) && lv %in% names(labels)) {
        col_name <- paste0("cl_", cv, "_", labels[lv])
      } else {
        col_name <- paste0("cl_", cv, "_lev", lv)
      }
      # Make syntactically valid R name
      col_name <- make.names(col_name)
      
      # Compute per-cluster proportion for this level
      lv_props <- model_data %>%
        group_by(cluster_psu) %>%
        summarise(!!col_name := mean(.data[[cv]] == as.numeric(lv), na.rm = TRUE),
                  .groups = "drop")
      
      cluster_data <- left_join(cluster_data, lv_props, by = "cluster_psu")
      created_cols <- c(created_cols, col_name)
    }
    
    # Also create the reference level proportion for descriptive purposes
    # (not used in the model, but useful for summary tables)
    if (!is.null(labels) && ref_level %in% names(labels)) {
      ref_col_name <- paste0("cl_", cv, "_", labels[1], "_REF")
    } else {
      ref_col_name <- paste0("cl_", cv, "_lev", ref_level, "_REF")
    }
    ref_col_name <- make.names(ref_col_name)
    ref_props <- model_data %>%
      group_by(cluster_psu) %>%
      summarise(!!ref_col_name := mean(.data[[cv]] == as.numeric(ref_level), na.rm = TRUE),
                .groups = "drop")
    cluster_data <- left_join(cluster_data, ref_props, by = "cluster_psu")
    
    prop_cols_created[[cv]] <- list(
      ref_level      = ref_level,
      ref_label      = if (!is.null(labels) && ref_level %in% names(labels)) labels[ref_level] else ref_level,
      ref_col        = ref_col_name,
      non_ref_cols   = created_cols,
      all_levels     = as.character(levs)
    )
    
    cat("  ", cv, ": ref='", ref_col_name, "' (dropped) | ",
        length(created_cols), " proportion cols created\n", sep = "")
  }
  
  cat("\nTotal categorical proportion columns created:",
      sum(sapply(prop_cols_created, function(x) length(x$non_ref_cols))), "\n")
  
  
  # ----- 5d. Remove the old misleading mean-of-code columns for categorical vars ”€
  #
  # Now that we have proper per-level proportions, REMOVE the old cl_<varname>
  # columns for categorical variables  these were the misleading mean-of-codes.
  # For example, cl_respondent_edu_level (mean of 0,1,2,3 codes) is now
  # replaced by cl_respondent_edu_level_primary, cl_..._secondary, etc.
  
  cat("\n----- 5d: Removing misleading mean-of-code columns for categorical vars -----\n")
  
  old_cat_cols <- paste0("cl_", categorical_in_data)
  old_cat_present <- intersect(old_cat_cols, names(cluster_data))
  
  if (length(old_cat_present) > 0) {
    cluster_data <- cluster_data %>% dplyr::select(-all_of(old_cat_present))
    cat("  Removed", length(old_cat_present), "mean-of-code columns:",
        paste(head(old_cat_present, 5), collapse = ", "),
        if (length(old_cat_present) > 5) "..." else "", "\n")
  }
  
  
  # ----- 5e. Standardise continuous covariates (z-scores) ------------------------------------------------------------”€
  #
  # Z-scores put all continuous variables on the same scale (mean=0, SD=1).
  # This is useful because:
  #   - Odds ratios become "per 1-SD increase" which is more interpretable
  #     than "per 1-unit increase" when units differ (e.g. years vs km)
  #   - It helps with model convergence in INLA
  #   - It allows direct comparison of effect sizes across covariates
  #
  # We standardise ACROSS ALL COUNTRIES here (global z-scores).
  # NOTE: we only z-score the truly continuous variables and GC variables,
  # NOT the binary proportions or categorical proportions (which are already
  # on a 0-1 scale and directly interpretable as "per 10% increase").
  
  cat("\n----- 5e: Standardising continuous covariates (z-scores) -----\n")
  
  # Identify continuous cl_ columns (NOT binary proportions, NOT categorical proportions)
  continuous_cl_cols <- c(
    paste0("cl_", intersect(continuous_in_data, names(model_data)))
  )
  # Fix: gc columns are named cl_gc_... from the across() call
  continuous_cl_cols <- intersect(
    grep("^cl_", names(cluster_data), value = TRUE),
    c(paste0("cl_", continuous_in_data),
      grep("^cl_gc_", names(cluster_data), value = TRUE),
      grep("^cl_travel_", names(cluster_data), value = TRUE))
  )
  
  for (v in continuous_cl_cols) {
    v_z <- paste0(v, "_z")
    x <- cluster_data[[v]]
    s <- sd(x, na.rm = TRUE)
    if (!is.na(s) && s > 1e-10) {
      cluster_data[[v_z]] <- (x - mean(x, na.rm = TRUE)) / s
    }
  }
  
  cat("Added", sum(grepl("_z$", names(cluster_data))), "z-score columns for continuous vars\n")
  
  
  # ----- 5f. Final cluster-level dataset summary -------------------------------------------------------------------------------------
  
  cat("\n----- 5f: Cluster-level dataset summary -----\n")
  
  cl_vars <- grep("^cl_", names(cluster_data), value = TRUE)
  cl_vars_model <- cl_vars[!grepl("_REF$", cl_vars)]  # exclude reference-level cols
  
  cat("Cluster-level dataset:", nrow(cluster_data), "clusters across",
      n_distinct(cluster_data$CountryName), "countries\n")
  cat("Total cl_ columns:", length(cl_vars), "\n")
  cat("  Model covariates (excl. reference):", length(cl_vars_model), "\n")
  cat("  Reference-level columns (descriptive only):", sum(grepl("_REF$", cl_vars)), "\n")
  
  # Print a summary of column types
  cat("\n  Binary proportion cols (mean of 0/1):\n")
  bin_cls <- paste0("cl_", binary_in_data)
  bin_cls_present <- intersect(bin_cls, names(cluster_data))
  cat("    ", paste(head(bin_cls_present, 6), collapse = ", "),
      if (length(bin_cls_present) > 6) paste0(" ... (", length(bin_cls_present), " total)") else "", "\n")
  
  cat("  Categorical proportion cols (per-level %):\n")
  cat_cls <- unlist(lapply(prop_cols_created, `[[`, "non_ref_cols"))
  cat("    ", paste(head(cat_cls, 6), collapse = ", "),
      if (length(cat_cls) > 6) paste0(" ... (", length(cat_cls), " total)") else "", "\n")
  
  cat("  Continuous cols (cluster mean):\n")
  cat("    ", paste(head(continuous_cl_cols, 6), collapse = ", "),
      if (length(continuous_cl_cols) > 6) paste0(" ... (", length(continuous_cl_cols), " total)") else "", "\n")
  
  
  # ----- 5g. Save cluster-level dataset ---------------------------------------------------------------------------------------------------------
  # data.table::fwrite(cluster_data, paste0(base_path, "/Data/DHS_final/cluster_data_v2.csv"))
  # cluster_data <- vroom::vroom(paste0(base_path, "/Data/DHS_final/cluster_data_v2.csv"))
  # 
  
  # ═════════════════════════════════════════════════════════════════
  # SECTION 6: DESCRIPTIVE SUMMARIES OF CLUSTER-LEVEL COVARIATE DISTRIBUTIONS
  # ═════════════════════════════════════════════════════════════════
  #
  # Now that we have proper cluster-level proportions, we print descriptive
  # summaries showing what a "typical cluster" looks like in the data.
  # This helps contextualise the model results: if 80% of clusters have
  # <10% with higher education, the model can't tell us much about the
  # effect of higher education.
  #
  # These summaries are at the CLUSTER LEVEL (not individual level), because
  # the cluster is the unit of analysis in the model.
  # ═════════════════════════════════════════════════════════════════
  
  cat("\n----- SECTION 6: Cluster-level covariate distribution summaries -----\n")
  
  # Print summary of binary proportions across clusters
  cat("\n----- Binary covariates (mean proportion across all clusters) -----\n")
  cat("  (Read as: 'on average, X% of individuals within a cluster have this trait')\n\n")
  
  for (bv in bin_cls_present) {
    x <- cluster_data[[bv]]
    if (!all(is.na(x))) {
      cat(sprintf("  %-45s mean=%.1f%%  median=%.1f%%  SD=%.1f%%\n",
                  bv,
                  mean(x, na.rm = TRUE) * 100,
                  median(x, na.rm = TRUE) * 100,
                  sd(x, na.rm = TRUE) * 100))
    }
  }
  
  # Print summary of categorical proportions across clusters
  cat("\n----- Categorical covariates (mean proportion across all clusters) -----\n")
  cat("  (Read as: 'on average, X% of individuals within a cluster are in this category')\n\n")
  
  for (cv in names(prop_cols_created)) {
    info <- prop_cols_created[[cv]]
    cat("  ", cv, " (reference = '", info$ref_label, "', dropped):\n", sep = "")
    
    # Print reference level mean
    if (info$ref_col %in% names(cluster_data)) {
      ref_x <- cluster_data[[info$ref_col]]
      if (!all(is.na(ref_x)))
        cat(sprintf("    %-42s mean=%.1f%% [REF  not in model]\n",
                    info$ref_col, mean(ref_x, na.rm = TRUE) * 100))
    }
    # Print non-reference levels
    for (nc in info$non_ref_cols) {
      if (nc %in% names(cluster_data)) {
        x <- cluster_data[[nc]]
        if (!all(is.na(x)))
          cat(sprintf("    %-42s mean=%.1f%%  SD=%.1f%%\n",
                      nc, mean(x, na.rm = TRUE) * 100, sd(x, na.rm = TRUE) * 100))
      }
    }
  }
  
  # Print summary of continuous variables across clusters
  cat("\n----- Continuous covariates (distribution across clusters) -----\n\n")
  for (cv in continuous_cl_cols) {
    if (cv %in% names(cluster_data)) {
      x <- cluster_data[[cv]]
      if (!all(is.na(x)))
        cat(sprintf("  %-45s mean=%.2f  SD=%.2f  range=[%.2f, %.2f]\n",
                    cv, mean(x, na.rm = TRUE), sd(x, na.rm = TRUE),
                    min(x, na.rm = TRUE), max(x, na.rm = TRUE)))
    }
  }
  
  
  # ═════════════════════════════════════════════════════════════════
  # SECTION 7: COVERAGE-BY-CLUSTER-SIZE TABLE (0/1 BOUNDARY DIAGNOSTICS)
  # ═════════════════════════════════════════════════════════════════
  #
  # This table shows why so many clusters have 0% or 100% coverage.
  # INTERPRETATION: clusters with few children (1-3) often show extreme
  # coverage simply because the denominator is too small.  For example,
  # if a cluster has only 2 children, coverage can only be 0%, 50%, or 100%.
  # This is a SAMPLING ARTEFACT, not a real signal.  The binomial likelihood
  # handles this naturally because 0/2 and 2/2 are perfectly valid binomial
  # outcomes  the model knows that a cluster with 2/2 vaccinated is NOT the
  # same as a cluster with 20/20.
  # ═════════════════════════════════════════════════════════════════
  
  cat("\n----- SECTION 7: Coverage by cluster size (0/1 boundary diagnostic) -----\n")
  
  cat("\nOverall coverage distribution:\n")
  cat("  Exactly 0%:       ", round(mean(cluster_data$coverage == 0) * 100, 1), "%\n")
  cat("  Exactly 100%:     ", round(mean(cluster_data$coverage == 1) * 100, 1), "%\n")
  cat("  Between 0% & 100%:", round(mean(cluster_data$coverage > 0 &
                                           cluster_data$coverage < 1) * 100, 1), "%\n")
  
  cov_by_size <- cluster_data %>%
    mutate(sz = cut(n_children, c(0, 3, 5, 10, 20, Inf),
                    labels = c("1-3", "4-5", "6-10", "11-20", "21+"))) %>%
    group_by(sz) %>%
    summarise(
      n_clusters    = n(),
      pct_zero      = round(mean(coverage == 0) * 100, 1),
      pct_one       = round(mean(coverage == 1) * 100, 1),
      pct_extreme   = round(mean(coverage == 0 | coverage == 1) * 100, 1),
      mean_coverage = round(mean(coverage) * 100, 1),
      .groups       = "drop"
    )
  
  cat("\nCoverage extremes by cluster size:\n")
  print(cov_by_size, n = Inf)
  
  cat("\nINTERPRETATION:\n")
  cat("  Small clusters (1-3 children) will often show 0% or 100%\n")
  cat("  simply because the denominator is too small for intermediate values.\n")
  cat("  The binomial likelihood handles this: 0/3 and 3/3 are valid.\n")
  cat("  The beta-binomial (Model 4) tests if overdispersion remains.\n")
  
  
  # ═════════════════════════════════════════════════════════════════
  # SECTION 8: COUNTRY-BY-COUNTRY MODELLING LOOP
  # ═════════════════════════════════════════════════════════════════
  #
  # For each country with ‰¥ 50 geo-located clusters, we run the full
  # Utazi et al. (2022) pipeline:
  #
  #   Step 1: Missingness screening (<5% rule)
  #   Step 2: Bivariate screening (crude ORs, keep p < 0.2)
  #   Step 3: Multicollinearity check (|r| > 0.8 then iterative GVIF)
  #   Step 4: Model 1  covariates-only binomial GLM
  #   Step 5: Model 2  spatial-only INLA-SPDE
  #   Step 6: Model 3  full INLA-SPDE (covariates + spatial + nugget)
  #   Step 7: Model 4  beta-binomial INLA-SPDE
  #   Step 8: Model comparison + diagnostic plots + Excel output
  # ═════════════════════════════════════════════════════════════════
  
  cat("\n----- SECTION 8: Beginning country-by-country modelling loop -----\n")
  
  # Identify countries with enough data
  countries <- cluster_data %>%
    count(CountryName) %>%
    filter(n >= 50) %>%
    arrange(desc(n)) %>%
    pull(CountryName)
  
  cat("\n", length(countries), "countries with >= 50 geo-located clusters:\n")
  cat(paste("  ", countries), sep = "\n")
  

  
}






# ═════════════════════════════════════════════════════════════════•════
# ═══  HELPER FUNCTION: Prediction metrics                                      ═══
# ═════════════════════════════════════════════════════════════════•════
#
# Computes three numbers to judge how well a model predicts coverage:
#   r   = Pearson correlation (do predicted and observed track together?)
#   R²  = proportion of variance explained (0 = useless, 1 = perfect)
#   MAE = mean absolute error in coverage units (e.g. 0.10 = 10pp off)


safe_file_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  substr(x, 1, 80)
}

safe_sheet_name <- function(x, used = character()) {
  x <- gsub("[][\\:*?/]+", "_", x)
  x <- substr(x, 1, 31)
  base <- x
  i <- 1
  while (x %in% used) {
    suffix <- paste0("_", i)
    x <- paste0(substr(base, 1, 31 - nchar(suffix)), suffix)
    i <- i + 1
  }
  x
}

make_residual_covariate_plots <- function(data, vars, country_name) {
  resid_col <- if ("resid_full" %in% names(data)) "resid_full" else "resid_best"
  plot_vars <- vars[vars %in% names(data)]
  plot_vars <- plot_vars[sapply(data[plot_vars], function(x) sd(x, na.rm = TRUE) > 1e-8)]
  if (length(plot_vars) == 0 || !(resid_col %in% names(data))) return(NULL)
  purrr::map(plot_vars, function(v) {
    ggplot(data, aes(x = .data[[v]], y = .data[[resid_col]])) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_point(aes(size = n_children), alpha = 0.35, colour = "steelblue") +
      geom_smooth(method = "loess", se = TRUE, colour = "firebrick", linewidth = 0.7) +
      scale_size_continuous(range = c(0.6, 3), guide = "none") +
      labs(
        title = paste(country_name, "- residuals vs", v),
        subtitle = "Loess pattern suggests possible non-linearity or remaining structure",
        x = v, y = "Residual (observed - predicted)"
      ) +
      theme_minimal(base_size = 9)
  })
}


raw_covariate_name <- function(v) {
  ifelse(grepl("_z$", v), sub("_z$", "", v), v)
}


covariate_map_info <- function(v, data_names) {
  raw_v <- raw_covariate_name(v)
  use_raw <- grepl("_z$", v) && raw_v %in% data_names
  map_var <- if (use_raw) raw_v else v
  is_z <- grepl("_z$", v) && !use_raw
  
  label <- map_var %>%
    gsub("^cl_", "", .) %>%
    gsub("^gc_", "", .) %>%
    gsub("_", " ", .)
  label <- tools::toTitleCase(label)
  
  legend <- case_when(
    grepl("travel_time", map_var, ignore.case = TRUE) ~ "Minutes",
    grepl("elevation", map_var, ignore.case = TRUE) ~ "Metres",
    grepl("nightlights", map_var, ignore.case = TRUE) ~ "Nightlights",
    grepl("_prop|_primary|_secondary|_rich|_middle|_skilled|media|edu|wealth|household",
          map_var, ignore.case = TRUE) ~ "Proportion",
    is_z ~ "Standardised z-score",
    TRUE ~ "Value"
  )
  
  subtitle <- if (is_z) {
    paste0("Mapped as z-score because raw column ", raw_v, " is not available here.")
  } else if (use_raw) {
    paste0("Model used ", v, "; map shows original-unit column ", raw_v, ".")
  } else {
    "Mapped on the same scale used in the model."
  }
  
  list(
    model_var = v,
    map_var = map_var,
    label = label,
    legend = legend,
    subtitle = subtitle
  )
}


needs_percentile_truncation <- function(v) {
  grepl(
    paste(c(
      "travel_time_to_HC_walk",
      "travel_time_to_city",
      "travel_time_to_HC_motor",
      "respondent_occupation_professional",
      "nightlights_composite"
    ), collapse = "|"),
    v,
    ignore.case = TRUE
  )
}

get_plot_limits <- function(x, v) {
  if (needs_percentile_truncation(v)) {
    quantile(x, c(0.02, 0.98), na.rm = TRUE)
  } else {
    range(x, na.rm = TRUE)
  }
}

limit_note <- function(v) {
  if (needs_percentile_truncation(v)) {
    "Colour scale truncated at 2nd-98th percentiles to reduce extreme-value dominance."
  } else {
    NULL
  }
}

make_covariate_histograms <- function(sub_m, vars_final) {
  map_infos <- purrr::map(vars_final, covariate_map_info, data_names = names(sub_m))
  
  purrr::map(map_infos, function(info) {
    x <- sub_m[[info$map_var]]
    
    lims <- get_plot_limits(x, info$map_var)
    if (any(!is.finite(lims)) || lims[1] == lims[2]) lims <- NULL
    
    subtitle <- limit_note(info$map_var)
    if (is.null(subtitle)) subtitle <- "Cluster-level distribution."
    
    ggplot(sub_m, aes(x = .data[[info$map_var]])) +
      geom_histogram(bins = 40, fill = "grey55", colour = "white") +
      coord_cartesian(xlim = lims) +
      labs(
        title = paste("Histogram:", info$label),
        subtitle = subtitle,
        x = info$legend,
        y = "Clusters"
      ) +
      theme_minimal(base_size = 9)
  })
}

make_partial_residual_plots <- function(sub_m, vars_final, feff = NULL, glm_co = NULL) {
  resid_col <- if ("resid_full" %in% names(sub_m)) "resid_full" else "resid_best"
  
  coef_lookup <- NULL
  
  if (!is.null(feff)) {
    coef_lookup <- feff %>%
      dplyr::select(term, beta = mean)
  } else if (!is.null(glm_co)) {
    coef_lookup <- glm_co %>%
      dplyr::select(term, beta = Estimate)
  } else {
    return(NULL)
  }
  
  map_infos <- purrr::map(vars_final, covariate_map_info, data_names = names(sub_m))
  
  purrr::map(map_infos, function(info) {
    coef_row <- coef_lookup %>% filter(term == info$model_var)
    if (nrow(coef_row) != 1) return(NULL)
    
    beta <- coef_row$beta[1]
    x_model <- sub_m[[info$model_var]]
    x_plot <- sub_m[[info$map_var]]
    
    partial_resid <- sub_m[[resid_col]] + beta * x_model
    
    lims <- get_plot_limits(x_plot, info$map_var)
    if (any(!is.finite(lims)) || lims[1] == lims[2]) lims <- NULL
    
    subtitle <- limit_note(info$map_var)
    if (is.null(subtitle)) subtitle <- "Loess curve checks for non-linearity after accounting for the fitted linear effect."
    
    ggplot(
      tibble(x_plot = x_plot, partial_resid = partial_resid, n_children = sub_m$n_children),
      aes(x = x_plot, y = partial_resid)
    ) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
      geom_point(aes(size = n_children), alpha = 0.35, colour = "steelblue") +
      geom_smooth(method = "loess", se = TRUE, colour = "firebrick", linewidth = 0.7) +
      coord_cartesian(xlim = lims) +
      scale_size_continuous(range = c(0.6, 3), guide = "none") +
      labs(
        title = paste("Partial residual:", info$label),
        subtitle = subtitle,
        x = info$legend,
        y = "Partial residual"
      ) +
      theme_minimal(base_size = 9)
  }) %>%
    purrr::compact()
}




# make_country_surface_diagnostics <- function(country_name, adm1_sf, sub_m, vars_final,
#                                              inla_fu, mesh, spde, lcrs, cell_km = 5) {
#   if (is.null(adm1_sf) || is.null(inla_fu) || is.null(mesh) || is.null(spde) || is.null(lcrs)) {
#     return(NULL)
#   }
#   
#   country_poly_ll <- tryCatch(
#     adm1_sf %>%
#       st_transform(4326) %>%
#       st_union() %>%
#       st_as_sf(),
#     error = function(e) NULL
#   )
#   if (is.null(country_poly_ll)) return(NULL)
#   
#   country_poly_km <- tryCatch(
#     st_transform(country_poly_ll, crs = lcrs),
#     error = function(e) NULL
#   )
#   if (is.null(country_poly_km)) return(NULL)
#   
#   grid_km <- tryCatch({
#     g <- st_make_grid(
#       country_poly_km,
#       cellsize = cell_km,
#       square = TRUE,
#       what = "polygons"
#     )
#     g <- st_as_sf(g) %>%
#       st_filter(country_poly_km, .predicate = st_intersects)
#     g$cell_id <- seq_len(nrow(g))
#     g
#   }, error = function(e) NULL)
#   
#   if (is.null(grid_km) || nrow(grid_km) < 10) return(NULL)
#   
#   grid_cent_km <- st_centroid(grid_km)
#   grid_cent_ll <- st_transform(grid_cent_km, 4326)
#   
#   grid_ll <- st_coordinates(grid_cent_ll)
#   grid_cp <- st_coordinates(grid_cent_km)
#   
#   A_grid <- INLA::inla.spde.make.A(mesh = mesh, loc = grid_cp)
#   spatial_mean <- as.vector(A_grid %*% inla_fu$summary.random$spatial$mean)
#   
#   grid_df <- tibble(
#     cell_id = grid_km$cell_id,
#     LONGNUM = grid_ll[, 1],
#     LATNUM = grid_ll[, 2],
#     spatial_field = spatial_mean
#   )
#   
#   map_infos <- purrr::map(vars_final, covariate_map_info, data_names = names(sub_m))
#   needed_cols <- unique(c(vars_final, purrr::map_chr(map_infos, "map_var")))
#   needed_cols <- needed_cols[needed_cols %in% names(sub_m)]
#   
#   obs_sp <- sub_m %>%
#     dplyr::select(LONGNUM, LATNUM, all_of(needed_cols))
#   
#   for (v in needed_cols) {
#     obs_v <- obs_sp %>%
#       dplyr::select(LONGNUM, LATNUM, all_of(v)) %>%
#       drop_na()
#     
#     names(obs_v)[3] <- "value"
#     
#     if (nrow(obs_v) >= 5 && sd(obs_v$value, na.rm = TRUE) > 1e-8) {
#       idw_fit <- tryCatch(
#         gstat::idw(
#           value ~ 1,
#           locations = ~ LONGNUM + LATNUM,
#           data = obs_v,
#           newdata = as.data.frame(grid_df[, c("LONGNUM", "LATNUM")]),
#           idp = 2
#         ),
#         error = function(e) NULL
#       )
#       
#       grid_df[[v]] <- if (!is.null(idw_fit)) {
#         idw_fit$var1.pred
#       } else {
#         median(obs_v$value, na.rm = TRUE)
#       }
#     } else if (nrow(obs_v) > 0) {
#       grid_df[[v]] <- median(obs_v$value, na.rm = TRUE)
#     } else {
#       grid_df[[v]] <- NA_real_
#     }
#   }
#   
#   fix <- inla_fu$summary.fixed$mean
#   names(fix) <- rownames(inla_fu$summary.fixed)
#   
#   lp <- rep(
#     ifelse("intercept" %in% names(fix), fix[["intercept"]], 0),
#     nrow(grid_df)
#   )
#   
#   for (v in vars_final) {
#     if (v %in% names(fix) && v %in% names(grid_df)) {
#       lp <- lp + fix[[v]] * grid_df[[v]]
#     }
#   }
#   
#   grid_df$predicted_coverage <- plogis(lp + grid_df$spatial_field)
#   
#   grid_map_km <- left_join(grid_km, grid_df, by = "cell_id")
#   boundary_km <- st_boundary(country_poly_km)
#   
#   map_theme <- theme_minimal(base_size = 9) +
#     theme(
#       axis.title = element_blank(),
#       panel.grid = element_line(colour = "grey88", linewidth = 0.2),
#       legend.position = "right",
#       plot.title = element_text(face = "bold", size = 10),
#       plot.subtitle = element_text(size = 7, colour = "grey35")
#     )
#   
#   p_pred <- ggplot() +
#     geom_sf(data = grid_map_km, aes(fill = predicted_coverage), colour = NA) +
#     geom_sf(data = boundary_km, colour = "grey25", linewidth = 0.25) +
#     scale_fill_viridis_c(
#       name = "Predicted\ncoverage",
#       limits = c(0, 1)
#     ) +
#     labs(
#       title = paste(country_name, "- predicted MCV1 coverage"),
#       subtitle = paste0(cell_km, " km x ", cell_km, " km grid; full fitted model")
#     ) +
#     coord_sf(datum = NA) +
#     map_theme
#   
#   p_spatial <- ggplot() +
#     geom_sf(data = grid_map_km, aes(fill = spatial_field), colour = NA) +
#     geom_sf(data = boundary_km, colour = "grey25", linewidth = 0.25) +
#     scale_fill_gradient2(
#       name = "Spatial\nfield",
#       low = "#d73027",
#       mid = "white",
#       high = "#4575b4",
#       midpoint = 0
#     ) +
#     labs(
#       title = paste(country_name, "- posterior mean spatial field"),
#       subtitle = paste0(cell_km, " km x ", cell_km, " km grid; logit scale")
#     ) +
#     coord_sf(datum = NA) +
#     map_theme
#   
#   covariate_plots <- purrr::map(map_infos, function(info) {
#     d <- grid_map_km
#     d$map_value <- d[[info$map_var]]
#     
#     is_prop <- identical(info$legend, "Proportion")
#     
#     lims <- get_plot_limits(d$map_value, info$map_var)
#     if (any(!is.finite(lims)) || lims[1] == lims[2]) lims <- NULL
#     
#     subtitle <- paste(
#       info$subtitle,
#       limit_note(info$map_var),
#       sep = ifelse(is.null(limit_note(info$map_var)), "", " ")
#     )
#     
#     ggplot() +
#       geom_sf(data = d, aes(fill = map_value), colour = NA) +
#       geom_sf(data = boundary_km, colour = "grey25", linewidth = 0.25) +
#       scale_fill_viridis_c(
#         name = info$legend,
#         limits = if (is_prop && is.null(lims)) c(0, 1) else lims,
#         labels = if (is_prop) scales::label_percent(accuracy = 1) else waiver(),
#         oob = scales::squish,
#         na.value = "grey90"
#       ) +
#       labs(
#         title = info$label,
#         subtitle = subtitle
#       ) +
#       coord_sf(datum = NA) +
#       map_theme
#   })
#   
#   names(covariate_plots) <- vars_final
#   
#   list(
#     grid = st_drop_geometry(grid_map_km),
#     surface_plot = p_pred / p_spatial,
#     covariate_plots = covariate_plots
#   )
# }
# 


# ════════════════════════════════════════════════════════════════════════════════
# REPLACEMENT: make_country_surface_diagnostics()
# ════════════════════════════════════════════════════════════════════════════════
# REPLACES: lines 1271–1460 in your original script
# (from "make_country_surface_diagnostics <- function(" to the closing "}" 
#  just before "get_spatial_var <- function(")
#
# TWO METHODS included — toggle via the 'method' argument:
#   method = "posterior_sample"  ← DEFAULT, RECOMMENDED
#       Uses inla.posterior.sample() to draw from the joint posterior.
#       Covariates at grid cells are interpolated via ordinary kriging.
#       Predictions = plogis(Xβ + S(s)), NO nugget.
#       Automatically produces uncertainty (SD, 95% CrI width).
#
#   method = "prediction_stack"
#       Builds a second INLA stack for prediction locations and re-runs
#       INLA with both estimation + prediction stacks. Covariates at grid
#       cells are interpolated via ordinary kriging. cluster_id = NA so
#       the nugget is excluded. INLA returns fitted values + uncertainty
#       at grid cells directly.
#       SLOWER (re-fits the model) but fully Bayesian.
#
# BOTH methods replace the old IDW interpolation with ordinary kriging
# for DHS-derived covariates, giving smooth covariate surfaces.
# ════════════════════════════════════════════════════════════════════════════════

make_country_surface_diagnostics <- function(country_name, adm1_sf, sub_m, vars_final,
                                             inla_fu, mesh, spde, lcrs,
                                             stk_est = NULL,    # needed for prediction_stack method
                                             f3 = NULL,         # needed for prediction_stack method
                                             iid_spec = NULL,   # needed for prediction_stack method
                                             cell_km = 5,
                                             n_posterior_samples = 500,
                                             method = "posterior_sample") {
  
  # ---- Guard clauses ----
  if (is.null(adm1_sf) || is.null(inla_fu) || is.null(mesh) || is.null(spde) || is.null(lcrs)) {
    return(NULL)
  }
  
  # ---- Country boundary ----
  country_poly_ll <- tryCatch(
    adm1_sf %>% st_transform(4326) %>% st_union() %>% st_as_sf(),
    error = function(e) NULL
  )
  if (is.null(country_poly_ll)) return(NULL)
  
  country_poly_km <- tryCatch(
    st_transform(country_poly_ll, crs = lcrs),
    error = function(e) NULL
  )
  if (is.null(country_poly_km)) return(NULL)
  
  # ---- Build prediction grid ----
  grid_km <- tryCatch({
    g <- st_make_grid(country_poly_km, cellsize = cell_km, square = TRUE, what = "polygons")
    g <- st_as_sf(g) %>% st_filter(country_poly_km, .predicate = st_intersects)
    g$cell_id <- seq_len(nrow(g))
    g
  }, error = function(e) NULL)
  
  if (is.null(grid_km) || nrow(grid_km) < 10) return(NULL)
  
  grid_cent_km <- st_centroid(grid_km)
  grid_cent_ll <- st_transform(grid_cent_km, 4326)
  grid_ll <- st_coordinates(grid_cent_ll)
  grid_cp <- st_coordinates(grid_cent_km)
  
  n_grid <- nrow(grid_km)
  cat("    Prediction grid:", n_grid, "cells at", cell_km, "km resolution\n")
  
  # ---- Projection matrix for spatial field at grid locations ----
  A_grid <- INLA::inla.spde.make.A(mesh = mesh, loc = grid_cp)
  
  # ---- Initialise grid data frame ----
  grid_df <- tibble(
    cell_id = grid_km$cell_id,
    LONGNUM = grid_ll[, 1],
    LATNUM  = grid_ll[, 2]
  )
  
  
  # ══════════════════════════════════════════════════════════════════
  # INTERPOLATE COVARIATES TO GRID: Ordinary Kriging (replaces IDW)
  # ══════════════════════════════════════════════════════════════════
  #
  # For each covariate, we fit a variogram to the cluster-level values
  # and use ordinary kriging to predict at grid locations. This gives:
  #   (a) smooth surfaces (no bullseye/dotty artefacts)
  #   (b) kriging variance at each cell (used for covariate uncertainty maps)
  #
  # Falls back to IDW if kriging fails (e.g. too few points for variogram).
  
  cat("    Kriging covariates to prediction grid...\n")
  
  map_infos <- purrr::map(vars_final, covariate_map_info, data_names = names(sub_m))
  needed_cols <- unique(c(vars_final, purrr::map_chr(map_infos, "map_var")))
  needed_cols <- needed_cols[needed_cols %in% names(sub_m)]
  
  obs_sp <- sub_m %>% dplyr::select(LONGNUM, LATNUM, all_of(needed_cols))
  
  # Store kriging variance for covariate uncertainty maps
  krige_variance <- list()
  
  for (v in needed_cols) {
    obs_v <- obs_sp %>%
      dplyr::select(LONGNUM, LATNUM, all_of(v)) %>%
      drop_na()
    names(obs_v)[3] <- "value"
    
    if (nrow(obs_v) >= 30 && sd(obs_v$value, na.rm = TRUE) > 1e-8) {
      
      # Try ordinary kriging
      krige_result <- tryCatch({
        # Fit variogram
        coords_sf <- st_as_sf(obs_v, coords = c("LONGNUM", "LATNUM"), crs = 4326)
        coords_proj <- st_transform(coords_sf, crs = lcrs)
        coords_xy <- st_coordinates(coords_proj)
        obs_v_proj <- obs_v
        obs_v_proj$X <- coords_xy[, 1]
        obs_v_proj$Y <- coords_xy[, 2]
        
        emp_vgm <- gstat::variogram(value ~ 1, locations = ~ X + Y, data = obs_v_proj)
        fit_vgm <- gstat::fit.variogram(emp_vgm, gstat::vgm("Sph"))
        
        # Prediction locations in projected coords
        newdata_proj <- data.frame(X = grid_cp[, 1], Y = grid_cp[, 2])
        
        # Krige
        k <- gstat::krige(
          value ~ 1,
          locations = ~ X + Y,
          data = obs_v_proj,
          newdata = newdata_proj,
          model = fit_vgm
        )
        
        list(pred = k$var1.pred, var = k$var1.var)
      }, error = function(e) NULL)
      
      if (!is.null(krige_result)) {
        grid_df[[v]] <- krige_result$pred
        krige_variance[[v]] <- krige_result$var
        cat("      ", v, ": kriging OK\n")
      } else {
        # Fallback to IDW if kriging fails
        idw_fit <- tryCatch(
          gstat::idw(value ~ 1, locations = ~ LONGNUM + LATNUM,
                     data = obs_v,
                     newdata = as.data.frame(grid_df[, c("LONGNUM", "LATNUM")]),
                     idp = 2),
          error = function(e) NULL
        )
        grid_df[[v]] <- if (!is.null(idw_fit)) idw_fit$var1.pred else median(obs_v$value, na.rm = TRUE)
        krige_variance[[v]] <- NULL
        cat("      ", v, ": kriging failed, used IDW fallback\n")
      }
      
    } else if (nrow(obs_v) >= 5 && sd(obs_v$value, na.rm = TRUE) > 1e-8) {
      # Too few points for kriging, use IDW
      idw_fit <- tryCatch(
        gstat::idw(value ~ 1, locations = ~ LONGNUM + LATNUM,
                   data = obs_v,
                   newdata = as.data.frame(grid_df[, c("LONGNUM", "LATNUM")]),
                   idp = 2),
        error = function(e) NULL
      )
      grid_df[[v]] <- if (!is.null(idw_fit)) idw_fit$var1.pred else median(obs_v$value, na.rm = TRUE)
      krige_variance[[v]] <- NULL
      cat("      ", v, ": <30 pts, used IDW\n")
    } else if (nrow(obs_v) > 0) {
      grid_df[[v]] <- median(obs_v$value, na.rm = TRUE)
      krige_variance[[v]] <- NULL
    } else {
      grid_df[[v]] <- NA_real_
      krige_variance[[v]] <- NULL
    }
  }
  
  
  # ══════════════════════════════════════════════════════════════════
  # METHOD A: POSTERIOR SAMPLING (DEFAULT)
  # ══════════════════════════════════════════════════════════════════
  #
  # Draw n_posterior_samples from the joint posterior of all latent
  # effects and hyperparameters. For each draw, compute:
  #   linear_predictor = intercept + Σ(β_j × cov_j) + A_grid × spatial
  #   predicted_coverage = plogis(linear_predictor)
  #
  # The nugget (cluster_id iid) is NOT included because we are
  # predicting at new locations, not at observed clusters.
  #
  # From the matrix of samples we get mean, SD, and quantiles.
  # ══════════════════════════════════════════════════════════════════
  
  if (method == "posterior_sample") {
    
    cat("    Drawing", n_posterior_samples, "posterior samples...\n")
    
    samps <- tryCatch(
      INLA::inla.posterior.sample(n_posterior_samples, inla_fu),
      error = function(e) { cat("      Posterior sampling failed:", e$message, "\n"); NULL }
    )
    
    if (!is.null(samps)) {
      
      # Build the covariate design matrix at grid locations
      # Columns: intercept, then each covariate in vars_final
      X_grid <- matrix(1, nrow = n_grid, ncol = 1)  # intercept column
      colnames(X_grid) <- "intercept"
      
      for (v in vars_final) {
        if (v %in% names(grid_df)) {
          X_grid <- cbind(X_grid, grid_df[[v]])
          colnames(X_grid)[ncol(X_grid)] <- v
        }
      }
      
      # Pre-allocate prediction matrix: n_grid × n_samples
      pred_matrix <- matrix(NA_real_, nrow = n_grid, ncol = length(samps))
      
      # Identify which rows in the latent vector correspond to fixed effects
      # and to the spatial field
      latent_names <- rownames(samps[[1]]$latent)
      
      # Fixed effect names in INLA latent: "intercept:1", "cl_varname:1", etc.
      fixed_names_inla <- paste0(colnames(X_grid), ":1")
      
      # Spatial field indices
      spatial_idx <- grep("^spatial:", latent_names)
      
      for (s in seq_along(samps)) {
        lat <- samps[[s]]$latent[, 1]  # latent vector for this sample
        
        # Extract fixed effects for this sample
        beta_s <- lat[match(fixed_names_inla, latent_names)]
        
        # If any fixed effect not found, skip this sample
        if (any(is.na(beta_s))) {
          # Try without the :1 suffix (some INLA versions)
          beta_s <- lat[match(colnames(X_grid), latent_names)]
        }
        if (any(is.na(beta_s))) next
        
        # Extract spatial field and project to grid
        spatial_s <- lat[spatial_idx]
        spatial_at_grid <- as.vector(A_grid %*% spatial_s)
        
        # Linear predictor: X β + spatial (NO nugget)
        lp_s <- as.vector(X_grid %*% beta_s) + spatial_at_grid
        
        # Transform to probability scale
        pred_matrix[, s] <- plogis(lp_s)
      }
      
      # Remove any all-NA columns (failed samples)
      good_cols <- !apply(pred_matrix, 2, function(x) all(is.na(x)))
      pred_matrix <- pred_matrix[, good_cols, drop = FALSE]
      cat("      Usable samples:", ncol(pred_matrix), "of", length(samps), "\n")
      
      if (ncol(pred_matrix) >= 50) {
        grid_df$predicted_coverage <- rowMeans(pred_matrix, na.rm = TRUE)
        grid_df$pred_sd            <- apply(pred_matrix, 1, sd, na.rm = TRUE)
        grid_df$pred_lo            <- apply(pred_matrix, 1, quantile, 0.025, na.rm = TRUE)
        grid_df$pred_hi            <- apply(pred_matrix, 1, quantile, 0.975, na.rm = TRUE)
        grid_df$pred_iqr           <- grid_df$pred_hi - grid_df$pred_lo
        
        cat("      Mean predicted coverage:", round(mean(grid_df$predicted_coverage, na.rm = TRUE), 3), "\n")
        cat("      Mean uncertainty (SD):", round(mean(grid_df$pred_sd, na.rm = TRUE), 3), "\n")
      } else {
        cat("      Too few usable samples — falling back to point estimate.\n")
        method <- "point_fallback"
      }
      
      # Also extract the posterior mean spatial field for the spatial field map
      grid_df$spatial_field <- as.vector(A_grid %*% inla_fu$summary.random$spatial$mean)
      
    } else {
      cat("      Posterior sampling not available — falling back to point estimate.\n")
      method <- "point_fallback"
    }
    
    # Point estimate fallback (same logic as original but with kriged covariates)
    if (method == "point_fallback") {
      fix <- inla_fu$summary.fixed$mean
      names(fix) <- rownames(inla_fu$summary.fixed)
      
      lp <- rep(ifelse("intercept" %in% names(fix), fix[["intercept"]], 0), n_grid)
      for (v in vars_final) {
        if (v %in% names(fix) && v %in% names(grid_df)) {
          lp <- lp + fix[[v]] * grid_df[[v]]
        }
      }
      
      grid_df$spatial_field       <- as.vector(A_grid %*% inla_fu$summary.random$spatial$mean)
      grid_df$predicted_coverage  <- plogis(lp + grid_df$spatial_field)
      grid_df$pred_sd             <- NA_real_
      grid_df$pred_lo             <- NA_real_
      grid_df$pred_hi             <- NA_real_
      grid_df$pred_iqr            <- NA_real_
    }
  }
  
  
  # ══════════════════════════════════════════════════════════════════
  # METHOD B: PREDICTION STACK (comment in to use instead)
  # ══════════════════════════════════════════════════════════════════
  #
  # Builds a second INLA stack for the prediction grid, combines it
  # with the estimation stack, and re-runs INLA. INLA then directly
  # produces fitted values (with uncertainty) at grid locations.
  #
  # cluster_id = NA → nugget is NOT predicted at new locations.
  #
  # SLOWER than posterior sampling (re-fits the model) but the
  # uncertainty estimates are exact (no Monte Carlo noise).
  #
  # To use this method, you must pass stk_est, f3, and iid_spec
  # from the main loop into this function.
  # ══════════════════════════════════════════════════════════════════
  
  if (method == "prediction_stack") {
    
    cat("    Building prediction stack and re-running INLA...\n")
    
    if (is.null(stk_est) || is.null(f3) || is.null(iid_spec)) {
      cat("      Missing stk_est/f3/iid_spec — cannot use prediction_stack method.\n")
      return(NULL)
    }
    
    # Build covariate data frame for grid locations
    cov_grid <- data.frame(intercept = rep(1, n_grid))
    for (v in vars_final) {
      if (v %in% names(grid_df)) {
        cov_grid[[v]] <- grid_df[[v]]
      } else {
        cov_grid[[v]] <- NA_real_
      }
    }
    # cluster_id = NA → tells INLA "no nugget at these locations"
    cov_grid$cluster_id <- NA_integer_
    
    # Prediction stack
    stk_pred <- INLA::inla.stack(
      data    = list(y = rep(NA, n_grid), n_trial = rep(NA, n_grid)),
      A       = list(A_grid, 1),
      effects = list(
        spatial = 1:spde$n.spde,
        cov_grid
      ),
      tag = "prediction"
    )
    
    # Combine estimation + prediction
    stk_full <- INLA::inla.stack(stk_est, stk_pred)
    
    # Re-run INLA with the combined stack
    inla_pred <- tryCatch(
      INLA::inla(
        formula = f3,
        family  = "binomial",
        Ntrials = INLA::inla.stack.data(stk_full)$n_trial,
        data    = INLA::inla.stack.data(stk_full, spde = spde),
        control.predictor = list(
          A = INLA::inla.stack.A(stk_full), compute = TRUE, link = 1
        ),
        control.compute = list(dic = FALSE, waic = FALSE, config = FALSE),
        control.inla = list(strategy = "adaptive"),
        control.mode = list(theta = inla_fu$mode$theta, restart = FALSE),
        verbose = FALSE
      ),
      error = function(e) { cat("      Prediction INLA failed:", e$message, "\n"); NULL }
    )
    
    if (!is.null(inla_pred)) {
      # Extract predictions at grid locations
      idx_pred <- INLA::inla.stack.index(stk_full, "prediction")$data
      
      grid_df$predicted_coverage <- inla_pred$summary.fitted.values$mean[idx_pred]
      grid_df$pred_sd            <- inla_pred$summary.fitted.values$sd[idx_pred]
      grid_df$pred_lo            <- inla_pred$summary.fitted.values$`0.025quant`[idx_pred]
      grid_df$pred_hi            <- inla_pred$summary.fitted.values$`0.975quant`[idx_pred]
      grid_df$pred_iqr           <- grid_df$pred_hi - grid_df$pred_lo
      
      cat("      Mean predicted coverage:", round(mean(grid_df$predicted_coverage, na.rm = TRUE), 3), "\n")
      cat("      Mean uncertainty (SD):", round(mean(grid_df$pred_sd, na.rm = TRUE), 3), "\n")
    } else {
      return(NULL)
    }
    
    # Spatial field (posterior mean from the ORIGINAL fit, not the re-run)
    grid_df$spatial_field <- as.vector(A_grid %*% inla_fu$summary.random$spatial$mean)
  }
  
  
  # ══════════════════════════════════════════════════════════════════
  # BUILD PLOTS (same for both methods)
  # ══════════════════════════════════════════════════════════════════
  
  grid_map_km  <- left_join(grid_km, grid_df, by = "cell_id")
  boundary_km  <- st_boundary(country_poly_km)
  
  map_theme <- theme_minimal(base_size = 9) +
    theme(
      axis.title    = element_blank(),
      panel.grid    = element_line(colour = "grey88", linewidth = 0.2),
      legend.position = "right",
      plot.title    = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 7, colour = "grey35")
    )
  
  # ---- Predicted coverage surface ----
  p_pred <- ggplot() +
    geom_sf(data = grid_map_km, aes(fill = predicted_coverage), colour = NA) +
    geom_sf(data = boundary_km, colour = "grey25", linewidth = 0.25) +
    scale_fill_gradient2(
      name = "Predicted\ncoverage",
      low = "#b6382d",
      mid = "#f2eac6",
      high = "#496fa4",
      midpoint = 0.5,
      limits = c(0, 1)
    ) +
    labs(
      title    = paste(country_name, "- predicted MCV1 coverage"),
      subtitle = paste0(cell_km, " km grid; full model (no nugget)")
    ) +
    coord_sf(datum = NA) + map_theme
  
  # ---- Spatial field surface ----
  p_spatial <- ggplot() +
    geom_sf(data = grid_map_km, aes(fill = spatial_field), colour = NA) +
    geom_sf(data = boundary_km, colour = "grey25", linewidth = 0.25) +
    scale_fill_gradient2(
      name = "Spatial\nfield",
      low = "#b6382d",
      mid = "#f2eac6",
      high = "#496fa4",
      midpoint = 0,
      limits = c(-1, 1)
    ) +
    labs(
      title    = paste(country_name, "- posterior mean spatial field"),
      subtitle = paste0(cell_km, " km grid; logit scale")
    ) +
    coord_sf(datum = NA) + map_theme
  
  # ---- Uncertainty maps ----
  uncertainty_plots <- list()
  
  # MCV1 coverage uncertainty (SD)
  if ("pred_sd" %in% names(grid_df) && !all(is.na(grid_df$pred_sd))) {
    
    p_unc_sd <- ggplot() +
      geom_sf(data = grid_map_km, aes(fill = pred_sd), colour = NA) +
      geom_sf(data = boundary_km, colour = "grey25", linewidth = 0.25) +
      scale_fill_gradient2(
        name = "Posterior\nSD",
        low = "white",
        high = "#781707",
        limits = c(0, NA)
      ) +
      labs(
        title    = paste(country_name, "- MCV1 prediction uncertainty (SD)"),
        subtitle = paste0(cell_km, " km grid; higher = less certain")
      ) +
      coord_sf(datum = NA) + map_theme
    
    uncertainty_plots[["mcv1_sd"]] <- p_unc_sd
  }
  
  # MCV1 coverage uncertainty (95% CrI width)
  if ("pred_iqr" %in% names(grid_df) && !all(is.na(grid_df$pred_iqr))) {
    
    p_unc_ci <- ggplot() +
      geom_sf(data = grid_map_km, aes(fill = pred_iqr), colour = NA) +
      geom_sf(data = boundary_km, colour = "grey25", linewidth = 0.25) +
      scale_fill_gradient2(
        name = "95% CrI\nwidth",
        low = "white",
        high = "#781707",
        limits = c(0, NA)
      ) +
      labs(
        title    = paste(country_name, "- MCV1 prediction uncertainty (95% CrI width)"),
        subtitle = paste0(cell_km, " km grid; plogis scale; wider = less certain")
      ) +
      coord_sf(datum = NA) + map_theme
    
    uncertainty_plots[["mcv1_ci_width"]] <- p_unc_ci
  }
  
  # Covariate kriging uncertainty maps (where kriging variance is available)
  for (v in names(krige_variance)) {
    kv <- krige_variance[[v]]
    if (is.null(kv) || all(is.na(kv))) next
    
    grid_map_km$krige_sd <- sqrt(kv)
    
    info <- covariate_map_info(v, names(sub_m))
    
    p_cov_unc <- ggplot() +
      geom_sf(data = grid_map_km, aes(fill = krige_sd), colour = NA) +
      geom_sf(data = boundary_km, colour = "grey25", linewidth = 0.25) +
      scale_fill_gradient2(
        name = "Kriging\nSD",
        low = "white",
        high = "#781707",
        limits = c(0, NA)
      ) +
      labs(
        title    = paste(info$label, "- kriging uncertainty"),
        subtitle = "Higher values = fewer nearby clusters = less certain interpolation"
      ) +
      coord_sf(datum = NA) + map_theme
    
    uncertainty_plots[[paste0("cov_", v)]] <- p_cov_unc
    grid_map_km$krige_sd <- NULL  # clean up temp column
  }
  
  # ---- Covariate surface maps (same logic as original, but now using kriged values) ----
  covariate_plots <- purrr::map(map_infos, function(info) {
    d <- grid_map_km
    d$map_value <- d[[info$map_var]]
    
    is_prop <- identical(info$legend, "Proportion")
    lims <- get_plot_limits(d$map_value, info$map_var)
    if (any(!is.finite(lims)) || lims[1] == lims[2]) lims <- NULL
    
    subtitle <- paste(
      info$subtitle,
      limit_note(info$map_var),
      sep = ifelse(is.null(limit_note(info$map_var)), "", " ")
    )
    
    ggplot() +
      geom_sf(data = d, aes(fill = map_value), colour = NA) +
      geom_sf(data = boundary_km, colour = "grey25", linewidth = 0.25) +
      scale_fill_viridis_c(
        name   = info$legend,
        limits = if (is_prop && is.null(lims)) c(0, 1) else lims,
        labels = if (is_prop) scales::label_percent(accuracy = 1) else waiver(),
        oob    = scales::squish,
        na.value = "grey90"
      ) +
      labs(title = info$label, subtitle = subtitle) +
      coord_sf(datum = NA) + map_theme
  })
  names(covariate_plots) <- vars_final
  
  
  # ---- Return everything ----
  list(
    grid              = st_drop_geometry(grid_map_km),
    surface_plot      = p_pred / p_spatial,
    uncertainty_plots = uncertainty_plots,
    covariate_plots   = covariate_plots
  )
}










get_spatial_var <- function(res, spde, tag = "spatial") {
  if (is.null(res)) return(list(range = NA, var = NA))
  
  r <- INLA::inla.spde2.result(res, tag, spde, do.transf = TRUE)
  
  range_km <- INLA::inla.emarginal(
    function(x) x,
    r$marginals.range.nominal[[1]]
  )
  
  variance <- INLA::inla.emarginal(
    function(x) x,
    r$marginals.variance.nominal[[1]]
  )
  
  list(range = range_km, var = variance)
}


# ════════════════════════════════════════════════════════════════════════════════
# MODELS 5 & 6: ZERO-AND-ONE INFLATED BINOMIAL (ZOIB) GEOSTATISTICAL MODELS
# ════════════════════════════════════════════════════════════════════════════════
#
# PURPOSE
# ────────
# These models extend Models 3-4 by explicitly accounting for BOTH boundaries
# of the coverage distribution (0% and 100%).  The standard binomial and
# beta-binomial likelihoods treat all zeros and ones as arising from the same
# data-generating process.  But in DHS cluster data, zeros and ones can be:
#
#   (a) STRUCTURAL: a cluster truly has zero/complete vaccination because of
#       unmeasured local factors (e.g. programme failure, complete outreach)
#   (b) SAMPLING: a small cluster (e.g. n=3) happens to show 0/3 or 3/3
#       by chance, even if the true underlying rate is moderate
#
# The Zero-and-One Inflated Binomial (ZOIB) mixture model separates these:
#
#   P(y_i = 0) = π₀ + (1 - π₀ - π₁) × Binom(0 | n_i, p_i)
#   P(y_i = n_i) = π₁ + (1 - π₀ - π₁) × Binom(n_i | n_i, p_i)
#   P(y_i = k) = (1 - π₀ - π₁) × Binom(k | n_i, p_i)  for 0 < k < n_i
#
# where:
#   π₀ = probability of a "structural zero" (excess zero-coverage cluster)
#   π₁ = probability of a "structural one" (excess full-coverage cluster)
#   p_i = binomial success probability from covariates + spatial field
#
#
# TWO VARIANTS
# ─────────────
# Model 5: FIXED ZOIB
#   π₀ and π₁ are FIXED per country, estimated from the observed proportions
#   of 0% and 100% clusters relative to what a standard binomial would predict.
#   The spatial + covariate structure is the same as Model 3.
#
# Model 6: COVARIATE-DEPENDENT ZOIB
#   π₀ᵢ and π₁ᵢ VARY by cluster, modelled as functions of covariates:
#     logit(π₀ᵢ) = γ₀ + spatial_zero(sᵢ)
#     logit(π₁ᵢ) = δ₀ + spatial_one(sᵢ)
#   This uses INLA's multi-likelihood trick: three response columns
#   (is_zero, is_one, y_interior) with separate linear predictors sharing
#   the mesh but with independent spatial fields.
#
#
# IMPLEMENTATION STRATEGY
# ────────────────────────
# INLA does not have a built-in zero-AND-one inflated binomial family.
# It only supports zero-inflation (zeroinflatedbinomial0/1).
#
# MODEL 5 uses a data-augmentation approach:
#   - Estimate π₀_excess and π₁_excess from comparing observed vs expected
#     boundary proportions under a fitted binomial
#   - Fix the zero-inflation parameter in INLA's zeroinflatedbinomial1 to
#     the estimated π₀_excess
#   - Transform the data (flip y → n-y) and run a second ZIB to estimate
#     the one-inflation component
#   - OR (preferred): use inla.rgeneric for a custom three-component likelihood
#   - SIMPLEST (implemented here): use the "classification trick" where we
#     pre-classify boundary observations as structural vs sampling based on
#     cluster size, then fit a standard binomial on the "sampling" subset
#     with the structural boundaries contributing through adjusted weights.
#
#   ACTUALLY IMPLEMENTED: We use the most robust approach within INLA:
#   Fix π₀ and π₁ at country-level estimates, then fit the binomial INLA-SPDE
#   on the FULL data but with a custom likelihood via zeroinflatedbinomial1
#   for the zero-inflation part. For the one-inflation, we use the algebraic
#   equivalence: transform the problem by fitting TWO linked models.
#
#   FINAL APPROACH (Model 5): We estimate π₀ and π₁ from observed data,
#   construct the ZOIB log-likelihood contributions manually, and feed them
#   into INLA via the rgeneric mechanism.
#
#   PRACTICAL APPROACH (Model 5 — what we actually do):
#   Since rgeneric is complex and fragile, we use a SIMPLER but valid method:
#   (1) Estimate excess π₀ and π₁ from data
#   (2) Use INLA's zeroinflatedbinomial1 with π₀ FIXED at the estimated value
#   (3) Compute one-inflation-adjusted predictions post-hoc
#   This gives us zero-inflation within INLA + one-inflation as a post-hoc
#   mixture correction.
#
# MODEL 6 uses INLA's multi-likelihood / joint model approach:
#   Three response columns sharing effects through copy/replicate.
#   Component 1: Bernoulli — is cluster a structural zero?
#   Component 2: Bernoulli — is cluster a structural one? (given not struct. zero)
#   Component 3: Binomial  — coverage among "interior" clusters
#
#
# LITERATURE BASIS
# ─────────────────
# Diallo et al. (2012). PLOS ONE. Bayesian ZIB geostatistical modelling.
#   - Used ZIB for malaria with MCMC; single π₀ parameter; no one-inflation.
#   - Found ZIB had better predictive ability than standard binomial.
#
# Sweeney, Haslett & Parnell (2014/2025). Biometrika.
#   - Developed zero-&-N-inflated binomial distribution.
#   - Showed covariates in inflation components improve understanding.
#   - Proved the distribution arises naturally from constrained Poisson.
#
# Asmarian et al. (2019). Int J Environ Res Public Health.
#   - Bayesian spatial joint model for ZI data with R-INLA.
#   - Used two-likelihood trick: Bernoulli + truncated Poisson sharing
#     spatial effects via BYM2.
#
# Sugasawa et al. (2025). arXiv:2508.05041.
#   - Boundary-inflated binomial (BIB) with spatio-temporal components.
#   - Three-component mixture: Dirac(0) + Dirac(N) + Binomial.
#   - Used Pólya-Gamma augmentation + Gaussian predictive processes.
#
# Hall (2000). Biometrics.
#   - Original zero-inflated binomial formulation.
#
# Lambert (1992). Technometrics.
#   - Zero-inflated Poisson (foundational for all ZI models).
# ════════════════════════════════════════════════════════════════════════════════


# ════════════════════════════════════════════════════════════════════════════════
# HELPER: Estimate excess boundary probabilities
# ════════════════════════════════════════════════════════════════════════════════
#
# Under a standard Binomial(n_i, p_i), the probability of observing y_i = 0
# is (1 - p_i)^n_i and the probability of y_i = n_i is p_i^n_i.
# 
# If the OBSERVED proportion of zeros exceeds what the fitted binomial predicts,
# the excess is attributed to "structural zeros".
#
# We estimate p_i from the GLM (Model 1) fitted values, then:
#   E[fraction of zeros under binomial] = mean( (1-p_hat_i)^n_i )
#   observed fraction of zeros = mean( y_i == 0 )
#   π₀_excess = max(0, observed - expected)
#
# Similarly for ones.

# ════════════════════════════════════════════════════════════════════════════════
# MODELS 5 & 6: ZERO-INFLATED BINOMIAL MODELS
# Fully comparable to binomial / beta-binomial models via WAIC, DIC, R2, MAE, RMSE
# ════════════════════════════════════════════════════════════════════════════════


# ════════════════════════════════════════════════════════════════════════════════
# HELPER: Estimate excess zero probability
# ════════════════════════════════════════════════════════════════════════════════

estimate_excess_zero_prob <- function(sub_m, pred_col = "pred_glm") {
  
  n <- sub_m$n_children
  y <- sub_m$n_vaccinated
  
  if (pred_col %in% names(sub_m) && !all(is.na(sub_m[[pred_col]]))) {
    p_hat <- sub_m[[pred_col]]
  } else {
    p_hat <- rep(sum(y) / sum(n), length(n))
  }
  
  p_hat <- pmax(pmin(p_hat, 1 - 1e-10), 1e-10)
  
  expected_zero <- mean((1 - p_hat)^n)
  observed_zero <- mean(y == 0)
  
  pi0_excess <- max(0, observed_zero - expected_zero)
  pi0_excess <- min(pi0_excess, 0.95)
  
  list(
    pi0_excess    = pi0_excess,
    pi1_excess    = 0,
    observed_zero = observed_zero,
    expected_zero = expected_zero,
    observed_one  = mean(y == n),
    expected_one  = mean(p_hat^n),
    n_clusters    = length(n)
  )
}


# Backwards-compatible alias if your main script already calls this name
estimate_excess_boundary_probs <- estimate_excess_zero_prob


# ════════════════════════════════════════════════════════════════════════════════
# HELPER: Comparable prediction metrics
# ════════════════════════════════════════════════════════════════════════════════

compute_comparable_metrics <- function(y, n, pred, inla_result = NULL) {
  
  obs <- y / n
  pred <- pmax(pmin(pred, 1), 0)
  
  keep <- is.finite(obs) & is.finite(pred)
  obs <- obs[keep]
  pred <- pred[keep]
  
  mae  <- mean(abs(obs - pred))
  rmse <- sqrt(mean((obs - pred)^2))
  
  sse <- sum((obs - pred)^2)
  sst <- sum((obs - mean(obs))^2)
  r2 <- ifelse(sst > 0, 1 - sse / sst, NA_real_)
  
  out <- list(
    n_obs = length(obs),
    mae = mae,
    rmse = rmse,
    r2 = r2
  )
  
  if (!is.null(inla_result)) {
    out$dic <- if (!is.null(inla_result$dic$dic)) inla_result$dic$dic else NA_real_
    out$waic <- if (!is.null(inla_result$waic$waic)) inla_result$waic$waic else NA_real_
    
    if (!is.null(inla_result$cpo$cpo)) {
      cpo <- inla_result$cpo$cpo
      cpo <- cpo[is.finite(cpo) & cpo > 0]
      out$mean_log_cpo <- ifelse(length(cpo) > 0, mean(log(cpo)), NA_real_)
      out$neg_log_cpo  <- ifelse(length(cpo) > 0, -sum(log(cpo)), NA_real_)
    } else {
      out$mean_log_cpo <- NA_real_
      out$neg_log_cpo  <- NA_real_
    }
  }
  
  out
}


# ════════════════════════════════════════════════════════════════════════════════
# HELPER: Extract zero-inflation probability from INLA object
# ════════════════════════════════════════════════════════════════════════════════

extract_zib_pi0 <- function(inla_result) {
  
  if (is.null(inla_result$summary.hyperpar)) return(NA_real_)
  
  hp <- inla_result$summary.hyperpar
  rn <- rownames(hp)
  
  theta_row <- grep("zero|inflation|prob|theta", rn, ignore.case = TRUE)
  
  if (length(theta_row) == 0) return(NA_real_)
  
  theta_mean <- hp$mean[theta_row[1]]
  
  # INLA reports the internal parameter on logit scale for this family.
  pi0 <- plogis(theta_mean)
  pmax(pmin(pi0, 1), 0)
}


# ════════════════════════════════════════════════════════════════════════════════
# HELPER: Build robust model formula
# ════════════════════════════════════════════════════════════════════════════════

build_zib_formula <- function(cov_terms, iid_spec, spde_name = "spde") {
  
  rhs_terms <- c("intercept")
  
  if (!is.null(cov_terms) && nzchar(cov_terms)) {
    rhs_terms <- c(rhs_terms, cov_terms)
  }
  
  rhs_terms <- c(rhs_terms, paste0("f(spatial, model = ", spde_name, ")"))
  
  if (!is.null(iid_spec) && nzchar(iid_spec)) {
    rhs_terms <- c(rhs_terms, iid_spec)
  }
  
  as.formula(paste("y ~ -1 +", paste(rhs_terms, collapse = " + ")))
}


# ═══════════════════════════════════════════════════════════════════════
# ROBUST HELPER: Get INLA stack data and Ntrials safely
# ═══════════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════════
# ROBUST HELPER: Get INLA stack data and Ntrials safely
# ═══════════════════════════════════════════════════════════════════════

get_inla_stack_data_and_ntrials <- function(stk_template, spde, sub_m) {
  
  stk_data_spde <- INLA::inla.stack.data(stk_template, spde = spde)
  stk_data_plain <- INLA::inla.stack.data(stk_template)
  
  # INLA stack data is a list, so use the response length, not nrow().
  n_stack <- NULL
  
  if ("y" %in% names(stk_data_spde)) {
    n_stack <- length(stk_data_spde$y)
  } else if ("Y" %in% names(stk_data_spde)) {
    n_stack <- nrow(stk_data_spde$Y)
    if (is.null(n_stack)) n_stack <- length(stk_data_spde$Y)
  } else if ("y" %in% names(stk_data_plain)) {
    n_stack <- length(stk_data_plain$y)
  } else if ("Y" %in% names(stk_data_plain)) {
    n_stack <- nrow(stk_data_plain$Y)
    if (is.null(n_stack)) n_stack <- length(stk_data_plain$Y)
  } else {
    n_stack <- nrow(sub_m)
  }
  
  n_trials <- NULL
  
  if ("n_trial" %in% names(stk_data_plain)) {
    n_trials <- stk_data_plain$n_trial
  } else if ("n_trial" %in% names(stk_data_spde)) {
    n_trials <- stk_data_spde$n_trial
  } else if ("n_trials" %in% names(stk_data_plain)) {
    n_trials <- stk_data_plain$n_trials
  } else if ("n_trials" %in% names(stk_data_spde)) {
    n_trials <- stk_data_spde$n_trials
  } else if ("Ntrials" %in% names(stk_data_plain)) {
    n_trials <- stk_data_plain$Ntrials
  } else if ("Ntrials" %in% names(stk_data_spde)) {
    n_trials <- stk_data_spde$Ntrials
  } else {
    n_trials <- sub_m$n_children
  }
  
  n_trials <- as.integer(n_trials)
  
  if (length(n_trials) == 0) {
    n_trials <- as.integer(sub_m$n_children)
  }
  
  if (length(n_trials) != n_stack) {
    if (length(sub_m$n_children) == n_stack) {
      n_trials <- as.integer(sub_m$n_children)
    } else {
      stop(
        paste0(
          "Ntrials length does not match INLA stack rows. ",
          "length(Ntrials)=", length(n_trials),
          "; stack rows=", n_stack,
          "; nrow(sub_m)=", nrow(sub_m)
        )
      )
    }
  }
  
  if (any(is.na(n_trials)) || any(n_trials <= 0)) {
    stop("Ntrials contains NA or non-positive values.")
  }
  
  list(
    data = stk_data_spde,
    ntrials = n_trials,
    n_stack = n_stack
  )
}


# ═══════════════════════════════════════════════════════════════════════
# MODEL 5: FIXED ZERO-INFLATED BINOMIAL
# ═══════════════════════════════════════════════════════════════════════

fit_model5_zib_fixed <- function(sub_m, vars_final, mesh, spde, A, stk_template,
                                 iid_spec, cov_terms, boundary_probs,
                                 run_inla_fn = NULL) {
  
  pi0 <- boundary_probs$pi0_excess
  pi0 <- pmax(pmin(pi0, 0.95), 0)
  
  cat("\n  Model 5: Zero-Inflated Binomial (FIXED inflation)\n")
  cat("    Estimated structural zero probability (pi0):", round(pi0, 4), "\n")
  cat("    Binomial component weight (1 - pi0):        ", round(1 - pi0, 4), "\n")
  
  theta_zero <- if (pi0 < 1e-8) {
    -20
  } else if (pi0 > 1 - 1e-8) {
    20
  } else {
    qlogis(pi0)
  }
  
  cat("    Fixed zero-inflation logit(theta):", round(theta_zero, 3), "\n")
  
  f5 <- as.formula(paste0(
    "y ~ -1 + intercept + ", cov_terms,
    " + f(spatial, model = spde) + ", iid_spec
  ))
  
  sd_nt <- get_inla_stack_data_and_ntrials(stk_template, spde, sub_m)
  stk_data <- sd_nt$data
  n_trials <- sd_nt$ntrials
  
  cat("    Ntrials length:", length(n_trials),
      "| stack rows:", sd_nt$n_stack,
      "| min/max:", min(n_trials), "/", max(n_trials), "\n")
  
  inla_m5 <- tryCatch(
    INLA::inla(
      formula = f5,
      family  = "zeroinflatedbinomial1",
      Ntrials = n_trials,
      data    = stk_data,
      control.predictor = list(
        A = INLA::inla.stack.A(stk_template),
        compute = TRUE,
        link = 1
      ),
      control.compute = list(
        dic = TRUE,
        waic = TRUE,
        cpo = TRUE,
        config = TRUE
      ),
      control.inla = list(strategy = "adaptive"),
      control.family = list(
        hyper = list(
          theta = list(
            initial = theta_zero,
            fixed = TRUE
          )
        )
      ),
      verbose = FALSE
    ),
    error = function(e) {
      cat("    Model 5 ZIB failed:", e$message, "\n")
      NULL
    }
  )
  
  if (is.null(inla_m5)) return(NULL)
  
  idx <- INLA::inla.stack.index(stk_template, "estimation")$data
  
  p_zib <- inla_m5$summary.fitted.values$mean[idx]
  p_zib <- pmax(pmin(p_zib, 1), 0)
  
  p_binom_component <- if (pi0 < 1 - 1e-8) p_zib / (1 - pi0) else p_zib
  p_binom_component <- pmax(pmin(p_binom_component, 1), 0)
  
  metrics <- calc_pred_metrics(sub_m$coverage, p_zib)
  
  cat("    DIC:", round(inla_m5$dic$dic, 1),
      "| WAIC:", round(inla_m5$waic$waic, 1),
      "| R2:", metrics$r2,
      "| MAE:", metrics$mae, "\n")
  
  feff_m5 <- inla_m5$summary.fixed %>%
    mutate(
      term = rownames(.),
      sig  = ifelse(`0.025quant` > 0 | `0.975quant` < 0, "*", ""),
      OR    = exp(mean),
      OR_lo = exp(`0.025quant`),
      OR_hi = exp(`0.975quant`)
    ) %>%
    dplyr::select(term, mean, sd, `0.025quant`, `0.975quant`,
                  OR, OR_lo, OR_hi, sig)
  
  cat("\n  Fixed effects (Model 5 — binomial component):\n")
  print(feff_m5 %>% mutate(across(where(is.numeric), ~ round(.x, 3))),
        row.names = FALSE)
  
  sp5 <- get_spatial_var(inla_m5, spde)
  
  cat(sprintf("\n  Spatial range: %.1f km | Spatial variance: %.4f\n",
              sp5$range, sp5$var))
  
  list(
    inla_result     = inla_m5,
    predictions     = p_zib,
    pred_zib_raw    = p_zib,
    pred_binom_comp = p_binom_component,
    fixed_effects   = feff_m5,
    spatial_range   = sp5$range,
    spatial_var     = sp5$var,
    pi0             = pi0,
    pi1             = 0,
    dic             = inla_m5$dic$dic,
    waic            = inla_m5$waic$waic,
    r2              = metrics$r2,
    mae             = metrics$mae,
    rmse            = metrics$rmse,
    metrics         = metrics
  )
}


# ═══════════════════════════════════════════════════════════════════════
# MODEL 6: ESTIMATED ZERO-INFLATED BINOMIAL
# ═══════════════════════════════════════════════════════════════════════

fit_model6_zib_estimated <- function(sub_m, vars_final, mesh, spde, A, stk_template,
                                     iid_spec, cov_terms, boundary_probs = NULL,
                                     run_inla_fn = NULL) {
  
  cat("\n  Model 6: Zero-Inflated Binomial (ESTIMATED inflation)\n")
  cat("    Estimating zero-inflation probability inside INLA\n")
  
  if (!is.null(boundary_probs) && !is.null(boundary_probs$pi0_excess)) {
    pi0_init <- pmax(pmin(boundary_probs$pi0_excess, 0.95), 1e-6)
  } else {
    pi0_init <- mean(sub_m$n_vaccinated == 0)
    pi0_init <- pmax(pmin(pi0_init, 0.95), 1e-6)
  }
  
  theta_init <- qlogis(pi0_init)
  
  cat("    Initial pi0:", round(pi0_init, 4),
      "| initial logit(theta):", round(theta_init, 3), "\n")
  
  f6 <- as.formula(paste0(
    "y ~ -1 + intercept + ", cov_terms,
    " + f(spatial, model = spde) + ", iid_spec
  ))
  
  sd_nt <- get_inla_stack_data_and_ntrials(stk_template, spde, sub_m)
  stk_data <- sd_nt$data
  n_trials <- sd_nt$ntrials
  
  cat("    Ntrials length:", length(n_trials),
      "| stack rows:", sd_nt$n_stack,
      "| min/max:", min(n_trials), "/", max(n_trials), "\n")
  
  inla_m6 <- tryCatch(
    INLA::inla(
      formula = f6,
      family  = "zeroinflatedbinomial1",
      Ntrials = n_trials,
      data    = stk_data,
      control.predictor = list(
        A = INLA::inla.stack.A(stk_template),
        compute = TRUE,
        link = 1
      ),
      control.compute = list(
        dic = TRUE,
        waic = TRUE,
        cpo = TRUE,
        config = TRUE
      ),
      control.inla = list(strategy = "adaptive"),
      control.family = list(
        hyper = list(
          theta = list(
            initial = theta_init,
            fixed = FALSE
          )
        )
      ),
      verbose = FALSE
    ),
    error = function(e) {
      cat("    Model 6 ZIB failed:", e$message, "\n")
      NULL
    }
  )
  
  if (is.null(inla_m6)) return(NULL)
  
  idx <- INLA::inla.stack.index(stk_template, "estimation")$data
  
  p_zib <- inla_m6$summary.fitted.values$mean[idx]
  p_zib <- pmax(pmin(p_zib, 1), 0)
  
  pi0_hat <- extract_zib_pi0(inla_m6)
  if (is.na(pi0_hat)) pi0_hat <- pi0_init
  
  p_binom_component <- if (pi0_hat < 1 - 1e-8) p_zib / (1 - pi0_hat) else p_zib
  p_binom_component <- pmax(pmin(p_binom_component, 1), 0)
  
  metrics <- calc_pred_metrics(sub_m$coverage, p_zib)
  
  cat("    Estimated pi0:", round(pi0_hat, 4), "\n")
  cat("    DIC:", round(inla_m6$dic$dic, 1),
      "| WAIC:", round(inla_m6$waic$waic, 1),
      "| R2:", metrics$r2,
      "| MAE:", metrics$mae, "\n")
  
  feff_m6 <- inla_m6$summary.fixed %>%
    mutate(
      term = rownames(.),
      sig  = ifelse(`0.025quant` > 0 | `0.975quant` < 0, "*", ""),
      OR    = exp(mean),
      OR_lo = exp(`0.025quant`),
      OR_hi = exp(`0.975quant`)
    ) %>%
    dplyr::select(term, mean, sd, `0.025quant`, `0.975quant`,
                  OR, OR_lo, OR_hi, sig)
  
  cat("\n  Fixed effects (Model 6 — binomial component):\n")
  print(feff_m6 %>% mutate(across(where(is.numeric), ~ round(.x, 3))),
        row.names = FALSE)
  
  sp6 <- get_spatial_var(inla_m6, spde)
  
  cat(sprintf("\n  Spatial range: %.1f km | Spatial variance: %.4f\n",
              sp6$range, sp6$var))
  
  list(
    inla_result      = inla_m6,
    predictions      = p_zib,
    pred_zib_raw     = p_zib,
    pred_binom_comp  = p_binom_component,
    pred_pi0         = rep(pi0_hat, length(p_zib)),
    pred_pi1         = rep(0, length(p_zib)),
    pred_binom       = p_binom_component,
    fixed_effects    = feff_m6,
    spatial_range    = sp6$range,
    spatial_var      = sp6$var,
    pi0              = pi0_hat,
    pi1              = 0,
    dic              = inla_m6$dic$dic,
    waic             = inla_m6$waic$waic,
    r2               = metrics$r2,
    mae              = metrics$mae,
    rmse             = metrics$rmse,
    metrics          = metrics,
    method           = "estimated_global_zero_inflation"
  )
}
# ═════════════════════════════════════════════════════════════════════
# ═══  HELPER FUNCTION: Prediction metrics                         ═══
# ═════════════════════════════════════════════════════════════════════
#
# Computes prediction metrics for observed vs predicted coverage:
#   r    = Pearson correlation
#   R²   = proportion of variance explained
#   MAE  = mean absolute error
#   RMSE = root mean squared error

calc_pred_metrics <- function(observed, predicted) {
  ok <- !is.na(observed) & !is.na(predicted)
  
  if (sum(ok) < 10) {
    return(list(r = NA, r2 = NA, mae = NA, rmse = NA, n = sum(ok)))
  }
  
  obs  <- observed[ok]
  pred <- predicted[ok]
  
  r <- cor(obs, pred)
  
  ss_res <- sum((obs - pred)^2)
  ss_tot <- sum((obs - mean(obs))^2)
  
  r2 <- ifelse(ss_tot > 0, 1 - ss_res / ss_tot, NA)
  mae <- mean(abs(obs - pred))
  rmse <- sqrt(mean((obs - pred)^2))
  
  list(
    r = round(r, 4),
    r2 = round(r2, 4),
    mae = round(mae, 4),
    rmse = round(rmse, 4),
    n = sum(ok)
  )
}


make_cv_folds <- function(sub_m, k = 3,
                          scheme = c("random", "spatial"),
                          seed = 123) {
  scheme <- match.arg(scheme)
  set.seed(seed)
  
  n <- nrow(sub_m)
  
  if (scheme == "random") {
    return(sample(rep(seq_len(k), length.out = n)))
  }
  
  coords <- cbind(sub_m$LONGNUM, sub_m$LATNUM)
  coords_scaled <- scale(coords)
  
  km <- stats::kmeans(coords_scaled, centers = k, nstart = 50)
  km$cluster
}


crps_empirical <- function(samples, observed) {
  samples <- samples[is.finite(samples)]
  
  if (length(samples) < 2 || !is.finite(observed)) {
    return(NA_real_)
  }
  
  term1 <- mean(abs(samples - observed))
  
  xs <- sort(samples)
  s <- length(xs)
  term2 <- sum((2 * seq_len(s) - s - 1) * xs) / (s^2)
  
  term1 - term2
}


sample_fitted_values <- function(inla_fit, idx, n_samples = 1000) {
  out <- matrix(NA_real_, nrow = length(idx), ncol = n_samples)
  
  if (is.null(inla_fit$marginals.fitted.values)) {
    return(out)
  }
  
  for (i in seq_along(idx)) {
    marg <- inla_fit$marginals.fitted.values[[idx[i]]]
    
    if (!is.null(marg)) {
      out[i, ] <- INLA::inla.rmarginal(n_samples, marg)
    }
  }
  
  out
}


extract_cv_zib_pi0 <- function(inla_fit, fallback = NA_real_) {
  
  if (is.null(inla_fit) || is.null(inla_fit$summary.hyperpar)) {
    return(fallback)
  }
  
  hp <- inla_fit$summary.hyperpar
  hp_names <- rownames(hp)
  
  theta_idx <- grep("zero|inflation|theta|prob", hp_names, ignore.case = TRUE)
  
  if (length(theta_idx) == 0) {
    return(fallback)
  }
  
  theta_mean <- hp$mean[theta_idx[1]]
  pi0 <- stats::plogis(theta_mean)
  
  pmax(pmin(pi0, 0.999), 0)
}


estimate_cv_excess_zero_prob <- function(train_m, pred_col = NULL) {
  
  n <- train_m$n_children
  y <- train_m$n_vaccinated
  
  if (!is.null(pred_col) &&
      pred_col %in% names(train_m) &&
      !all(is.na(train_m[[pred_col]]))) {
    p_hat <- train_m[[pred_col]]
  } else {
    p_hat <- rep(sum(y) / sum(n), length(n))
  }
  
  p_hat <- pmax(pmin(p_hat, 1 - 1e-10), 1e-10)
  
  expected_zero <- mean((1 - p_hat)^n)
  observed_zero <- mean(y == 0)
  
  pi0_excess <- max(0, observed_zero - expected_zero)
  pi0_excess <- min(pi0_excess, 0.95)
  
  list(
    pi0_excess = pi0_excess,
    pi1_excess = 0,
    observed_zero = observed_zero,
    expected_zero = expected_zero,
    observed_one = mean(y == n),
    expected_one = mean(p_hat^n),
    n_clusters = length(n)
  )
}


summarise_cv_predictions <- function(cv_pred_df) {
  cv_pred_df %>%
    dplyr::group_by(scheme, model) %>%
    dplyr::summarise(
      avg_bias = mean(pred - obs, na.rm = TRUE),
      mae = mean(abs(pred - obs), na.rm = TRUE),
      rmse = sqrt(mean((obs - pred)^2, na.rm = TRUE)),
      r = cor(obs, pred, use = "complete.obs"),
      r2 = 1 - sum((obs - pred)^2, na.rm = TRUE) /
        sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE),
      crps = mean(crps, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(scheme, mae)
}


run_cv_models_2_to_6_utazi_style <- function(sub_m, vars_final,
                                             k = 3,
                                             schemes = c("random", "spatial"),
                                             seed = 123,
                                             boundary_probs = NULL,
                                             n_crps_samples = 1000) {
  
  all_preds <- list()
  
  for (scheme in schemes) {
    cat("\n---", toupper(scheme), k, "-fold cross-validation ---\n")
    
    sub_m$cv_fold <- make_cv_folds(
      sub_m = sub_m,
      k = k,
      scheme = scheme,
      seed = seed
    )
    
    for (fold in seq_len(k)) {
      cat("  Fold", fold, "of", k, "\n")
      
      train_m <- sub_m[sub_m$cv_fold != fold, , drop = FALSE]
      test_m  <- sub_m[sub_m$cv_fold == fold, , drop = FALSE]
      
      coords_train <- cbind(train_m$LONGNUM, train_m$LATNUM)
      coords_test  <- cbind(test_m$LONGNUM, test_m$LATNUM)
      
      cent <- colMeans(coords_train)
      
      lcrs <- paste0(
        "+proj=aeqd +lat_0=", cent[2],
        " +lon_0=", cent[1],
        " +x_0=0 +y_0=0 +datum=WGS84 +units=km"
      )
      
      cp_train <- tryCatch(
        sf::sf_project(sf::st_crs(4326)$proj4string, lcrs, coords_train),
        error = function(e) NULL
      )
      
      cp_test <- tryCatch(
        sf::sf_project(sf::st_crs(4326)$proj4string, lcrs, coords_test),
        error = function(e) NULL
      )
      
      if (is.null(cp_train) || is.null(cp_test)) next
      
      ext <- max(diff(range(cp_train[, 1])), diff(range(cp_train[, 2])))
      rg <- max(ext / 3, 20)
      
      mesh <- tryCatch(
        INLA::inla.mesh.2d(
          loc = cp_train,
          max.edge = c(rg / 5, rg),
          cutoff = max(rg / 20, 5),
          offset = c(rg / 2, rg * 1.5)
        ),
        error = function(e) NULL
      )
      
      if (is.null(mesh)) next
      
      spde <- INLA::inla.spde2.pcmatern(
        mesh = mesh,
        alpha = 2,
        prior.range = c(rg, 0.5),
        prior.sigma = c(1, 0.01)
      )
      
      A_train <- INLA::inla.spde.make.A(mesh = mesh, loc = cp_train)
      A_test  <- INLA::inla.spde.make.A(mesh = mesh, loc = cp_test)
      
      cov_train <- train_m %>% dplyr::select(dplyr::all_of(vars_final))
      cov_test  <- test_m %>% dplyr::select(dplyr::all_of(vars_final))
      
      cov_terms <- paste(vars_final, collapse = " + ")
      
      iid_spec <- paste0(
        "f(cluster_id, model = 'iid', ",
        "hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))"
      )
      
      stk_train <- INLA::inla.stack(
        data = list(
          y = train_m$n_vaccinated,
          n_trial = train_m$n_children
        ),
        A = list(A_train, 1),
        effects = list(
          spatial = 1:spde$n.spde,
          data.frame(
            intercept = 1,
            cluster_id = seq_len(nrow(train_m)),
            cov_train
          )
        ),
        tag = "train"
      )
      
      stk_test <- INLA::inla.stack(
        data = list(
          y = NA,
          n_trial = test_m$n_children
        ),
        A = list(A_test, 1),
        effects = list(
          spatial = 1:spde$n.spde,
          data.frame(
            intercept = 1,
            cluster_id = NA_integer_,
            cov_test
          )
        ),
        tag = "test"
      )
      
      stk_cv <- INLA::inla.stack(stk_train, stk_test)
      stk_data <- INLA::inla.stack.data(stk_cv, spde = spde)
      
      idx_test <- INLA::inla.stack.index(stk_cv, "test")$data
      
      obs_cluster <- test_m$n_vaccinated / test_m$n_children
      
      f2 <- as.formula(paste0(
        "y ~ -1 + intercept + f(spatial, model = spde) + ", iid_spec
      ))
      
      f3 <- as.formula(paste0(
        "y ~ -1 + intercept + ", cov_terms,
        " + f(spatial, model = spde) + ", iid_spec
      ))
      
      run_cv_inla <- function(formula, family, control.family = NULL) {
        tryCatch(
          INLA::inla(
            formula = formula,
            family = family,
            Ntrials = stk_data$n_trial,
            data = stk_data,
            control.predictor = list(
              A = INLA::inla.stack.A(stk_cv),
              compute = TRUE,
              link = 1
            ),
            control.compute = list(
              dic = TRUE,
              waic = TRUE,
              cpo = TRUE,
              config = TRUE
            ),
            control.inla = list(strategy = "adaptive"),
            control.family = control.family,
            verbose = FALSE
          ),
          error = function(e) {
            cat("    CV INLA failed:", e$message, "\n")
            NULL
          }
        )
      }
      
      # Model 5 fixed pi0:
      # Use supplied boundary_probs if provided; otherwise estimate within
      # each training fold only, to avoid leakage from held-out data.
      fold_boundary_probs <- if (!is.null(boundary_probs)) {
        boundary_probs
      } else {
        estimate_cv_excess_zero_prob(train_m)
      }
      
      pi0_fixed <- fold_boundary_probs$pi0_excess
      pi0_fixed <- pmax(pmin(pi0_fixed, 0.95), 0)
      
      theta_fixed <- if (pi0_fixed < 1e-8) {
        -20
      } else if (pi0_fixed > 1 - 1e-8) {
        20
      } else {
        qlogis(pi0_fixed)
      }
      
      # Model 6 estimated pi0:
      # Start from the training-fold excess-zero estimate, then estimate pi0.
      pi0_init <- pmax(pmin(pi0_fixed, 0.95), 1e-6)
      theta_init <- qlogis(pi0_init)
      
      fit_list <- list(
        `Model 2: Spatial-only` = run_cv_inla(f2, "binomial"),
        `Model 3: Full INLA` = run_cv_inla(f3, "binomial"),
        `Model 4: Beta-binomial` = run_cv_inla(f3, "betabinomial"),
        `Model 5: ZIB-Fixed` = run_cv_inla(
          f3,
          "zeroinflatedbinomial1",
          control.family = list(
            hyper = list(
              theta = list(initial = theta_fixed, fixed = TRUE)
            )
          )
        ),
        `Model 6: ZIB-Estimated` = run_cv_inla(
          f3,
          "zeroinflatedbinomial1",
          control.family = list(
            hyper = list(
              theta = list(initial = theta_init, fixed = FALSE)
            )
          )
        )
      )
      
      for (model_name in names(fit_list)) {
        fit <- fit_list[[model_name]]
        if (is.null(fit)) next
        
        pred <- fit$summary.fitted.values$mean[idx_test]
        pred <- pmax(pmin(pred, 1), 0)
        
        fitted_samples <- sample_fitted_values(
          inla_fit = fit,
          idx = idx_test,
          n_samples = n_crps_samples
        )
        fitted_samples <- pmax(pmin(fitted_samples, 1), 0)
        
        model_pi0 <- NA_real_
        
        if (model_name == "Model 5: ZIB-Fixed") {
          model_pi0 <- pi0_fixed
        }
        
        if (model_name == "Model 6: ZIB-Estimated") {
          model_pi0 <- extract_cv_zib_pi0(fit, fallback = pi0_init)
        }
        
        crps <- vapply(
          seq_along(obs_cluster),
          function(i) crps_empirical(fitted_samples[i, ], obs_cluster[i]),
          numeric(1)
        )
        
        all_preds[[length(all_preds) + 1]] <- data.frame(
          scheme = scheme,
          fold = fold,
          model = model_name,
          obs = obs_cluster,
          pred = pred,
          crps = crps,
          pi0 = model_pi0,
          n_children = test_m$n_children,
          n_vaccinated = test_m$n_vaccinated
        )
      }
    }
  }
  
  cv_pred_df <- dplyr::bind_rows(all_preds)
  
  cv_summary <- summarise_cv_predictions(cv_pred_df) %>%
    dplyr::left_join(
      cv_pred_df %>%
        dplyr::group_by(scheme, model) %>%
        dplyr::summarise(
          mean_pi0 = mean(pi0, na.rm = TRUE),
          sd_pi0 = sd(pi0, na.rm = TRUE),
          .groups = "drop"
        ),
      by = c("scheme", "model")
    )
  
  cv_plot <- ggplot2::ggplot(
    cv_summary,
    ggplot2::aes(x = reorder(model, mae), y = mae, fill = scheme)
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Out-of-sample predictive performance",
      subtitle = paste0(k, "-fold random and spatial cross-validation"),
      x = NULL,
      y = "Mean absolute error"
    ) +
    ggplot2::theme_minimal(base_size = 11)
  
  cv_r2_plot <- ggplot2::ggplot(
    cv_summary,
    ggplot2::aes(x = reorder(model, r2), y = r2, fill = scheme)
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Out-of-sample variance explained",
      subtitle = paste0(k, "-fold random and spatial cross-validation"),
      x = NULL,
      y = "Cross-validated R2"
    ) +
    ggplot2::theme_minimal(base_size = 11)
  
  list(
    predictions = cv_pred_df,
    summary = cv_summary,
    plot = cv_plot,
    r2_plot = cv_r2_plot
  )
}


# Backwards-compatible alias if your main script already calls the old function name
run_cv_models_2_to_5_utazi_style <- run_cv_models_2_to_6_utazi_style
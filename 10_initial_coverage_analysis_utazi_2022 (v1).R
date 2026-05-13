# ════════════════════════════════════════════════════════════════════════════════
# 10_initial_coverage_analysis.R
# ════════════════════════════════════════════════════════════════════════════════
#
# PURPOSE
# ────────
# This script analyses MCV1 (measles first-dose) vaccination coverage across
# sub-Saharan African countries using Demographic and Health Survey (DHS)
# cluster-level data.  It follows the methodology of Utazi et al. (2018, 2020,
# 2022) to:
#
#   1. Estimate OBSERVED coverage at each survey cluster (village/neighbourhood)
#   2. Predict coverage using a geostatistical model that combines socioeconomic
#      and demographic risk factors with spatial location
#   3. Identify WHERE pockets of low coverage remain after accounting for known
#      risk factors (residual spatial variation)
#   4. Identify WHICH population sub-groups are most at risk (covariate effects)
#
# The analysis produces country-level results with:
#   - Maps of observed vs predicted coverage and residuals
#   - Summary statistics on factor proportions (e.g. % with primary education,
#     % unemployed, % in urban slums) to describe the study population
#   - Odds ratios showing how each risk factor relates to vaccination
#   - Model comparison (non-spatial vs spatial, binomial vs beta-binomial)
#
#
# PIPELINE OVERVIEW (runs once per country)
# ──────────────────────────────────────────
#   A. LOAD data & covariates
#   B. PICK the GC (Geospatial Covariates) year closest to the survey
#   C. FILTER for age-eligible children (12–23 months with MCV1 data & GPS)
#   D. AGGREGATE individual data to cluster level:
#        - coverage = n_vaccinated / n_children
#        - n_children per cluster
#        - means of continuous covariates (+ standardised z-scores)
#        - proportions of each level for factor/categorical covariates
#          (e.g. % with primary education, % unemployed, etc.)
#   E. CALCULATE coverage-by-cluster-size table (shows the 0/1 boundary problem)
#   F. GENERATE RESULTS per country:
#        Step 1 — missingness screening (<5% rule)
#        Step 2 — bivariate screening (crude ORs, p < 0.2)
#        Step 3 — multicollinearity check (|r| > 0.8 + GVIF)
#        Step 4 — Model 1: covariates-only GLM (non-spatial baseline)
#        Step 5 — Model 2: spatial-only INLA-SPDE
#        Step 6 — Model 3: full INLA-SPDE (covariates + spatial) ← PRIMARY
#        Step 7 — Model 4: beta-binomial INLA-SPDE (overdispersion check)
#        Step 8 — model comparison (DIC, WAIC, AUC, R², MAE)
#   G. GENERATE MAPS with ADM1 boundaries:
#        - Observed MCV1 coverage
#        - Model-predicted coverage
#        - Residuals (observed − predicted)
#        - Obs vs pred scatter, coverage histogram, residual histogram
#   H. SAVE outputs in PDF (plots) and Excel (tables) form
#   I. RUN the loop for each country with ≥ 50 clusters
#
#
# WHY CLUSTER-LEVEL?
# ──────────────────
# DHS surveys sample clusters of ~25 households.  Each cluster has GPS
# coordinates (jittered ±2 km urban, ±5 km rural for confidentiality).
# The cluster is the natural unit for geostatistical modelling because:
#   (a) GPS coordinates attach to the cluster, not the individual child
#   (b) Children in the same cluster share unmeasured local factors
#       (nearby health facility quality, local leaders, road access)
# This is standard in the field — Utazi et al. (2018), Dong & Wakefield
# (2021), and Giorgi et al. (2018) all model cluster-level proportions.
#
#
# THE 0/1 BOUNDARY PROBLEM
# ─────────────────────────
# Many DHS clusters are small (5–15 children).  With so few children,
# coverage can only take a limited set of values (e.g. 0/5, 1/5, ..., 5/5).
# This means lots of clusters show exactly 0% or 100% coverage — not because
# they truly have zero or perfect vaccination, but because the sample is too
# small to detect intermediate values.  This is a SAMPLING ARTEFACT.
#
# We address this via four mechanisms:
#   (a) Binomial likelihood: handles 0/n and n/n naturally (unlike logit-OLS
#       which produces −Inf and +Inf at the boundaries)
#   (b) Empirical logit: Haldane-Anscombe correction log((y+0.5)/(n−y+0.5))
#       for any exploratory logit-based analysis
#   (c) iid cluster nugget in INLA: absorbs small-sample overdispersion
#   (d) Beta-binomial model: explicitly tests for extra-binomial variation
#
#
# KEY REFERENCES
# ──────────────
# Utazi CE et al. (2022). PLOS Global Public Health.
# Utazi CE et al. (2020). Lancet Digit Health 2: e536-e544.
# Utazi CE et al. (2018). Vaccine 36: 1583-1591.
# Dong TQ, Wakefield J (2021). Vaccine 39: 2557-2569.
# Mosser JF et al. (2019). Lancet 393: 1843-1855.
# Tessema GA et al. (2024). Geospatial inequalities in zero-dose.
# Acharya P et al. (2018). Health inequities and clustering.
# Diggle PJ, Giorgi E (2019). Model-based Geostatistics. Chapman & Hall.
# Diggle PJ, Giorgi E (2021). J R Soc Interface 18: 20210104.
# Wilson K, Wakefield J (2020). Biostatistics 21: e17-e32.
# Fuglstad GA et al. (2019). J Am Stat Assoc 114: 445-452.
# Lindgren F et al. (2011). J R Stat Soc B 73: 423-498.
# Fox J, Monette G (1992). J Am Stat Assoc 87: 178-183.
# Hosmer DW et al. (2013). Applied Logistic Regression. 3rd ed. Wiley.
# Krainski ET et al. (2019). Advanced Spatial Modeling with SPDE. CRC.
# ════════════════════════════════════════════════════════════════════════════════

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

final_DHS_data <- readRDS(
  file.path(base_path, "Data/DHS_final/final_DHS_data.rds")
)


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
all_surveys <- bind_rows(final_DHS_data, .id = "survey_idx")
rm(final_DHS_data)  # free memory

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

cat("\n----- SECTION 4: Filtering for eligible children -----\n")

eligible <- all_surveys %>%
  filter(
    !is.na(child_age_months), child_age_months >= 12, child_age_months <= 23,
    !is.na(child_vacc_measles),
    !is.na(cluster_psu),
    !is.na(LATNUM), LATNUM != 0,
    !is.na(LONGNUM), LONGNUM != 0
  ) %>%
  mutate(
    # Binary MCV1 indicator: vaccinated from any source
    # (card dated = 1, mother's report = 2, card marked = 3 in DHS coding)
    mcv1_received = as.integer(child_vacc_measles %in% c(1, 2, 3))
  )

cat("Eligible children (12-23 months with MCV1 data & GPS):", nrow(eligible), "\n")
cat("Countries:", n_distinct(eligible$CountryName), "\n")

# Save intermediate file for reproducibility
data.table::fwrite(eligible, paste0(base_path, "/Data/DHS_final/model_data_v2.csv"))
model_data <- vroom::vroom(paste0(base_path, "/Data/DHS_final/model_data_v2.csv"))






#model_data <- model_data %>% filter(CountryName == "Nigeria")





rm(eligible, all_surveys)  # free memory

# Add year-matched GC columns
gc_col_names <- paste0("gc_", tolower(gc_vars_base))
gc_results   <- lapply(gc_vars_base, function(gc) pick_gc_year(model_data, gc))
names(gc_results) <- gc_col_names
model_data <- cbind(model_data, as.data.frame(gc_results))
gc_available <- gc_col_names[gc_col_names %in% names(model_data)]

# Save with GC columns added
data.table::fwrite(model_data, paste0(base_path, "/Data/DHS_final/model_data_v3.csv"))
model_data <- vroom::vroom(paste0(base_path, "/Data/DHS_final/model_data_v3.csv"))
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
data.table::fwrite(cluster_data, paste0(base_path, "/Data/DHS_final/cluster_data_v2.csv"))
cluster_data <- vroom::vroom(paste0(base_path, "/Data/DHS_final/cluster_data_v2.csv"))


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



# ═════════════════════════════════════════════════════════════════•════
# ═══  HELPER FUNCTION: Prediction metrics                                      ═══
# ═════════════════════════════════════════════════════════════════•════
#
# Computes three numbers to judge how well a model predicts coverage:
#   r   = Pearson correlation (do predicted and observed track together?)
#   R²  = proportion of variance explained (0 = useless, 1 = perfect)
#   MAE = mean absolute error in coverage units (e.g. 0.10 = 10pp off)

calc_pred_metrics <- function(observed, predicted) {
  ok <- !is.na(observed) & !is.na(predicted)
  if (sum(ok) < 10) return(list(r = NA, r2 = NA, mae = NA))
  obs  <- observed[ok]
  pred <- predicted[ok]
  r    <- cor(obs, pred)
  ss_res <- sum((obs - pred)^2)
  ss_tot <- sum((obs - mean(obs))^2)
  r2   <- ifelse(ss_tot > 0, 1 - ss_res / ss_tot, NA)
  mae  <- mean(abs(obs - pred))
  list(r = round(r, 4), r2 = round(r2, 4), mae = round(mae, 4))
}
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




make_country_surface_diagnostics <- function(country_name, adm1_sf, sub_m, vars_final,
                                             inla_fu, mesh, spde, lcrs, cell_km = 5) {
  if (is.null(adm1_sf) || is.null(inla_fu) || is.null(mesh) || is.null(spde) || is.null(lcrs)) {
    return(NULL)
  }
  
  country_poly_ll <- tryCatch(
    adm1_sf %>%
      st_transform(4326) %>%
      st_union() %>%
      st_as_sf(),
    error = function(e) NULL
  )
  if (is.null(country_poly_ll)) return(NULL)
  
  country_poly_km <- tryCatch(
    st_transform(country_poly_ll, crs = lcrs),
    error = function(e) NULL
  )
  if (is.null(country_poly_km)) return(NULL)
  
  grid_km <- tryCatch({
    g <- st_make_grid(
      country_poly_km,
      cellsize = cell_km,
      square = TRUE,
      what = "polygons"
    )
    g <- st_as_sf(g) %>%
      st_filter(country_poly_km, .predicate = st_intersects)
    g$cell_id <- seq_len(nrow(g))
    g
  }, error = function(e) NULL)
  
  if (is.null(grid_km) || nrow(grid_km) < 10) return(NULL)
  
  grid_cent_km <- st_centroid(grid_km)
  grid_cent_ll <- st_transform(grid_cent_km, 4326)
  
  grid_ll <- st_coordinates(grid_cent_ll)
  grid_cp <- st_coordinates(grid_cent_km)
  
  A_grid <- INLA::inla.spde.make.A(mesh = mesh, loc = grid_cp)
  spatial_mean <- as.vector(A_grid %*% inla_fu$summary.random$spatial$mean)
  
  grid_df <- tibble(
    cell_id = grid_km$cell_id,
    LONGNUM = grid_ll[, 1],
    LATNUM = grid_ll[, 2],
    spatial_field = spatial_mean
  )
  
  map_infos <- purrr::map(vars_final, covariate_map_info, data_names = names(sub_m))
  needed_cols <- unique(c(vars_final, purrr::map_chr(map_infos, "map_var")))
  needed_cols <- needed_cols[needed_cols %in% names(sub_m)]
  
  obs_sp <- sub_m %>%
    dplyr::select(LONGNUM, LATNUM, all_of(needed_cols))
  
  for (v in needed_cols) {
    obs_v <- obs_sp %>%
      dplyr::select(LONGNUM, LATNUM, all_of(v)) %>%
      drop_na()
    
    names(obs_v)[3] <- "value"
    
    if (nrow(obs_v) >= 5 && sd(obs_v$value, na.rm = TRUE) > 1e-8) {
      idw_fit <- tryCatch(
        gstat::idw(
          value ~ 1,
          locations = ~ LONGNUM + LATNUM,
          data = obs_v,
          newdata = as.data.frame(grid_df[, c("LONGNUM", "LATNUM")]),
          idp = 2
        ),
        error = function(e) NULL
      )
      
      grid_df[[v]] <- if (!is.null(idw_fit)) {
        idw_fit$var1.pred
      } else {
        median(obs_v$value, na.rm = TRUE)
      }
    } else if (nrow(obs_v) > 0) {
      grid_df[[v]] <- median(obs_v$value, na.rm = TRUE)
    } else {
      grid_df[[v]] <- NA_real_
    }
  }
  
  fix <- inla_fu$summary.fixed$mean
  names(fix) <- rownames(inla_fu$summary.fixed)
  
  lp <- rep(
    ifelse("intercept" %in% names(fix), fix[["intercept"]], 0),
    nrow(grid_df)
  )
  
  for (v in vars_final) {
    if (v %in% names(fix) && v %in% names(grid_df)) {
      lp <- lp + fix[[v]] * grid_df[[v]]
    }
  }
  
  grid_df$predicted_coverage <- plogis(lp + grid_df$spatial_field)
  
  grid_map_km <- left_join(grid_km, grid_df, by = "cell_id")
  boundary_km <- st_boundary(country_poly_km)
  
  map_theme <- theme_minimal(base_size = 9) +
    theme(
      axis.title = element_blank(),
      panel.grid = element_line(colour = "grey88", linewidth = 0.2),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 7, colour = "grey35")
    )
  
  p_pred <- ggplot() +
    geom_sf(data = grid_map_km, aes(fill = predicted_coverage), colour = NA) +
    geom_sf(data = boundary_km, colour = "grey25", linewidth = 0.25) +
    scale_fill_viridis_c(
      name = "Predicted\ncoverage",
      limits = c(0, 1)
    ) +
    labs(
      title = paste(country_name, "- predicted MCV1 coverage"),
      subtitle = paste0(cell_km, " km x ", cell_km, " km grid; full fitted model")
    ) +
    coord_sf(datum = NA) +
    map_theme
  
  p_spatial <- ggplot() +
    geom_sf(data = grid_map_km, aes(fill = spatial_field), colour = NA) +
    geom_sf(data = boundary_km, colour = "grey25", linewidth = 0.25) +
    scale_fill_gradient2(
      name = "Spatial\nfield",
      low = "#d73027",
      mid = "white",
      high = "#4575b4",
      midpoint = 0
    ) +
    labs(
      title = paste(country_name, "- posterior mean spatial field"),
      subtitle = paste0(cell_km, " km x ", cell_km, " km grid; logit scale")
    ) +
    coord_sf(datum = NA) +
    map_theme
  
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
        name = info$legend,
        limits = if (is_prop && is.null(lims)) c(0, 1) else lims,
        labels = if (is_prop) scales::label_percent(accuracy = 1) else waiver(),
        oob = scales::squish,
        na.value = "grey90"
      ) +
      labs(
        title = info$label,
        subtitle = subtitle
      ) +
      coord_sf(datum = NA) +
      map_theme
  })
  
  names(covariate_plots) <- vars_final
  
  list(
    grid = st_drop_geometry(grid_map_km),
    surface_plot = p_pred / p_spatial,
    covariate_plots = covariate_plots
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

# Storage for results
all_results   <- list()
all_sheets    <- list()
all_comp      <- list()
all_pdf_paths <- list()
all_xl_paths  <- list()

# Per-country PDF files are opened inside the country loop after each model is fitted.


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# BEGIN COUNTRY LOOP
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

for (cname in countries[1]) {
  
  cat("\n\n")
  cat("═════════════════════════════════════════════════════════\n")
  cat("═══  COUNTRY:", formatC(cname, width = 48, flag = "-"),         "═══\n")
  cat("═════════════════════════════════════════════════════════\n")
  
  
  # ══════════════════════════════════════════════════════════════════
  # ═══  8.0  SUBSET & COUNTRY-LEVEL FACTOR PROPORTIONS                        ═══
  # ══════════════════════════════════════════════════════════════════
  #
  # First, extract this country's data and print descriptive statistics
  # including the proportions of key factor variables.
  # This tells us about the study population: how educated, wealthy,
  # urban, etc. the mothers in this country are.
  
  sub <- cluster_data %>% filter(CountryName == cname, n_children > 0)
  
  cat("Total clusters:", nrow(sub), "\n")
  cat("Total children:", sum(sub$n_children), "\n")
  cat("Overall MCV1 coverage:",
      round(sum(sub$n_vaccinated) / sum(sub$n_children) * 100, 1), "%\n")
  cat("Clusters at 0%:", sum(sub$coverage == 0),
      "(", round(mean(sub$coverage == 0) * 100, 1), "%)\n")
  cat("Clusters at 100%:", sum(sub$coverage == 1),
      "(", round(mean(sub$coverage == 1) * 100, 1), "%)\n")
  
  if (nrow(sub) < 50) { cat("  Too few clusters  skipping.\n"); next }
  
  # ----- Country-level coverage by cluster size --------------------------------------------------------------------------------
  cat("\n  Coverage by cluster size for", cname, ":\n")
  cov_by_size_country <- sub %>%
    mutate(sz = cut(n_children, c(0, 3, 5, 10, 20, Inf),
                    labels = c("1-3", "4-5", "6-10", "11-20", "21+"))) %>%
    group_by(sz) %>%
    summarise(
      n_clusters    = n(),
      pct_zero      = round(mean(coverage == 0) * 100, 1),
      pct_one       = round(mean(coverage == 1) * 100, 1),
      mean_coverage = round(mean(coverage) * 100, 1),
      .groups       = "drop"
    )
  print(cov_by_size_country, n = Inf)
  
  
  # ══════════════════════════════════════════════════════════════════
  # ═══  STEP 1: MISSINGNESS SCREENING                                          ═══
  # ══════════════════════════════════════════════════════════════════
  #
  # Utazi et al. (2022) Supplementary: "We excluded variables with greater
  # than 5% missing data from the analyses due to the possibility that
  # missingness is biased in an unknown way."
  #
  # WHY THIS MATTERS (for non-statisticians):
  # If a variable is missing for >5% of clusters, we can't be sure the
  # missingness is random.  For example, if only urban clusters have data
  # on internet use, including it would bias the model toward urban areas.
  # Dropping high-missingness variables avoids this unknown bias.
  #
  # We also drop variables with zero variance (all identical values in
  # this country  e.g. if ALL clusters are rural, urban_rural is useless).
  
  cat("\n--- Step 1: Missingness screening (Utazi 5% rule) ---\n")
  
  cols_to_remove <- c("cl_womens_autonomy_has_a_say", "cl_gc_nightlights_composite",  "cl_gc_elevation",
                      "cl_travel_time_to_city",       "cl_travel_time_to_HC_motor",
                     "cl_travel_time_to_HC_walk", "cl_gc_travel_times") # unstandardised
  
  
  avail_vars <- cl_vars_model[cl_vars_model %in% names(sub)]
  avail_vars <- avail_vars[!(avail_vars %in% cols_to_remove)]
  
  miss_pct <- sapply(sub[avail_vars], function(x) mean(is.na(x)) * 100)
  var_sd   <- sapply(sub[avail_vars], function(x) sd(x, na.rm = TRUE))
  
  vars_pass <- names(miss_pct[miss_pct <= 5 & !is.na(var_sd) & var_sd > 1e-8])
  vars_fail_miss <- names(miss_pct[miss_pct > 5])
  vars_fail_var  <- names(var_sd[!is.na(var_sd) & var_sd <= 1e-8])
  
  cat("  Passed (missingness + nonzero variance):",
      length(vars_pass), "of", length(avail_vars), "\n")
  if (length(vars_fail_miss) > 0)
    cat("  Dropped (>5% missing):", length(vars_fail_miss), "",
        paste(head(vars_fail_miss, 5), collapse = ", "),
        if (length(vars_fail_miss) > 5) "..." else "", "\n")
  if (length(vars_fail_var) > 0)
    cat("  Dropped (zero variance):", length(vars_fail_var), "",
        paste(head(vars_fail_var, 5), collapse = ", "),
        if (length(vars_fail_var) > 5) "..." else "", "\n")
  
  if (length(vars_pass) < 3) {
    cat("  Too few variables  skipping country.\n"); next
  }
  
  # Build an audit trail so we can track exactly what happened to each variable
  var_audit <- tibble(
    variable        = avail_vars,
    missingness_pct = round(miss_pct[avail_vars], 2),
    sd              = round(var_sd[avail_vars], 6),
    step1_status    = case_when(
      miss_pct[avail_vars] > 5 ~ "DROPPED: >5% missing",
      is.na(var_sd[avail_vars]) | var_sd[avail_vars] <= 1e-8 ~
        "DROPPED: zero/near-zero variance",
      TRUE ~ "PASSED"
    )
  )
  
  
  # ══════════════════════════════════════════════════════════════════
  # ═══  STEP 2: BIVARIATE SCREENING (crude odds ratios, keep p < 0.2)         ═══
  # ══════════════════════════════════════════════════════════════════
  #
  # For each surviving variable, we fit a simple one-variable model:
  #   "Is this variable, BY ITSELF, associated with MCV1 coverage?"
  #
  # We compute the CRUDE ODDS RATIO (OR):
  #   OR > 1 †’ higher values of the variable †’ MORE vaccination
  #   OR < 1 †’ higher values of the variable †’ LESS vaccination
  #   OR = 1 †’ no association
  #
  # We keep variables with p < 0.2 (a deliberately LIBERAL threshold).
  # This follows Hosmer & Lemeshow (2013, Applied Logistic Regression):
  # "Use a significance level of 0.20-0.25 for entering variables into the
  # model, as use of a more traditional level (e.g. 0.05) often fails to
  # identify variables known to be important."
  #
  # The idea is: we'd rather include a marginally important variable now
  # and let the multivariate model sort it out, than exclude it prematurely.
  
  cat("\n--- Step 2: Bivariate screening (crude ORs, keep p < 0.2) ---\n")
  
  biv <- map_dfr(vars_pass, function(v) {
    x <- sub[[v]]
    if (all(is.na(x)) || sd(x, na.rm = TRUE) < 1e-10)
      return(tibble(variable = v, estimate = NA, p_value = NA, status = "zero_var"))
    tryCatch({
      # Fit a binomial GLM: is coverage associated with this one variable?
      fit <- glm(cbind(n_vaccinated, n_children - n_vaccinated) ~ x,
                 family = binomial(link = "logit"), data = sub)
      s <- summary(fit)$coefficients
      if (nrow(s) < 2)
        return(tibble(variable = v, estimate = NA, p_value = NA, status = "singular"))
      tibble(
        variable = v,
        estimate = s[2, 1],
        se       = s[2, 2],
        z        = s[2, 3],
        p_value  = s[2, 4],
        OR       = exp(s[2, 1]),              # Odds ratio
        OR_lo    = exp(s[2, 1] - 1.96 * s[2, 2]),  # Lower 95% CI
        OR_hi    = exp(s[2, 1] + 1.96 * s[2, 2]),  # Upper 95% CI
        status   = "ok"
      )
    }, error = function(e) {
      tibble(variable = v, estimate = NA, p_value = NA,
             status = paste("error:", e$message))
    })
  })
  
  biv_pass <- biv %>%
    filter(!is.na(p_value), p_value < 0.2) %>%
    arrange(p_value) %>%
    pull(variable)
  
  cat("  Variables significant at p < 0.2:", length(biv_pass),
      "of", sum(!is.na(biv$p_value)), "tested\n")
  
  if (nrow(biv %>% filter(!is.na(p_value), p_value < 0.05)) > 0) {
    cat("  Top variables (p < 0.05):\n")
    print(biv %>% filter(!is.na(p_value), p_value < 0.05) %>%
            dplyr::select(variable, estimate, p_value, OR) %>%
            mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
            head(15),
          n = Inf)
  }
  
  if (length(biv_pass) < 2) {
    cat("  Fewer than 2 significant  using all passing vars.\n")
    biv_pass <- vars_pass
  }
  
  # Update audit trail
  var_audit <- var_audit %>%
    left_join(
      biv %>% dplyr::select(variable, estimate, p_value, OR, status) %>%
        rename(biv_estimate = estimate, biv_p_value = p_value,
               biv_OR = OR, biv_status = status),
      by = "variable"
    ) %>%
    mutate(
      step2_status = case_when(
        step1_status != "PASSED" ~ "N/A (dropped at Step 1)",
        is.na(biv_p_value) ~ paste0("DROPPED: model failed (", biv_status, ")"),
        biv_p_value >= 0.2 ~ paste0("DROPPED: p=", round(biv_p_value, 4), " >= 0.2"),
        TRUE ~ "PASSED"
      )
    )
  
  
  # ══════════════════════════════════════════════════════════════════
  # ═══  STEP 3: MULTICOLLINEARITY CHECK                                        ═══
  # ══════════════════════════════════════════════════════════════════
  #
  # If two variables measure essentially the same thing (e.g. wealth_score
  # and nighttime lights are both proxies for economic development), including
  # both in the model causes problems:
  #   - The model can't tell which one is "really" driving coverage
  #   - Coefficient estimates become unstable (huge standard errors)
  #   - Interpretation becomes misleading
  #
  # We check for this in two ways:
  #
  # Step 3a: PAIRWISE CORRELATION
  #   If two variables have |correlation| > 0.8, we drop the one with the
  #   weaker bivariate association (higher p-value).
  #   Example: if wealth_score and nightlights have r = 0.85, and wealth_score
  #   had a stronger crude association with coverage, we keep wealth_score
  #   and drop nightlights.
  #
  # Step 3b: GENERALISED VARIANCE INFLATION FACTOR (GVIF)
  #   Even after removing high pairwise correlations, there can be
  #   "multicollinearity" where a variable is predicted by the COMBINATION
  #   of other variables.  GVIF measures this.  Following Utazi et al. (2022):
  #   "We addressed the problem of (multi)collinearity by computing the
  #   generalized variance inflation factors (GVIFs) and excluded variables
  #   that had high GVIFs (> 2, on the scale that ensures comparability
  #   across the covariates as recommended by Fox and Monette)."
  #   We iteratively remove the worst offender until all GVIF^(1/(2*df)) ‰¤ 2.
  
  cat("\n--- Step 3: Multicollinearity check ---\n")
  
  sub_mat <- sub %>% dplyr::select(all_of(biv_pass)) %>% drop_na()
  cor_drop_log <- tibble()
  vif_drop_log <- tibble()
  
  if (nrow(sub_mat) < 30 || ncol(sub_mat) < 2) {
    vars_final <- biv_pass
    cat("  Insufficient data for correlation check  using all bivariate vars.\n")
  } else {
    # Step 3a: pairwise correlations
    cor_mat <- cor(sub_mat, use = "complete.obs")
    
    drop_cor <- c()
    for (i in 1:(ncol(cor_mat) - 1)) {
      for (j in (i + 1):ncol(cor_mat)) {
        if (abs(cor_mat[i, j]) > 0.8) {
          v1 <- colnames(cor_mat)[i]; v2 <- colnames(cor_mat)[j]
          p1 <- biv$p_value[biv$variable == v1]
          p2 <- biv$p_value[biv$variable == v2]
          p1 <- ifelse(length(p1) == 0 || is.na(p1), 1, p1)
          p2 <- ifelse(length(p2) == 0 || is.na(p2), 1, p2)
          drop_v <- ifelse(p1 > p2, v1, v2)
          keep_v <- ifelse(p1 > p2, v2, v1)
          drop_cor <- c(drop_cor, drop_v)
          cor_drop_log <- bind_rows(cor_drop_log, tibble(
            var1 = v1, var2 = v2,
            correlation = round(cor_mat[i, j], 3),
            dropped = drop_v, kept = keep_v,
            reason = paste0("|r|=", round(abs(cor_mat[i, j]), 3),
                            "; kept ", keep_v, " (lower p)")
          ))
          cat("  |r| > 0.8:", v1, "&", v2,
              "(r =", round(cor_mat[i, j], 3), ") drop", drop_v, "\n")
        }
      }
    }
    vars_after_cor <- setdiff(biv_pass, unique(drop_cor))
    
    # Step 3b: iterative GVIF
    vars_final <- vars_after_cor
    if (length(vars_final) >= 2) {
      repeat {
        sub_vif <- sub %>%
          dplyr::select(n_vaccinated, n_children, all_of(vars_final)) %>%
          drop_na()
        if (nrow(sub_vif) < 30) break
        
        vif_form <- as.formula(paste(
          "cbind(n_vaccinated, n_children - n_vaccinated) ~",
          paste(vars_final, collapse = " + ")
        ))
        vif_fit <- tryCatch(glm(vif_form, binomial, sub_vif), error = function(e) NULL)
        if (is.null(vif_fit)) break
        
        vif_vals <- tryCatch(car::vif(vif_fit), error = function(e) NULL)
        if (is.null(vif_vals)) break
        
        if (is.matrix(vif_vals)) {
          gvif_adj <- vif_vals[, "GVIF^(1/(2*Df))"]
        } else {
          gvif_adj <- sqrt(vif_vals)
        }
        
        if (max(gvif_adj) <= 2) break
        
        drop_var <- names(which.max(gvif_adj))
        vif_drop_log <- bind_rows(vif_drop_log, tibble(
          variable       = drop_var,
          gvif_adj_value = round(max(gvif_adj), 3),
          reason         = paste0("GVIF^(1/2df) = ", round(max(gvif_adj), 3),
                                  " > 2.0 threshold")
        ))
        cat("  Dropping (GVIF^(1/2df) =", round(max(gvif_adj), 2), "):",
            drop_var, "\n")
        vars_final <- setdiff(vars_final, drop_var)
        
        if (length(vars_final) < 2) break
      }
    }
  }
  
  # Update audit trail
  var_audit <- var_audit %>%
    mutate(
      step3_status = case_when(
        step1_status != "PASSED" ~ "N/A (dropped at Step 1)",
        step2_status != "PASSED" ~ "N/A (dropped at Step 2)",
        variable %in% unique(drop_cor) ~
          paste0("DROPPED: pairwise |r| > 0.8 with ",
                 cor_drop_log$kept[match(variable, cor_drop_log$dropped)]),
        variable %in% vif_drop_log$variable ~
          paste0("DROPPED: ", vif_drop_log$reason[match(variable,
                                                        vif_drop_log$variable)]),
        variable %in% vars_final ~ "PASSED IN FINAL MODEL",
        TRUE ~ "DROPPED (reason unclear)"
      )
    )
  
  cat("\n  Final covariate set for", cname, "(", length(vars_final), "variables):\n")
  cat(paste("    ", vars_final), sep = "\n")
  
  
  # ----- Prepare complete-case modelling dataset --------------------------------------------------------------------------------
  raw_vars_for_mapping <- unique(raw_covariate_name(vars_final))
  
  sub_m <- sub %>%
    dplyr::select(
      LONGNUM, LATNUM, n_children, n_vaccinated, coverage,
      coverage_emp_logit,
      all_of(vars_final),
      any_of(raw_vars_for_mapping)
    ) %>%
    drop_na() %>%
    filter(n_children > 0)
  
  
  cat("\n  Complete cases for modelling:", nrow(sub_m), "\n")
  if (nrow(sub_m) < 30) { cat("   ï¸ Too few  skipping.\n"); next }
  
  
  # ══════════════════════════════════════════════════════════════════
  # ═══  MODEL 1: COVARIATES-ONLY BINOMIAL GLM (non-spatial baseline)          ═══
  # ══════════════════════════════════════════════════════════════════
  #
  # This is a standard logistic regression on the cluster-level data.
  # It asks: "Can we predict vaccination coverage from socioeconomic and
  # demographic factors alone, ignoring WHERE the cluster is located?"
  #
  # The model equation is:
  #   logit(p_i) = Î²‚€ + Î²‚ Ã— education + Î²‚‚ Ã— wealth + ... + Î²‚– Ã— conflict
  #
  # where p_i = probability of vaccination in cluster i, and the Î²s tell us
  # how each factor relates to coverage (as odds ratios, OR = exp(Î²)).
  #
  # This is the NON-SPATIAL model  it treats each cluster as independent,
  # ignoring that nearby clusters tend to have similar coverage.
  # It serves as a BASELINE for comparison with the spatial models.
  #
  # WHY START NON-SPATIAL?
  # Diggle & Giorgi (2021) recommend building covariates first without spatial
  # effects.  If we added the spatial field right away, it could "absorb"
  # the signal that should be attributed to covariates, making them look
  # less important than they really are.
  
  cat("\n--- Model 1: Covariates-only binomial GLM ---\n")
  
  glm_form <- as.formula(paste(
    "cbind(n_vaccinated, n_children - n_vaccinated) ~",
    paste(vars_final, collapse = " + ")
  ))
  
  fit_glm <- tryCatch(
    glm(glm_form, family = binomial(link = "logit"), data = sub_m),
    error = function(e) { cat("  GLM failed:", e$message, "\n"); NULL }
  )
  
  glm_co <- NULL; auc_g <- NA
  if (!is.null(fit_glm)) {
    glm_co <- as.data.frame(summary(fit_glm)$coefficients) %>%
      tibble::rownames_to_column("term") %>%
      mutate(
        OR    = exp(Estimate),
        OR_lo = exp(Estimate - 1.96 * `Std. Error`),
        OR_hi = exp(Estimate + 1.96 * `Std. Error`),
        sig   = case_when(
          `Pr(>|z|)` < 0.001 ~ "***",
          `Pr(>|z|)` < 0.01  ~ "**",
          `Pr(>|z|)` < 0.05  ~ "*",
          `Pr(>|z|)` < 0.1   ~ ".",
          TRUE                ~ ""
        )
      )
    
    sub_m$pred_glm  <- predict(fit_glm, type = "response")
    sub_m$resid_glm <- sub_m$coverage - sub_m$pred_glm
    
    auc_g <- tryCatch({
      roc_obj <- pROC::roc(
        response  = as.numeric(sub_m$coverage > 0.5),
        predictor = sub_m$pred_glm, quiet = TRUE)
      as.numeric(roc_obj$auc)
    }, error = function(e) NA_real_)
    
    cat("  AIC:", round(AIC(fit_glm), 1), "\n")
    cat("  AUC:", round(auc_g, 3), "\n")
    cat("  Significant covariates (p < 0.05):\n")
    sig_co <- glm_co %>% filter(`Pr(>|z|)` < 0.05, term != "(Intercept)")
    if (nrow(sig_co) > 0) {
      print(sig_co %>% dplyr::select(term, Estimate, OR, sig) %>%
              mutate(across(where(is.numeric), ~ round(.x, 3))),
            row.names = FALSE)
    }
  }
  
  
  # ══════════════════════════════════════════════════════════════════
  # ═══  MODELS 2-4: INLA-SPDE GEOSTATISTICAL MODELS                           ═══
  # ══════════════════════════════════════════════════════════════════
  #
  # Now we add GEOGRAPHY to the model.  The key insight from geostatistical
  # modelling is that nearby clusters tend to have similar coverage  this
  # is called SPATIAL AUTOCORRELATION.  It arises because nearby clusters
  # share unmeasured factors like:
  #   - Quality of the local health facility
  #   - Road infrastructure
  #   - Local community leaders' attitudes
  #   - Active conflict in the area
  #
  # The INLA-SPDE approach (Lindgren et al. 2011) approximates a continuous
  # spatial surface using a triangular mesh.  Think of it as laying a
  # flexible rubber sheet over the map and letting it warp to fit the data.
  #
  # The full model equation (Utazi et al. 2018, Equation 1):
  #
  #   y_i | p(s_i) ~ Binomial(n_i, p(s_i))
  #   logit(p(s_i)) = X_i Î² + S(s_i) + Z_i
  #
  # where:
  #   y_i = number of vaccinated children in cluster i
  #   n_i = total number of eligible children in cluster i
  #   p(s_i) = probability of vaccination at location s_i
  #   X_i Î² = covariate effects (education, wealth, etc.)
  #   S(s_i) = spatial random effect (the "rubber sheet")
  #   Z_i = iid nugget (non-spatial noise per cluster)
  #
  # WHY THREE SPATIAL MODELS?
  #   Model 2 (spatial-only): just the rubber sheet + nugget, NO covariates
  #     †’ tells us "how much does location alone explain?"
  #   Model 3 (full): covariates + rubber sheet + nugget
  #     †’ the PRIMARY model, gives us both "why" and "where"
  #   Model 4 (beta-binomial): same as Model 3 but allows extra variation
  #     †’ tests if the 0/1 boundary problem is severe
  #
  # KEY SPATIAL PARAMETERS:
  #   Range: the distance (in km) over which spatial correlation drops to
  #     ~13%.  A range of 100 km means clusters >100 km apart are
  #     essentially independent.  Large range †’ regional-scale patterns;
  #     small range †’ local-scale patterns.
  #   Spatial variance: how much variation the spatial field explains on
  #     the logit scale.  High †’ covariates miss important spatial patterns;
  #     low †’ covariates capture most geographic variation.
  
  
  # IN VERSION 2 V2
  # So the order I€™d try:
  # 
  # Relax iid prior: param = c(2, 0.05) or c(2, 0.1)
  # Relax spatial sigma: prior.sigma = c(2, 0.05)
  # If still too smooth, reduce prior.range
  
  
  
  
  inla_sp <- NULL; inla_fu <- NULL; inla_bb <- NULL
  rp <- NA; sp <- NA; feff <- NULL
  lcrs <- NULL; mesh <- NULL; spde <- NULL
  
  if (nrow(sub_m) >= 50 && inla_available) {
    
    cat("\n--- Models 2-4: INLA-SPDE geostatistical models ---\n")
    
    # Project GPS coordinates to a local coordinate system (in km)
    coords <- cbind(sub_m$LONGNUM, sub_m$LATNUM)
    cent   <- colMeans(coords)
    lcrs   <- paste0("+proj=aeqd +lat_0=", cent[2], " +lon_0=", cent[1],
                     " +x_0=0 +y_0=0 +datum=WGS84 +units=km")
    
    cp <- tryCatch(
      sf::sf_project(st_crs(4326)$proj4string, lcrs, coords),
      error = function(e) { cat("  Projection failed:", e$message, "\n"); NULL }
    )
    
    if (!is.null(cp)) {
      
      # Build the SPDE mesh
      ext <- max(diff(range(cp[, 1])), diff(range(cp[, 2])))

      # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      # Sensitivity analysis on spatial parameters
      # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      
      
      # cutoff   = max(rg / 20, 5),
      # offset   = c(rg / 2, rg * 1.5)          
      
      
      # V1
      
      # version <- "_(v1)"
      # 
      # output_folder <- file.path(base_path, paste0("TABLES/Statical_Model_Results", version))
      # if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
      # 
      # rg <- max(ext / 3, 20)
      # max_edge <- c(rg / 5, rg)
      # offset_value <- c(rg / 2, rg * 1.5)
      # cutoff_value <- max(rg / 20, 5)
      # prior_on_range <- c(rg, 0.5)
      # prior_on_variance <- c(1, 0.01)
      # nugget_effect <- "c(1, 0.01)"
      
      # V2
      
      version <- "_(v2)"

      output_folder <- file.path(base_path, paste0("TABLES/Statical_Model_Results", version))
      if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

      rg <- max(ext / 3, 20)
      max_edge <- c(rg / 5, rg)
      offset_value <- c(rg / 2, rg * 1.5)
      cutoff_value <- max(rg / 20, 5)
      prior_on_range <- c(rg/2, 0.5)
      prior_on_variance <- c(2, 0.05)
      nugget_effect <- "c(2, 0.05)"
      
      # V3
      
      
      # version <- "_(v3)"
      # 
      # output_folder <- file.path(base_path, paste0("TABLES/Statical_Model_Results", version))
      # if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
      # 
      # rg <- max(ext / 3, 20)
      # max_edge <- c(rg / 5, rg)
      # offset_value <- c(rg / 2, rg * 1.5)
      # cutoff_value <- max(rg / 20, 5)
      # prior_on_range <- c(rg*2, 0.5)
      # prior_on_variance <- c(1, 0.005)
      # nugget_effect <- "c(0.5, 0.1)"
      
      # V4
      
      # version <- "_(v4)"
      # 
      # output_folder <- file.path(base_path, paste0("TABLES/Statical_Model_Results", version))
      # if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
      # 
      # rg <- max(ext / 5, 10)
      # max_edge <- c(rg / 8, rg/2)
      # offset_value <- c(rg / 3, rg)
      # cutoff_value <- max(rg / 30, 2)
      # prior_on_range <- c(rg, 0.5)
      # prior_on_variance <- c(1, 0.01)
      # nugget_effect <- "c(1, 0.01)"
      
      # V5
      
      # version <- "_(v5)"
      # 
      # output_folder <- file.path(base_path, paste0("TABLES/Statical_Model_Results", version))
      # if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
      # 
      # rg <- max(ext / 2, 50)
      # max_edge <- c(rg / 3, rg)
      # offset_value <- c(rg, 2*rg)
      # cutoff_value <- max(rg / 10, 20)
      # prior_on_range <- c(rg, 0.5)
      # prior_on_variance <- c(1, 0.01)
      # nugget_effect <- "c(1, 0.01)"
      
      
      # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      # END OF ------ Sensitivity analysis on spatial parameters
      # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      
      
      mesh <- tryCatch(
        INLA::inla.mesh.2d(
          loc      = cp,
          max.edge = max_edge,
          cutoff   = cutoff_value,
          offset   = offset_value
        ),
        error = function(e) { cat("  Mesh failed:", e$message, "\n"); NULL }
      )
      
      if (!is.null(mesh)) {
        cat("  Mesh vertices:", mesh$n, "\n")
        
        # PC priors following Fuglstad et al. (2019):
        # "Constructing priors that penalize the complexity of Gaussian
        # random fields."  These priors shrink toward a simpler model
        # (no spatial effect) unless the data strongly support one.
        spde <- INLA::inla.spde2.pcmatern(
          mesh   = mesh,
          alpha  = 2,
          prior.range = prior_on_range,   # P(range < rg) = 0.5 AND IN V2  P(range < rg/2) = 0.25 TO SHORTEN SPATIAL RANGE IF THE FIELD IS TOO SMOOTH
          prior.sigma = prior_on_variance    # P(sigma > 1) = 0.01 AND TEST OUT IN V2 P(spatial sigma > 2) = 0.05
        )
        
        A <- INLA::inla.spde.make.A(mesh = mesh, loc = cp)
        cov_df <- sub_m %>% dplyr::select(all_of(vars_final))
        
        stk <- INLA::inla.stack(
          data = list(y = sub_m$n_vaccinated, n_trial = sub_m$n_children),
          A = list(A, 1),
          effects = list(
            spatial = 1:spde$n.spde,
            data.frame(intercept = 1, cluster_id = 1:nrow(sub_m), cov_df)
          ),
          tag = "estimation"
        )
        
        cov_terms <- paste(vars_final, collapse = " + ")
        iid_spec  <- paste0(
          "f(cluster_id, model = 'iid', ",
          "hyper = list(prec = list(prior = 'pc.prec', param = ", nugget_effect, ")))" # AND IN V2 param = c(2, 0.05)
        )
        
        run_inla <- function(formula, family) {
          tryCatch(
            INLA::inla(
              formula = formula,
              family  = family,
              Ntrials = INLA::inla.stack.data(stk)$n_trial,
              data    = INLA::inla.stack.data(stk, spde = spde),
              control.predictor = list(
                A = INLA::inla.stack.A(stk), compute = TRUE, link = 1
              ),
              control.compute = list(
                dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE
              ),
              control.inla = list(strategy = "adaptive"),
              verbose = FALSE
            ),
            error = function(e) {
              cat("    INLA failed:", e$message, "\n"); NULL
            }
          )
        }
        
        # MODEL 2: Spatial-only (no covariates  just location)
        cat("\n  Model 2: Spatial-only (intercept + Matern GP + iid)...\n")
        f2 <- as.formula(paste0(
          "y ~ -1 + intercept + f(spatial, model = spde) + ", iid_spec
        ))
        inla_sp <- run_inla(f2, "binomial")
        
        
        sp2 <- get_spatial_var(inla_sp, spde)
        
        sp_range <- sp2$range
        sp_var   <- sp2$var
        
        if (!is.null(inla_sp)) {
          cat("    DIC:", round(inla_sp$dic$dic, 1),
              "| WAIC:", round(inla_sp$waic$waic, 1), "\n")
          idx_sp <- INLA::inla.stack.index(stk, "estimation")$data
          sub_m$pred_sp  <- inla_sp$summary.fitted.values$mean[idx_sp]
          sub_m$resid_sp <- sub_m$coverage - sub_m$pred_sp
        }
        
        # MODEL 3: Full model (covariates + spatial + nugget) †═ PRIMARY
        cat("\n  Model 3: Full (covariates + Matern GP + iid)...\n")
        f3 <- as.formula(paste0(
          "y ~ -1 + intercept + ", cov_terms,
          " + f(spatial, model = spde) + ", iid_spec
        ))
        inla_fu <- run_inla(f3, "binomial")
        
        if (!is.null(inla_fu)) {
          cat("    DIC:", round(inla_fu$dic$dic, 1),
              "| WAIC:", round(inla_fu$waic$waic, 1), "\n")
          
          feff <- inla_fu$summary.fixed %>%
            mutate(
              term = rownames(.),
              sig  = ifelse(`0.025quant` > 0 | `0.975quant` < 0, "*", ""),
              OR    = exp(mean),
              OR_lo = exp(`0.025quant`),
              OR_hi = exp(`0.975quant`)
            ) %>%
            dplyr::select(term, mean, sd, `0.025quant`, `0.975quant`,
                          OR, OR_lo, OR_hi, sig)
          
          cat("\n  Fixed effects (full model):\n")
          print(feff %>% mutate(across(where(is.numeric), ~ round(.x, 3))),
                row.names = FALSE)
          
          # spr <- INLA::inla.spde2.result(inla_fu, "spatial", spde)
          # rp  <- INLA::inla.emarginal(function(x) x,
          #                             spr$marginals.range.nominal[[1]])
          # sp  <- INLA::inla.emarginal(function(x) x,
          #                             spr$marginals.variance.nominal[[1]])
          
          sp3 <- get_spatial_var(inla_fu, spde)
          
          rp <- sp3$range
          sp <- sp3$var
          
          
          cat(sprintf("\n  Spatial range: %.1f km\n", rp))
          cat(sprintf("  Spatial marginal variance: %.4f\n", sp))
          
          idx <- INLA::inla.stack.index(stk, "estimation")$data
          sub_m$pred_full  <- inla_fu$summary.fitted.values$mean[idx]
          sub_m$resid_full <- sub_m$coverage - sub_m$pred_full
        }
        
        # MODEL 4: Beta-binomial (overdispersion check)
        cat("\n  Model 4: Beta-binomial (overdispersion check)...\n")
        inla_bb <- run_inla(f3, "betabinomial")
        
        sp4 <- get_spatial_var(inla_bb, spde)
        
        bb_range <- sp4$range
        bb_var   <- sp4$var
        
        if (!is.null(inla_bb)) {
          cat("    DIC:", round(inla_bb$dic$dic, 1),
              "| WAIC:", round(inla_bb$waic$waic, 1), "\n")
          idx_bb <- INLA::inla.stack.index(stk, "estimation")$data
          sub_m$pred_bb  <- inla_bb$summary.fitted.values$mean[idx_bb]
          sub_m$resid_bb <- sub_m$coverage - sub_m$pred_bb
        }
        
      } # end if mesh
    } # end if projection
  } # end if INLA
  
  
  # ══════════════════════════════════════════════════════════════════
  # ═══  MODEL COMPARISON                                                        ═══
  # ══════════════════════════════════════════════════════════════════
  #
  # We compare all four models using metrics that a non-statistician can
  # interpret:
  #
  #   r (correlation): Do the model's predictions go up and down with the
  #     real data?  r close to 1 = good tracking.
  #
  #   R² (proportion of variance explained): What fraction of the variation
  #     in coverage does the model capture?  R² = 0.5 means the model
  #     explains half the variation; 0.7+ is considered good for this
  #     type of data.
  #
  #   MAE (mean absolute error): How far off are the predictions, on average?
  #     MAE = 0.10 means predictions are typically ±10 percentage points.
  #     MAE < 0.10 is good; < 0.15 is acceptable.
  #
  #   DIC / WAIC: Bayesian model comparison scores (lower = better).
  #     These balance goodness of fit against model complexity.
  #     A WAIC difference > 10 is strong evidence one model is better.
  #
  #   AUC: For the GLM only  how well does the model distinguish
  #     clusters with >50% coverage from those with <50%?
  #     AUC 0.8-0.9 = good; 0.9+ = excellent.
  
  cat("\n--- Model comparison ---\n")
  
  pm_glm  <- if (!is.null(fit_glm)) calc_pred_metrics(sub_m$coverage, sub_m$pred_glm)
  else list(r = NA, r2 = NA, mae = NA)
  pm_sp   <- if ("pred_sp" %in% names(sub_m)) calc_pred_metrics(sub_m$coverage, sub_m$pred_sp)
  else list(r = NA, r2 = NA, mae = NA)
  pm_full <- if ("pred_full" %in% names(sub_m)) calc_pred_metrics(sub_m$coverage, sub_m$pred_full)
  else list(r = NA, r2 = NA, mae = NA)
  pm_bb   <- if ("pred_bb" %in% names(sub_m)) calc_pred_metrics(sub_m$coverage, sub_m$pred_bb)
  else list(r = NA, r2 = NA, mae = NA)
  
  cat("\n  Prediction metrics (observed vs predicted cluster coverage):\n")
  pred_table <- tibble(
    Model = c("GLM (covariates only)", "Spatial-only INLA",
              "Full INLA (cov + spatial)", "Beta-binomial INLA"),
    r     = c(pm_glm$r, pm_sp$r, pm_full$r, pm_bb$r),
    R2    = c(pm_glm$r2, pm_sp$r2, pm_full$r2, pm_bb$r2),
    MAE   = c(pm_glm$mae, pm_sp$mae, pm_full$mae, pm_bb$mae)
  )
  print(pred_table, n = Inf)
  
  mc <- tibble(
    country      = cname,
    n_clusters   = nrow(sub_m),
    n_children   = sum(sub_m$n_children),
    n_covariates = length(vars_final),
    pct_zero     = round(mean(sub_m$coverage == 0) * 100, 1),
    pct_one      = round(mean(sub_m$coverage == 1) * 100, 1),
    glm_aic      = if (!is.null(fit_glm)) round(AIC(fit_glm), 1) else NA,
    glm_auc      = round(auc_g, 3),
    glm_r        = pm_glm$r,   
    glm_r2       = pm_glm$r2,   
    glm_mae      = pm_glm$mae,
    sp_dic       = if (!is.null(inla_sp)) round(inla_sp$dic$dic, 1) else NA,
    sp_waic      = if (!is.null(inla_sp)) round(inla_sp$waic$waic, 1) else NA,
    sp_r         = pm_sp$r,    
    sp_r2        = pm_sp$r2,    
    sp_mae       = pm_sp$mae,
    sp_range     = round(sp_range, 1),
    sp_var       = round(sp_var, 4),
    full_dic     = if (!is.null(inla_fu)) round(inla_fu$dic$dic, 1) else NA,
    full_waic    = if (!is.null(inla_fu)) round(inla_fu$waic$waic, 1) else NA,
    full_r       = pm_full$r,  
    full_r2      = pm_full$r2,  
    full_mae     = pm_full$mae,
    full_range   = round(rp, 1),
    full_var     = round(sp, 4),
    bb_dic       = if (!is.null(inla_bb)) round(inla_bb$dic$dic, 1) else NA,
    bb_waic      = if (!is.null(inla_bb)) round(inla_bb$waic$waic, 1) else NA,
    bb_r         = pm_bb$r,    
    bb_r2        = pm_bb$r2,    
    bb_mae       = pm_bb$mae,
    bb_range     = bb_range,
    bb_var       = bb_var,
    range_km     = round(rp, 1),
    spatial_var  = round(sp, 4)
  )
  mc$best_model <- case_when(
    !is.na(mc$full_waic) & !is.na(mc$bb_waic) & mc$bb_waic < mc$full_waic ~ "Beta-binomial",
    !is.na(mc$full_waic) ~ "Full INLA",
    !is.null(fit_glm) ~ "GLM",
    TRUE ~ "None"
  )
  
  all_comp[[cname]] <- mc
  cat("\n"); print(t(mc), quote = FALSE)
  
  
  # ══════════════════════════════════════════════════════════════════
  # ═══  DIAGNOSTIC PLOTS (PDF output  multiple pages per country)             ═══
  # ══════════════════════════════════════════════════════════════════
  
  # ----- Select best predictions for plotting --------------------------------------------------------------------------------
  if ("pred_full" %in% names(sub_m)) {
    sub_m$pred_best  <- sub_m$pred_full
    sub_m$resid_best <- sub_m$resid_full
    bl <- "Full INLA (Binomial + Matern GP + Nugget)"
  } else if ("pred_glm" %in% names(sub_m)) {
    sub_m$pred_best  <- sub_m$pred_glm
    sub_m$resid_best <- sub_m$resid_glm
    bl <- "Covariates-only Binomial GLM"
  } else {
    sub_m$pred_best <- NA; sub_m$resid_best <- NA; bl <- "None"
  }
  
  sub_sf <- st_as_sf(sub_m, coords = c("LONGNUM", "LATNUM"), crs = 4326)
  
  
  # ----- Get ADM1 boundaries from GADM ----------------------------------------------------------------------------------------------------
  country_iso <- tryCatch({
    countrycode::countrycode(cname, "country.name", "iso3c")
  }, error = function(e) NA_character_)
  
  adm1_sf <- NULL
  if (!is.na(country_iso)) {
    adm1_sf <- tryCatch(
      geodata::gadm(country = country_iso, level = 1, path = tempdir()) %>%
        st_as_sf(),
      error = function(e) NULL
    )
  }
  
  
  country_slug <- safe_file_name(cname)
  pdf_path <- file.path(output_folder, paste0(country_slug, "_geostatistical_models_unweighted", version, ".pdf"))
  pdf(pdf_path, width = 16, height = 20)
  surface_diag <- NULL
  residual_covariate_plots <- NULL
  
  # ----- PAGE 1: Country overview map ---------------------------------------------------------------------------------------------------------
  p_map <- ggplot()
  if (!is.null(adm1_sf)) {
    p_map <- p_map +
      geom_sf(data = adm1_sf, fill = NA, colour = "grey30", linewidth = 0.4)
  }
  p_map <- p_map +
    geom_sf(data = sub_sf, aes(colour = coverage), size = 2.0, alpha = 0.7) +
    scale_colour_viridis_c(name = "MCV1\nCoverage", limits = c(0, 1)) +
    labs(
      title    = paste(cname, " - DHS Cluster Locations with ADM1 Boundaries"),
      subtitle = paste0(nrow(sub_m), " clusters | Overall coverage: ",
                        round(sum(sub_m$n_vaccinated) / sum(sub_m$n_children) * 100, 1), "%"),
      caption  = "Source: GADM boundaries; DHS GPS coordinates (jittered)"
    ) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "right")
  print(p_map)
  
  
  # ----- PAGE 2: Model summary text page -----------------------------------------------------------------------------------------------
  # Build equation string
  if (!is.null(feff)) {
    eq_fixed <- paste0("logit(p(s_i)) = ",
                       round(feff$mean[feff$term == "intercept"], 3))
    for (rr in seq_len(nrow(feff))) {
      if (feff$term[rr] != "intercept") {
        sgn <- ifelse(feff$mean[rr] >= 0, " + ", " - ")
        eq_fixed <- paste0(eq_fixed, sgn, round(abs(feff$mean[rr]), 3),
                           " Ã— ", feff$term[rr])
      }
    }
    eq_full <- paste0(eq_fixed, "\n       + S(s_i) + Z_i")
    eq_label <- "Full INLA Geostatistical Model"
  } else if (!is.null(glm_co)) {
    eq_fixed <- paste0("logit(p_i) = ",
                       round(glm_co$Estimate[glm_co$term == "(Intercept)"], 3))
    for (rr in seq_len(nrow(glm_co))) {
      if (glm_co$term[rr] != "(Intercept)") {
        sgn <- ifelse(glm_co$Estimate[rr] >= 0, " + ", " - ")
        eq_fixed <- paste0(eq_fixed, sgn, round(abs(glm_co$Estimate[rr]), 3),
                           " Ã— ", glm_co$term[rr])
      }
    }
    eq_full <- eq_fixed
    eq_label <- "Covariates-only Binomial GLM"
  } else {
    eq_full <- "No model fitted"; eq_label <- "None"
  }
  
  summary_lines <- c(
    paste0("COUNTRY: ", cname),
    paste0("Best model: ", mc$best_model),
    "",
    "----- DATA SUMMARY -----",
    paste0("Clusters: ", nrow(sub_m), " | Children: ", sum(sub_m$n_children)),
    paste0("Clusters at 0%: ", mc$pct_zero, "% | Clusters at 100%: ", mc$pct_one, "%"),
    paste0("Overall coverage: ",
           round(sum(sub_m$n_vaccinated) / sum(sub_m$n_children) * 100, 1), "%"),
    "",
    "----- VARIABLE SELECTION PIPELINE -----",
    paste0("Step 1 (Missingness <5%): ",
           sum(var_audit$step1_status == "PASSED"), " of ",
           nrow(var_audit), " candidates"),
    paste0("Step 2 (Bivariate p<0.2): ",
           sum(var_audit$step2_status == "PASSED", na.rm = TRUE), " passed"),
    paste0("Step 3 (Multicollinearity): ", length(vars_final), " in final set"),
    "",
    "----- FINAL COVARIATES -----",
    paste(vars_final, collapse = ", "),
    "",
    "----- MODEL EQUATION -----",
    eq_label,
    eq_full,
    "",
    "----- MODEL FIT -----"
  )
  
  if (!is.null(fit_glm))
    summary_lines <- c(
      summary_lines,
      paste0(
        "GLM: AIC=", round(AIC(fit_glm), 1),
        " | AUC=", round(auc_g, 3)
      )
    )
  
  if (!is.null(inla_sp))
    summary_lines <- c(
      summary_lines,
      paste0(
        "Spatial-only: DIC=", mc$sp_dic,
        " | WAIC=", mc$sp_waic,
        " | Range=", mc$sp_range, " km",
        " | Variance=", mc$sp_var
      )
    )
  
  if (!is.null(inla_fu))
    summary_lines <- c(
      summary_lines,
      paste0(
        "Full INLA: DIC=", mc$full_dic,
        " | WAIC=", mc$full_waic,
        " | Range=", mc$full_range, " km",
        " | Variance=", mc$full_var
      )
    )
  
  if (!is.null(inla_bb))
    summary_lines <- c(
      summary_lines,
      paste0(
        "Beta-binomial: DIC=", mc$bb_dic,
        " | WAIC=", mc$bb_waic,
        " | Range=", round(mc$bb_range, 1), " km",
        " | Variance=", round(mc$bb_var, 4)
      )
    )
  
  # Prediction metrics
  summary_lines <- c(
    summary_lines, "",
    "----- PREDICTION METRICS -----",
    "  r = correlation | R² = variance explained | MAE = avg error",
    
    sprintf(
      "  %-30s  r=%s  R²=%s  MAE=%s",
      "GLM:",
      ifelse(is.na(pm_glm$r), " N/A", sprintf("%6.3f", pm_glm$r)),
      ifelse(is.na(pm_glm$r2), " N/A", sprintf("%6.3f", pm_glm$r2)),
      ifelse(is.na(pm_glm$mae), "N/A", sprintf("%5.3f", pm_glm$mae))
    ),
    
    sprintf(
      "  %-30s  r=%s  R²=%s  MAE=%s",
      "Spatial-only:",
      ifelse(is.na(pm_sp$r), " N/A", sprintf("%6.3f", pm_sp$r)),
      ifelse(is.na(pm_sp$r2), " N/A", sprintf("%6.3f", pm_sp$r2)),
      ifelse(is.na(pm_sp$mae), "N/A", sprintf("%5.3f", pm_sp$mae))
    ),
    
    sprintf(
      "  %-30s  r=%s  R²=%s  MAE=%s",
      "Full INLA:",
      ifelse(is.na(pm_full$r), " N/A", sprintf("%6.3f", pm_full$r)),
      ifelse(is.na(pm_full$r2), " N/A", sprintf("%6.3f", pm_full$r2)),
      ifelse(is.na(pm_full$mae), "N/A", sprintf("%5.3f", pm_full$mae))
    ),
    
    sprintf(
      "  %-30s  r=%s  R²=%s  MAE=%s",
      "Beta-binomial:",
      ifelse(is.na(pm_bb$r), " N/A", sprintf("%6.3f", pm_bb$r)),
      ifelse(is.na(pm_bb$r2), " N/A", sprintf("%6.3f", pm_bb$r2)),
      ifelse(is.na(pm_bb$mae), "N/A", sprintf("%5.3f", pm_bb$mae))
    )
  )
  
  # Model fit assessment
  summary_lines <- c(summary_lines, "", "----- MODEL FIT ASSESSMENT -----")
  
  all_r2 <- c(
    GLM = pm_glm$r2,
    Spatial = pm_sp$r2,
    Full = pm_full$r2,
    BetaBin = pm_bb$r2
  )
  
  best_r2_name <- names(which.max(all_r2[!is.na(all_r2)]))
  best_r2_val  <- max(all_r2, na.rm = TRUE)
  
  if (!is.na(best_r2_val)) {
    
    r2_assessment <- case_when(
      best_r2_val >= 0.7 ~ "Excellent (>70% variance explained)",
      best_r2_val >= 0.5 ~ "Good (50-70%)",
      best_r2_val >= 0.3 ~ "Moderate (30-50%)",
      best_r2_val >= 0.1 ~ "Weak (10-30%)",
      TRUE ~ "Poor (<10%)"
    )
    
    summary_lines <- c(
      summary_lines,
      paste0(
        "Best R² = ", round(best_r2_val, 3),
        " (", best_r2_name, ")  ",
        r2_assessment
      )
    )
  }
  
  if (!is.na(pm_glm$r2) && !is.na(pm_full$r2)) {
    summary_lines <- c(
      summary_lines,
      paste0(
        "Spatial field adds R² = ",
        round(pm_full$r2 - pm_glm$r2, 3)
      )
    )
  }
  
  if (!is.null(inla_fu) && !is.null(inla_bb)) {
    
    if (mc$bb_waic < mc$full_waic) {
      
      summary_lines <- c(
        summary_lines,
        "Beta-binomial outperforms binomial: overdispersion present"
      )
      
    } else {
      
      summary_lines <- c(
        summary_lines,
        "Binomial adequate - no strong overdispersion"
      )
    }
  }
  
  # Render summary text page
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    x = 0.05, y = 0.95, width = 0.9, height = 0.9,
    just = c("left", "top")
  ))
  grid::grid.text(
    label = paste(summary_lines, collapse = "\n"),
    x = 0, y = 1, just = c("left", "top"),
    gp = grid::gpar(fontfamily = "mono", fontsize = 7.5, lineheight = 1.15)
  )
  grid::popViewport()
  
  
  # ----- PAGE 3: Forest plot of fixed effects (odds ratios) --------------------------------------------------
  if (!is.null(feff) && nrow(feff %>% filter(term != "intercept")) > 0) {
    fe_plot <- feff %>%
      filter(term != "intercept") %>%
      mutate(term = factor(term, levels = rev(term))) %>%
      ggplot(aes(x = OR, y = term)) +
      geom_vline(xintercept = 1, linetype = "dashed", colour = "red") +
      geom_pointrange(aes(xmin = OR_lo, xmax = OR_hi,
                          colour = sig == "*"), size = 0.5) +
      scale_colour_manual(values = c("FALSE" = "grey50", "TRUE" = "steelblue"),
                          labels = c("Not significant", "Significant"),
                          name = "95% CrI excludes 1") +
      labs(title = paste(cname, " INLA Fixed Effects (Odds Ratios)"),
           subtitle = "Posterior mean and 95% credible interval",
           x = "Odds Ratio (exp(Î²))", y = NULL) +
      theme_minimal(base_size = 10) + theme(legend.position = "bottom")
    print(fe_plot)
  } else if (!is.null(glm_co) && nrow(glm_co %>% filter(term != "(Intercept)")) > 0) {
    fe_plot <- glm_co %>%
      filter(term != "(Intercept)") %>%
      mutate(term = factor(term, levels = rev(term))) %>%
      ggplot(aes(x = OR, y = term)) +
      geom_vline(xintercept = 1, linetype = "dashed", colour = "red") +
      geom_pointrange(aes(xmin = OR_lo, xmax = OR_hi,
                          colour = `Pr(>|z|)` < 0.05), size = 0.5) +
      scale_colour_manual(values = c("FALSE" = "grey50", "TRUE" = "steelblue"),
                          labels = c("p >= 0.05", "p < 0.05"),
                          name = "Significance") +
      labs(title = paste(cname, " GLM Fixed Effects (Odds Ratios)"),
           subtitle = "MLE estimate and 95% CI",
           x = "Odds Ratio (exp(Î²))", y = NULL) +
      theme_minimal(base_size = 10) + theme(legend.position = "bottom")
    print(fe_plot)
  }
  
  
  # ----- PAGE 4: Variable selection audit table --------------------------------------------------------------------------------
  audit_display <- var_audit %>%
    dplyr::select(variable, missingness_pct, step1_status,
                  step2_status, step3_status) %>%
    mutate(final = ifelse(grepl("FINAL MODEL", step3_status), "“", ""))
  
  audit_grob <- tryCatch(
    gridExtra::tableGrob(
      audit_display, rows = NULL,
      theme = gridExtra::ttheme_minimal(
        core    = list(fg_params = list(fontsize = 6)),
        colhead = list(fg_params = list(fontsize = 7, fontface = "bold"))
      )
    ),
    error = function(e) NULL
  )
  if (!is.null(audit_grob)) {
    grid::grid.newpage()
    grid::grid.text(
      paste(cname, " Variable Selection Audit Trail"),
      x = 0.5, y = 0.97, gp = grid::gpar(fontsize = 12, fontface = "bold")
    )
    grid::grid.draw(audit_grob)
  }
  
  
  # ----- PAGE 5: 6-panel diagnostic page -----------------------------------------------------------------------------------------------
  # Panel 1: Observed coverage map
  p1 <- ggplot()
  if (!is.null(adm1_sf))
    p1 <- p1 + geom_sf(data = adm1_sf, fill = NA, colour = "grey60", linewidth = 0.3)
  p1 <- p1 +
    geom_sf(data = sub_sf, aes(colour = coverage), size = 1.2, alpha = 0.7) +
    scale_colour_viridis_c(name = "Coverage", limits = c(0, 1)) +
    labs(title = paste(cname, " Observed MCV1 Coverage")) +
    theme_minimal(base_size = 10) + theme(legend.position = "bottom")
  
  # Panel 2: Predicted coverage map
  p2 <- ggplot()
  if (!is.null(adm1_sf))
    p2 <- p2 + geom_sf(data = adm1_sf, fill = NA, colour = "grey60", linewidth = 0.3)
  p2 <- p2 +
    geom_sf(data = sub_sf %>% mutate(pred_best = sub_m$pred_best),
            aes(colour = pred_best), size = 1.2, alpha = 0.7) +
    scale_colour_viridis_c(name = "Predicted", limits = c(0, 1)) +
    labs(title = paste(cname, " Model-Predicted MCV1 Coverage"),
         subtitle = bl) +
    theme_minimal(base_size = 10) + theme(legend.position = "bottom")
  
  # Panel 3: Residuals map (observed ˆ’ predicted)
  rl <- max(abs(sub_m$resid_best), na.rm = TRUE)
  p3 <- ggplot()
  if (!is.null(adm1_sf))
    p3 <- p3 + geom_sf(data = adm1_sf, fill = NA, colour = "grey60", linewidth = 0.3)
  p3 <- p3 +
    geom_sf(data = sub_sf %>% mutate(resid_best = sub_m$resid_best),
            aes(colour = resid_best), size = 1.2, alpha = 0.7) +
    scale_colour_gradient2(name = "Residual",
                           low = "#d73027", mid = "white", high = "#4575b4",
                           midpoint = 0, limits = c(-rl, rl)) +
    labs(title = paste(cname, " Residuals (Observed ˆ’ Predicted)")) +
    theme_minimal(base_size = 10) + theme(legend.position = "bottom")
  
  # Panel 4: Observed vs predicted scatter
  p4 <- ggplot(sub_m, aes(x = pred_best, y = coverage)) +
    geom_point(aes(size = n_children), alpha = 0.3, colour = "steelblue") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red") +
    scale_size_continuous(name = "N children", range = c(0.5, 4)) +
    labs(x = "Predicted coverage", y = "Observed coverage",
         title = paste(cname, " Observed vs Predicted")) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal(base_size = 10)
  
  # Panel 5: Coverage histogram
  p5 <- ggplot(sub_m, aes(x = coverage)) +
    geom_histogram(aes(fill = after_stat(x)), bins = 50,
                   colour = "white", linewidth = 0.2) +
    scale_fill_viridis_c() +
    geom_vline(xintercept = c(0, 1), linetype = "dashed", colour = "red") +
    labs(x = "MCV1 Coverage", y = "Clusters",
         title = paste(cname, " Coverage Distribution"),
         subtitle = paste0("At 0: ", round(mean(sub_m$coverage == 0) * 100, 1),
                           "% | At 1: ", round(mean(sub_m$coverage == 1) * 100, 1), "%")) +
    theme_minimal(base_size = 10) + guides(fill = "none")
  
  # Panel 6: Residual histogram
  p6 <- ggplot(sub_m, aes(x = resid_best)) +
    geom_histogram(bins = 50, fill = "grey60", colour = "white") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
    labs(x = "Residual", y = "Count",
         title = paste(cname, " Residual Distribution"),
         subtitle = paste0("Mean: ", round(mean(sub_m$resid_best, na.rm = TRUE), 3),
                           " | SD: ", round(sd(sub_m$resid_best, na.rm = TRUE), 3))) +
    theme_minimal(base_size = 10)
  
  combined_page <- (p1 | p2) / (p3 | p4) / (p5 | p6) +
    plot_annotation(
      title    = paste("Geostatistical Model Diagnostics:", cname),
      subtitle = paste("Model:", bl, "| Covariates:", length(vars_final)),
      theme    = theme(plot.title = element_text(face = "bold", size = 16),
                       plot.subtitle = element_text(size = 11, colour = "grey40"))
    )
  print(combined_page)
  
  # PAGE 6: Full-country prediction, spatial field, and covariate surfaces
  surface_diag <- make_country_surface_diagnostics(
    country_name = cname, adm1_sf = adm1_sf, sub_m = sub_m,
    vars_final = vars_final, inla_fu = inla_fu, mesh = mesh,
    spde = spde, lcrs = lcrs
  )
  if (!is.null(surface_diag)) {
    print(surface_diag$surface_plot)
    
    for (ii in seq(1, length(surface_diag$covariate_plots), by = 6)) {
      print(
        wrap_plots(
          surface_diag$covariate_plots[ii:min(ii + 5, length(surface_diag$covariate_plots))],
          ncol = 3
        ) +
          plot_annotation(
            title = paste(cname, "- 5 km covariate surfaces used in the model")
          )
      )
    }
    
    covariate_histograms <- make_covariate_histograms(sub_m, vars_final)
    
    for (ii in seq(1, length(covariate_histograms), by = 6)) {
      print(
        wrap_plots(
          covariate_histograms[ii:min(ii + 5, length(covariate_histograms))],
          ncol = 3
        ) +
          plot_annotation(
            title = paste(cname, "- cluster-level covariate histograms")
          )
      )
    }
    
    
    partial_residual_plots <- make_partial_residual_plots(
      sub_m = sub_m,
      vars_final = vars_final,
      feff = feff,
      glm_co = glm_co
    )
    
    if (!is.null(partial_residual_plots) && length(partial_residual_plots) > 0) {
      for (ii in seq(1, length(partial_residual_plots), by = 6)) {
        print(
          wrap_plots(
            partial_residual_plots[ii:min(ii + 5, length(partial_residual_plots))],
            ncol = 3
          ) +
            plot_annotation(
              title = paste(cname, "- partial residual plots by covariate")
            )
        )
      }
    }
    
    
    
    
    
  }
  
  # PAGE 7+: Residuals against each covariate, to flag possible non-linearity
  residual_covariate_plots <- make_residual_covariate_plots(sub_m, vars_final, cname)
  if (!is.null(residual_covariate_plots)) {
    for (ii in seq(1, length(residual_covariate_plots), by = 4)) {
      print(wrap_plots(residual_covariate_plots[ii:min(ii + 3, length(residual_covariate_plots))], ncol = 2) +
              plot_annotation(title = paste(cname, "- residual diagnostics by covariate")))
    }
  }
  
  dev.off()
  all_pdf_paths[[cname]] <- pdf_path
  cat("\n  PDF saved to:", pdf_path, "\n")  
  
  # ══════════════════════════════════════════════════════════════════
  # ═══  STORE RESULTS FOR EXCEL OUTPUT                                         ═══
  # ══════════════════════════════════════════════════════════════════
  
  sheets <- list()
  
  sheets[["Model_Summary"]] <- mc
  
  sheets[["Variables_Used"]] <- tibble(
    variable     = vars_final,
    bivariate_p  = biv$p_value[match(vars_final, biv$variable)],
    bivariate_OR = biv$OR[match(vars_final, biv$variable)],
    bivariate_beta = biv$estimate[match(vars_final, biv$variable)]
  ) %>% mutate(across(where(is.numeric), ~ round(.x, 4)))
  
  if (!is.null(glm_co))
    sheets[["GLM_Coefficients"]] <- glm_co %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))
  
  if (!is.null(feff))
    sheets[["INLA_Fixed_Effects"]] <- feff %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))
  
  if (!is.na(rp))
    sheets[["Spatial_Parameters"]] <- tibble(
      parameter = c("Spatial range (km)", "Spatial marginal variance",
                    "DIC (spatial-only)", "WAIC (spatial-only)",
                    "DIC (full model)", "WAIC (full model)",
                    "DIC (beta-binomial)", "WAIC (beta-binomial)"),
      value = c(round(rp, 1), round(sp, 4),
                mc$sp_dic, mc$sp_waic,
                mc$full_dic, mc$full_waic,
                mc$bb_dic, mc$bb_waic)
    )
  
  sheets[["Bivariate_Screening"]] <- biv %>%
    mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
    arrange(p_value)
  
  sheets[["Cluster_Predictions"]] <- sub_m %>%
    dplyr::select(LONGNUM, LATNUM, n_children, n_vaccinated, coverage,
                  any_of(c("pred_glm", "resid_glm",
                           "pred_sp", "resid_sp",
                           "pred_full", "resid_full",
                           "pred_bb", "resid_bb",
                           "pred_best", "resid_best"))) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))
  
  sheets[["Prediction_Metrics"]] <- tibble(
    model = c("Model 1: GLM", "Model 2: Spatial-only",
              "Model 3: Full INLA", "Model 4: Beta-binomial"),
    correlation_r = c(pm_glm$r, pm_sp$r, pm_full$r, pm_bb$r),
    R_squared     = c(pm_glm$r2, pm_sp$r2, pm_full$r2, pm_bb$r2),
    MAE           = c(pm_glm$mae, pm_sp$mae, pm_full$mae, pm_bb$mae),
    DIC           = c(NA, mc$sp_dic, mc$full_dic, mc$bb_dic),
    WAIC          = c(NA, mc$sp_waic, mc$full_waic, mc$bb_waic),
    AIC           = c(mc$glm_aic, NA, NA, NA),
    AUC_binary    = c(mc$glm_auc, NA, NA, NA)
  )
  
  sheets[["Variable_Audit"]] <- var_audit %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))
  
  if (nrow(cor_drop_log) > 0)
    sheets[["Correlation_Drops"]] <- cor_drop_log
  if (nrow(vif_drop_log) > 0)
    sheets[["VIF_Drops"]] <- vif_drop_log
  
  sheets[["Coverage_By_Size"]] <- cov_by_size_country
  
  sheets[["Model_Equation"]] <- tibble(
    item = c("Model type", "Formula (R syntax)",
             "Equation (fitted)", "Intercept",
             paste0("Covariate: ", vars_final)),
    value = c(
      eq_label,
      paste(deparse(glm_form), collapse = " "),
      eq_full,
      if (!is.null(feff)) as.character(round(feff$mean[feff$term == "intercept"], 4))
      else if (!is.null(glm_co)) as.character(round(glm_co$Estimate[glm_co$term == "(Intercept)"], 4))
      else "N/A",
      if (!is.null(feff)) {
        sapply(vars_final, function(v) {
          row <- feff %>% filter(term == v)
          if (nrow(row) == 1)
            paste0("Î²=", round(row$mean, 4), " (OR=", round(row$OR, 4),
                   ", 95% CrI: ", round(row$OR_lo, 4), "-", round(row$OR_hi, 4),
                   ")", ifelse(row$sig == "*", " *", ""))
          else "not in INLA model"
        })
      } else if (!is.null(glm_co)) {
        sapply(vars_final, function(v) {
          row <- glm_co %>% filter(term == v)
          if (nrow(row) == 1)
            paste0("Î²=", round(row$Estimate, 4), " (OR=", round(row$OR, 4),
                   ", 95% CI: ", round(row$OR_lo, 4), "-", round(row$OR_hi, 4),
                   ") ", row$sig)
          else "not in GLM"
        })
      } else rep("N/A", length(vars_final))
    )
  )
  
  if (!is.null(surface_diag)) {
    sheets[["Spatial_Surface_Grid"]] <- surface_diag$grid %>%
      mutate(across(where(is.numeric), ~ round(.x, 4)))
  }
  
  country_xl_path <- file.path(output_folder, paste0(country_slug, "_geostatistical_results_unweighted", version, ".xlsx"))
  writexl::write_xlsx(sheets, country_xl_path)
  all_xl_paths[[cname]] <- country_xl_path
  cat("  Excel workbook saved to:", country_xl_path, "\n")
  
  all_sheets[[cname]] <- sheets
  all_results[[cname]] <- list(
    vars_used   = vars_final,
    n_clusters  = nrow(sub_m),
    glm_coefs   = glm_co,
    inla_fixed  = feff,
    range_km    = rp,
    spatial_var = sp,
    best_model  = bl
  )
  
  cat("\n Country", cname, "complete.\n")
  
} 

# ----- END COUNTRY LOOP -----

# Per-country PDFs are closed inside the country loop.


# ═════════════════════════════════════════════════════════════════
# SECTION 9: SAVE EXCEL WORKBOOK
# ═════════════════════════════════════════════════════════════════

cat("\n----- SECTION 9: Saving Excel workbook -----\n")

xl <- list()

# ----- README sheet ------------------------------------------------------------------------------------------------------------------------------------------------------”€
readme_df <- tibble(
  Section = c(
    "OVERVIEW", "OVERVIEW", "OVERVIEW", "OVERVIEW",
    "VARIABLE SELECTION", "VARIABLE SELECTION", "VARIABLE SELECTION",
    "VARIABLE SELECTION", "VARIABLE SELECTION", "VARIABLE SELECTION",
    "MODELS", "MODELS", "MODELS", "MODELS",
    "INTERPRETING ORs", "INTERPRETING ORs", "INTERPRETING ORs",
    "INTERPRETING SPATIAL", "INTERPRETING SPATIAL",
    "INTERPRETING R²", "INTERPRETING R²", "INTERPRETING R²",
    "INTERPRETING MAE", "INTERPRETING MAE",
    "REFERENCES", "REFERENCES", "REFERENCES", "REFERENCES"
  ),
  Item = c(
    "Purpose", "Unit of analysis", "Outcome", "Framework",
    "Step 1: Missingness", "Step 1: Zero variance",
    "Step 2: Bivariate screening", "Step 2: Threshold",
    "Step 3a: Pairwise correlation", "Step 3b: GVIF",
    "Model 1 (GLM)", "Model 2 (Spatial-only)", "Model 3 (Full)", "Model 4 (Beta-binomial)",
    "OR > 1", "OR < 1", "OR = 1",
    "Spatial range", "Spatial variance",
    "R² > 0.5", "R² 0.3-0.5", "R² < 0.3",
    "MAE < 0.10", "MAE > 0.15",
    "Utazi 2022", "Fuglstad 2019", "Lindgren 2011", "Hosmer 2013"
  ),
  Description = c(
    "MCV1 coverage modelling across SSA using DHS cluster-level data, following Utazi et al. (2018, 2020, 2022).",
    "DHS survey cluster (GPS-jittered primary sampling unit) with aggregated child-level vaccination data.",
    "MCV1 coverage: proportion of children 12-23 months who received measles first-dose, modelled as y/n per cluster.",
    "Bayesian geostatistical modelling via R-INLA with SPDE mesh for approximating Gaussian processes.",
    "Variables with >5% missing data excluded (Utazi 2022: missingness may be biased in unknown ways).",
    "Variables with zero variance dropped (constant within country, e.g. all clusters rural).",
    "Crude odds ratios from univariate binomial GLMs: logit(p_i) = Î± + Î² Ã— x_i.",
    "Variables with p < 0.2 retained (Hosmer & Lemeshow 2013: liberal threshold to avoid premature exclusion).",
    "Pairs with |r| > 0.8: drop variable with weaker bivariate association.",
    "GVIF^(1/(2*Df)) > 2 threshold (Utazi 2022, Fox & Monette 1992). Iteratively remove worst offender.",
    "Non-spatial baseline: logit(p_i) = X_i Î². Tells us WHICH factors matter.",
    "Spatial-only: intercept + Matern GP + iid nugget, no covariates. Tells us how much location alone explains.",
    "Full model (PRIMARY): X_i Î² + S(s_i) + Z_i. Combines covariates (why) with spatial field (where).",
    "Same as Model 3 but beta-binomial likelihood. Tests for extra-binomial variation from the 0/1 boundary problem.",
    "Higher values †’ MORE vaccination. OR=1.25 means 25% higher odds per 1-unit increase.",
    "Higher values †’ LESS vaccination. OR=0.80 means 20% lower odds per 1-unit increase.",
    "No association: 95% CI/CrI crosses 1.",
    "Distance (km) where spatial correlation drops to ~13%. Large = regional patterns; small = local patterns.",
    "Variance of spatial field on logit scale. High = covariates miss important geographic patterns.",
    "Good: model explains >50% of coverage variation.",
    "Moderate: model explains 30-50% of variation.",
    "Weak: model explains <30%  substantial unexplained variation remains.",
    "Good: predictions within 10 percentage points on average.",
    "Poor: predictions off by >15 percentage points  limited practical value.",
    "Utazi CE et al. (2022). PLOS Global Public Health.",
    "Fuglstad G-A et al. (2019). J Am Stat Assoc 114: 445-452.",
    "Lindgren F et al. (2011). J R Stat Soc B 73: 423-498.",
    "Hosmer DW & Lemeshow S (2013). Applied Logistic Regression. 3rd ed. Wiley."
  )
)

xl[["ReadMe"]] <- readme_df

# Summary across all countries
comp <- bind_rows(all_comp)
xl[["Summary_All"]] <- comp

# Factor proportions (from Section 6)
if (exists("factor_props_all") && nrow(factor_props_all) > 0)
  xl[["Factor_Proportions"]] <- factor_props_all

if (exists("overall_props"))
  xl[["Overall_Proportions"]] <- as.data.frame(t(overall_props)) %>%
  tibble::rownames_to_column("metric")

xl[["Coverage_By_Size"]] <- cov_by_size

# Per-country sheets
for (cn in names(all_sheets)) {
  cs <- substr(gsub("[^A-Za-z0-9]", "", cn), 1, 15)
  for (st in names(all_sheets[[cn]])) {
    sn <- substr(paste0(cs, "_", substr(st, 1, 14)), 1, 31)
    xl[[sn]] <- all_sheets[[cn]][[st]]
  }
}

xl_path <- file.path(output_folder, paste0("country_geostatistical_results_unweighted", version, ".xlsx"))
writexl::write_xlsx(xl, xl_path)
cat(" Excel workbook saved to:", xl_path, "\n")


# ═════════════════════════════════════════════════════════════════
# SECTION 10: CROSS-COUNTRY SYNTHESIS
# ═════════════════════════════════════════════════════════════════

cat("
════════════════════════════════════════════════════════
═══                                                                  ═══
═══   CROSS-COUNTRY SYNTHESIS                                       ═══
═══                                                                  ═══
════════════════════════════════════════════════════════
")

if (nrow(comp) > 0) {
  
  cat("\n=== MODEL COMPARISON ACROSS COUNTRIES ===\n")
  print(comp, n = Inf, width = Inf)
  
  # Most commonly selected covariates
  cat("\n\n=== MOST COMMONLY SELECTED COVARIATES ===\n")
  cat("(Variables surviving the full screening pipeline in most countries)\n\n")
  
  all_vars <- unlist(lapply(all_results, `[[`, "vars_used"))
  var_freq <- sort(table(all_vars), decreasing = TRUE)
  var_freq_df <- data.frame(
    variable      = names(var_freq),
    n_countries   = as.integer(var_freq),
    pct_countries = round(as.integer(var_freq) / length(all_results) * 100, 1)
  )
  print(var_freq_df, row.names = FALSE)
  
  # Spatial vs non-spatial
  cat("\n\n=== SPATIAL vs NON-SPATIAL ===\n")
  sp_comp <- comp %>%
    filter(!is.na(full_waic)) %>%
    dplyr::select(country, n_clusters, glm_aic, full_waic, range_km, spatial_var)
  if (nrow(sp_comp) > 0) print(sp_comp, n = Inf)
  
  # Beta-binomial vs binomial
  cat("\n\n=== BETA-BINOMIAL vs BINOMIAL ===\n")
  bb_comp <- comp %>%
    filter(!is.na(full_waic), !is.na(bb_waic)) %>%
    mutate(
      delta_waic = round(full_waic - bb_waic, 1),
      better     = ifelse(bb_waic < full_waic, "Beta-binomial", "Binomial")
    ) %>%
    dplyr::select(country, full_waic, bb_waic, delta_waic, better)
  if (nrow(bb_comp) > 0) print(bb_comp, n = Inf)
}


# ═════════════════════════════════════════════════════════════════
# SECTION 11: NARRATIVE SUMMARY
# ═════════════════════════════════════════════════════════════════

cat("
════════════════════════════════════════════════════════
═══   ANALYSIS COMPLETE                                              ═══
════════════════════════════════════════════════════════

WHAT WE DID (per country):

  1. ASSEMBLED ~70 candidate covariates from DHS (Table A of Utazi 2022)
     plus geospatial/environmental data and ACLED conflict indicators.

  2. SCREENED for missingness (>5% †’ drop), bivariate association
     (p < 0.2 †’ keep), and multicollinearity (|r| > 0.8, GVIF > 2 †’ drop).

  3. COMPUTED FACTOR PROPORTIONS to describe the study population
     (e.g. % with primary education, % unemployed, % in urban slums).

  4. FITTED FOUR MODELS per country:
     Model 1: GLM (covariates only  the 'why')
     Model 2: Spatial-only INLA (location only  the 'where')
     Model 3: Full INLA (covariates + spatial  the complete picture)
     Model 4: Beta-binomial (overdispersion check for the 0/1 problem)

  5. GENERATED MAPS showing observed, predicted, and residual coverage
     with ADM1 administrative boundaries.

  6. ADDRESSED THE 0/1 PROBLEM via binomial likelihood, empirical logit,
     iid cluster nugget, and beta-binomial comparison.

WHY THIS MATTERS:

  The NON-SPATIAL model tells us WHICH factors drive coverage disparities
  (education, wealth, healthcare access, media exposure, etc.)  useful
  for designing targeted interventions.

  The SPATIAL model tells us WHERE unexplained pockets of low coverage
  remain  useful for targeting resources to areas with unmeasured local
  barriers (conflict, facility quality, programme failure).

  Together, they identify BOTH at-risk population sub-groups AND
  at-risk geographic areas.

KEY REFERENCES:
  Utazi CE et al. (2022). PLOS Global Public Health.
  Utazi CE et al. (2020). Lancet Digit Health.
  Utazi CE et al. (2018). Vaccine.
  Dong TQ, Wakefield J (2021). Vaccine 39: 2557-2569.
  Diggle PJ, Giorgi E (2021). J R Soc Interface 18: 20210104.
  Lindgren F et al. (2011). J R Stat Soc B 73: 423-498.
  Fox J, Monette G (1992). J Am Stat Assoc 87: 178-183.
  Hosmer DW et al. (2013). Applied Logistic Regression. 3rd ed. Wiley.
")

cat("\n Analysis complete. Output files:\n")
cat("   Per-country PDFs:\n")
cat(paste("     ", unlist(all_pdf_paths)), sep = "\n")
cat("\n   Per-country Excel workbooks:\n")
cat(paste("     ", unlist(all_xl_paths)), sep = "\n")
cat("\n   Combined Excel workbook:\n")
cat("     ", xl_path, "\n")










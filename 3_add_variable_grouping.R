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
  "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/all_KR_surveys_with_gps_HR_IR_maps_GC_ACLED.rds"
)

# ================================
# Helper: safe column accessor
# Returns NA if column does not exist in the data frame
# ================================
# col <- function(df, name) {
#   if (name %in% names(df)) as.numeric(df[[name]]) else rep(NA_real_, nrow(df))
# }
# 
# # ================================
# # Helper: first non-NA value across a list of column vectors
# # Used for "take the first available" coalesce-style logic
# # ================================
# first_nonNA <- function(...) {
#   vars <- list(...)
#   out <- rep(NA_real_, length(vars[[1]]))
#   for (v in vars) {
#     na_idx <- is.na(out)
#     out[na_idx] <- v[na_idx]
#   }
#   out
# }

# ================================
# Main function: adds all grouped variables to a single survey data frame
# ================================
add_grouped_variables <- function(df) {
  
  
  # Define helpers locally to avoid base::col() conflict
  get_col <- function(df, name) {
    if (name %in% names(df)) as.numeric(df[[name]]) else rep(NA_real_, nrow(df))
  }
  
  first_nonNA <- function(...) {
    vars <- list(...)
    out <- rep(NA_real_, length(vars[[1]]))
    for (v in vars) {
      na_idx <- is.na(out)
      out[na_idx] <- v[na_idx]
    }
    out
  }
  
  
  # ------------------------------------------------------------------
  # 1. respondent_edu_level
  # Maps v106. Raw DHS codes: 0=No education, 1=Primary, 2=Secondary,
  # 3=Higher, 9=Missing.
  # RECODE (Utazi et al. Table A): No education -> 0, Primary -> 1,
  # Secondary or Higher -> 2 (the paper collapses Secondary + Higher).
  # ------------------------------------------------------------------
  df$respondent_edu_level <- {
    x <-get_col(df, "v106")
    case_when(
      x == 0           ~ 0,   # no education
      x == 1           ~ 1,   # primary
      x == 2           ~ 2,   # secondary
      x == 3           ~ 3,   # higher
      TRUE             ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 2. respondent_edu_attainment
  # Maps v149 directly (six-category recode derived from v106 + v107).
  # Codes: 0=No edu, 1=Inc primary, 2=Comp primary, 3=Inc secondary,
  # 4=Comp secondary, 5=Higher, 9=Missing.
  # ------------------------------------------------------------------
  df$respondent_edu_attainment <- {
    x <-get_col(df, "v149")
    ifelse(x %in% 0:5, x, NA_real_)
  }
  
  # ------------------------------------------------------------------
  # 3. respondent_edu_years
  # Uses v133 (single years of education); falls back to v108 in older
  # surveys. Values 97 (inconsistent), 98 (DK) and 99 (missing) are
  # set to NA. Valid values are grouped into seven bins.
  # Bins: 0, 1-3, 4-6, 7-9, 10-12, 13-15, 16+
  # ------------------------------------------------------------------
  df$respondent_edu_years <- {
    x <- first_nonNA(get_col(df, "v133"),get_col(df, "v108"))
    x[x %in% c(97, 98, 99)] <- NA
    case_when(
      is.na(x)  ~ NA_real_,
      x == 0    ~ 0,
      x <= 3    ~ 1,
      x <= 6    ~ 2,
      x <= 9    ~ 3,
      x <= 12   ~ 4,
      x <= 15   ~ 5,
      x >= 16   ~ 6
    )
  }
  
  # ------------------------------------------------------------------
  # 4. respondent_edu_grade
  # Highest grade/year completed at the level given in v106.
  # Tries v107 first, then country-specific equivalents s109, s105,
  # s114, s904 in order. Values 98 (DK) and 99 (missing) set to NA.
  # ------------------------------------------------------------------
  df$respondent_edu_grade <- {
    x <- first_nonNA(
      get_col(df, "v107"),get_col(df, "s109"),
      get_col(df, "s105"),get_col(df, "s114"),
      get_col(df, "s904")
    )
    x[x %in% c(98, 99)] <- NA
    x
  }
  
  # ------------------------------------------------------------------
  # 5. partner_edu_level
  # Husband/partner's education level. Tries v701 first, then country-
  # specific equivalents s703, s704, s805. Codes: 0=No edu, 1=Primary,
  # 2=Secondary, 3=Higher, 8=DK, 9=Missing. DK and Missing set to NA.
  # ------------------------------------------------------------------
  df$partner_edu_level <- {
    x <- first_nonNA(
      get_col(df, "v701"),get_col(df, "s703"),
      get_col(df, "s704"),get_col(df, "s805")
    )
    x[x %in% c(8, 9)] <- NA
    ifelse(x %in% 0:3, x, NA_real_)
  }
  
  # ------------------------------------------------------------------
  # 6. partner_edu_years
  # Husband/partner's single years of education. Tries v715 first, then
  # s705 Values 97, 98, 99 set to NA. Grouped into seven bins.
  # Bins: 0, 1-3, 4-6, 7-9, 10-12, 13-15, 16+
  # ------------------------------------------------------------------
  df$partner_edu_years <- {
    x <- first_nonNA(
      get_col(df, "v715"),get_col(df, "s705")
    )
    x[x %in% c(97, 98, 99)] <- NA
    case_when(
      is.na(x)  ~ NA_real_,
      x == 0    ~ 0,
      x <= 3    ~ 1,
      x <= 6    ~ 2,
      x <= 9    ~ 3,
      x <= 12   ~ 4,
      x <= 15   ~ 5,
      x >= 16   ~ 6
    )
  }
  
  # ------------------------------------------------------------------
  # 7. partner_edu_attainment
  # Husband/partner's six-category educational attainment (v729),
  # analogous to respondent_edu_attainment. Codes: 0-5 and 9=Missing.
  # ------------------------------------------------------------------
  df$partner_edu_attainment <- {
    x <-get_col(df, "v729")
    ifelse(x %in% 0:5, x, NA_real_)
  }
  
  # ------------------------------------------------------------------
  # 8. respondent_literacy
  # Literacy of the respondent. Tries v155 first, then v108 (older
  # phases). Codes: 0=Cannot read, 1=Parts of sentence, 2=Whole
  # sentence, 3=No card, 4=Blind/visually impaired. 9=Missing -> NA.
  # ------------------------------------------------------------------
  df$respondent_literacy <- {
    x <- first_nonNA(get_col(df, "v155"),get_col(df, "v108"))
    x[x == 9] <- NA
    ifelse(x %in% 0:4, x, NA_real_)
  }
  
  # ------------------------------------------------------------------
  # 9. wealth_index_quintile
  # Household wealth quintile. Tries v190 (IR file) first, then hv270
  # (HR file), v190a / hv270a (urban/rural versions), country-specific
  # variants swlthir, shwlthir, swlthiu, shwlthiu, s052, sh052.
  # Raw DHS codes: 1=Poorest, 2=Poorer, 3=Middle, 4=Richer, 5=Richest.
  # RECODE (Utazi et al. Table A):
  # 0 = Poorest/Poorer (raw 1, 2)
  # 1 = Middle          (raw 3)
  # 2 = Richer/Richest  (raw 4, 5)
  # ------------------------------------------------------------------
  df$wealth_index_quintile <- {
    x <- first_nonNA(
      get_col(df, "v190"), get_col(df, "hv270"),
      get_col(df, "v190a"),get_col(df, "hv270a"),
      get_col(df, "swlthir"),get_col(df, "shwlthir"),
      get_col(df, "swlthiu"),get_col(df, "shwlthiu"),
      get_col(df, "s052"), get_col(df, "sh052")
    )
    case_when(
      x %in% c(1, 2) ~ 0,
      x == 3         ~ 1,
      x %in% c(4, 5) ~ 2,
      TRUE           ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 10. wealth_index_score
  # Continuous wealth factor score stored as an integer with 5 implied
  # decimal places; divide by 100,000 to obtain the true z-score.
  # Tries v191, hv271, v191a, hv271a, s051, sh051 in order.
  # ------------------------------------------------------------------
  df$wealth_index_score <- {
    x <- first_nonNA(
      get_col(df, "v191"), get_col(df, "hv271"),
      get_col(df, "v191a"),get_col(df, "hv271a"),
      get_col(df, "s051"), get_col(df, "sh051")
    )
    x / 100000
  }
  
  # ------------------------------------------------------------------
  # 11. urban_rural
  # De facto type of place of residence. Tries v025 first, then v102,
  # then hv025. Raw DHS codes: 1=Urban, 2=Rural.
  # RECODE (Utazi et al. Table A): Rural -> 0, Urban -> 1.
  # ------------------------------------------------------------------
  df$urban_rural <- {
    x <- first_nonNA(
      get_col(df, "v025"),get_col(df, "v102"),
      get_col(df, "hv025")
    )
    case_when(
      x == 2 ~ 0,   # Rural
      x == 1 ~ 1,   # Urban
      TRUE   ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 12. region
  # Region (first administrative level) of residence. Tries v101 first,
  # then v024 (copy added for analysis), then hv024. Country-specific
  # numeric codes; labels vary by survey. Raw codes retained here.
  # NOTE (Utazi et al. Table A / Fig E): "Administrative level one
  # regions in each country, grouped together in some cases" - any
  # collapsing of small adm-1 regions into larger model strata is
  # country-specific and must be done downstream.
  # ------------------------------------------------------------------
  df$region <- {
    first_nonNA(
      get_col(df, "v101"),get_col(df, "v024"),
      get_col(df, "hv024")
    )
  }
  
  # ------------------------------------------------------------------
  # 13. respondent_age
  # Respondent's current age in completed years from v012. Values 98
  # (DK) and 99 (missing) are set to NA. Typically 15-49.
  # RECODE (Utazi et al. Table A): 15-19 -> 0, 20-29 -> 1,
  #                                30-39 -> 2, 40-49 -> 3
  # ------------------------------------------------------------------
  df$respondent_age <- {
    x <-get_col(df, "v012")
    x[x %in% c(98, 99)] <- NA
    case_when(
      is.na(x)         ~ NA_real_,
      x >= 15 & x < 20 ~ 0,
      x >= 20 & x < 30 ~ 1,
      x >= 30 & x < 40 ~ 2,
      x >= 40 & x < 50 ~ 3,
      TRUE             ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 14. respondent_age_group
  # Respondent's age in standard 5-year groups from v013.
  # Codes: 1=15-19, 2=20-24, 3=25-29, 4=30-34, 5=35-39,
  #        6=40-44, 7=45-49.
  # ------------------------------------------------------------------
  df$respondent_age_group <- {
    x <-get_col(df, "v013")
    ifelse(x %in% 1:7, x, NA_real_)
  }
  
  # ------------------------------------------------------------------
  # 15. cluster_psu
  # Survey cluster / primary sampling unit identifier. Tries v001 (IR
  # file) first, then hv001 (HR file), then a generic 'cluster' column
  # used in some country-specific files.
  # ------------------------------------------------------------------
  df$cluster_psu <- {
    first_nonNA(
      get_col(df, "v001"),get_col(df, "hv001"),
      get_col(df, "cluster")
    )
  }
  
  # ------------------------------------------------------------------
  # 16. country_phase
  # Country code and DHS survey phase identifier. Tries v000 (IR file)
  # first, then hv000 (HR file). Stored as a 3-character string
  # (e.g. "TZ7"). Kept as character, not numeric.
  # ------------------------------------------------------------------
  df$country_phase <- {
    if ("v000" %in% names(df))  as.character(df[["v000"]])
    else if ("hv000" %in% names(df)) as.character(df[["hv000"]])
    else rep(NA_character_, nrow(df))
  }
  
  # ------------------------------------------------------------------
  # 17. total_children_ever_born
  # Total live births ever reported by respondent (v201). Value 99
  # (missing) set to NA. Grouped into bins of 3:
  # 0, 1-3, 4-6, 7-9, 10-12, 13+
  # ------------------------------------------------------------------
  df$total_children_ever_born <- {
    x <-get_col(df, "v201")
    x[x == 99] <- NA
    case_when(
      is.na(x) ~ NA_real_,
      x == 0   ~ 0,
      x <= 3   ~ 1,
      x <= 6   ~ 2,
      x <= 9   ~ 3,
      x <= 12  ~ 4,
      x >= 13  ~ 5
    )
  }
  
  # ------------------------------------------------------------------
  # 18. number_living_children
  # Total number of children currently living at home, derived as the
  # sum of v202 (sons at home) and v203 (daughters at home). Grouped
  # into bins of 3: 0, 1-3, 4-6, 7-9, 10-12, 13+
  # ------------------------------------------------------------------
  df$number_living_children <- {
    sons  <-get_col(df, "v202")
    daus  <-get_col(df, "v203")
    x     <- sons + daus
    x[is.na(sons) & is.na(daus)] <- NA
    case_when(
      is.na(x) ~ NA_real_,
      x == 0   ~ 0,
      x <= 3   ~ 1,
      x <= 6   ~ 2,
      x <= 9   ~ 3,
      x <= 12  ~ 4,
      x >= 13  ~ 5
    )
  }
  
  # ------------------------------------------------------------------
  # 19. child_sex
  # Sex of child from birth history (b4). Raw DHS codes: 1=Male, 2=Female.
  # RECODE (Utazi et al. Table A): Male=0, Female=1.
  # All other values set to NA.
  # ------------------------------------------------------------------
  df$child_sex <- {
    x <-get_col(df, "b4")
    case_when(
      x == 1 ~ 0,   # Male
      x == 2 ~ 1,   # Female
      TRUE   ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 20. child_age_months
  # Child's current age in completed months. Preferred source is b19
  # (DHS-7 standard: int((v008a - b18) / 30.4375)). If b19 is absent,
  # calculates from interview date CMC (v008) minus birth date CMC (b3).
  # If b18 and v008a are both present they are used for precision.
  # Values 98 and 99 set to NA.
  # ------------------------------------------------------------------
  df$child_age_months <- {
    b19  <-get_col(df, "b19")
    b3   <-get_col(df, "b3")
    v008 <-get_col(df, "v008")
    b18  <-get_col(df, "b18")
    v008a <-get_col(df, "v008a")
    
    # Method 1: use b19 if valid
    age <- ifelse(!is.na(b19) & !(b19 %in% c(98, 99)), b19, NA_real_)
    
    # Method 2: CMC subtraction v008 - b3
    age_cmc <- ifelse(is.na(age) & !is.na(b3) & !is.na(v008),
                      as.integer(v008 - b3), NA_real_)
    age <- ifelse(is.na(age), age_cmc, age)
    
    # Method 3: DHS-7 exact formula using v008a and b18 (decimal months)
    age_exact <- ifelse(is.na(age) & !is.na(b18) & !is.na(v008a),
                        as.integer((v008a - b18) / 30.4375), NA_real_)
    age <- ifelse(is.na(age), age_exact, age)
    age
  }
  
  # ------------------------------------------------------------------
  # 21. birth_order
  # Birth order (parity) of the child in the respondent's birth history.
  # Tries bord first, then bord98 (country-specific variant). Values 98
  # (DK) and 99 (missing) set to NA.
  # RECODE (Utazi et al. Table A): 1-2 -> 0, 3-5 -> 1, >5 -> 2
  # ------------------------------------------------------------------
  df$birth_order <- {
    x <- first_nonNA(get_col(df, "bord"),get_col(df, "bord98"))
    x[x %in% c(98, 99)] <- NA
    case_when(
      is.na(x) ~ NA_real_,
      x <= 2   ~ 0,   # 1st-2nd birth
      x <= 5   ~ 1,   # 3rd-5th
      x >= 6   ~ 2    # 6+
    )
  }
  
  # ------------------------------------------------------------------
  # 22. birth_interval_preceding
  # Preceding birth interval in months (b11). b12 (succeeding interval)
  # is included as a fallback in some older file structures. First births
  # have value 0 (not applicable) and are set to NA. Values >=998 set
  # to NA. Grouped into 6-month bins: 9-15, 16-22, 23-29, 30-36,
  # 37-43, 44-50, 51+
  # ------------------------------------------------------------------
  df$birth_interval_preceding <- {
    x <- first_nonNA(get_col(df, "b11"),get_col(df, "b12"))
    x[x == 0 | x >= 998] <- NA
    case_when(
      is.na(x) ~ NA_real_,
      x <= 15  ~ 1,   # 9-15 months
      x <= 22  ~ 2,   # 16-22 months
      x <= 29  ~ 3,   # 23-29 months
      x <= 36  ~ 4,   # 30-36 months
      x <= 43  ~ 5,   # 37-43 months
      x <= 50  ~ 6,   # 44-50 months
      x >= 51  ~ 7    # 51+ months
    )
  }
  
  # ------------------------------------------------------------------
  # 23. anc_visits_number
  # Total number of ANC visits for the most recent birth. Tries m14
  # first Values 97, 98, 99 set to NA.
  # RECODE (Utazi et al. Table A): No ANC/DK -> 0, 1-3 -> 1, 4+ -> 2
  # ------------------------------------------------------------------
  df$anc_visits_number <- {
    x <- first_nonNA(get_col(df, "m14"))
    x[x %in% c(97, 98, 99)] <- NA
    case_when(
      is.na(x) ~ 0,    # missing/DK collapsed with no ANC per paper
      x == 0   ~ 0,
      x <= 3   ~ 1,
      x >= 4   ~ 2
    )
  }
  
  # ------------------------------------------------------------------
  # 24. anc_first_visit_timing
  # Months pregnant at the first ANC visit (m13). Country-specific
  # fallbacks: p408, s306, s411. Values 97, 98, 99 set to NA.
  # Grouped into bins: 0, 1-2, 3-4, 5-6, 7-8, 9+
  # ------------------------------------------------------------------
  df$anc_first_visit_timing <- {
    x <- first_nonNA(
      get_col(df, "m13"),get_col(df, "p408"),
      get_col(df, "s306"),get_col(df, "s411")
    )
    x[x %in% c(97, 98, 99)] <- NA
    case_when(
      is.na(x) ~ NA_real_,
      x == 0   ~ 0,
      x <= 2   ~ 1,
      x <= 4   ~ 2,
      x <= 6   ~ 3,
      x <= 8   ~ 4,
      x >= 9   ~ 5
    )
  }
  
  # ------------------------------------------------------------------
  # 25. place_of_delivery
  # Place where the most recent birth was delivered (m15). The first
  # digit of the two-digit code indicates sector:
  # 1x = Home, 2x = Government/public facility, 3x = Private facility.
  # Values 96 (Other) and 99 (Missing) set to NA.
  # New codes: 1=Home (10-19), 2=Govt facility (20-29),
  #            3=Private facility (30-39)
  # ------------------------------------------------------------------
  df$place_of_delivery <- {
    x <-get_col(df, "m15")
    x[x %in% c(96, 99)] <- NA
    case_when(
      is.na(x)         ~ NA_real_,
      x >= 10 & x < 20 ~ 1,   # Home
      x >= 20 & x < 30 ~ 2,   # Government/public facility
      x >= 30 & x < 40 ~ 3,   # Private facility
      TRUE             ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 26. anc_provider
  # Type of ANC provider for the most recent pregnancy. Derived from
  # binary indicator variables m2a-m2n. Priority hierarchy:
  # 1 = Skilled (any of m2a doctor, m2b nurse/midwife, m2c aux midwife,
  #     m2d, m2e, m2f trained TBA = 1)
  # 2 = Traditional birth attendant (m2g = 1)
  # 3 = Relative (m2h = 1)
  # 4 = No one (m2n = 1 OR all m2a-m2m = 0)
  # ------------------------------------------------------------------
  df$anc_provider <- {
    m2a <-get_col(df, "m2a"); m2b <-get_col(df, "m2b")
    m2c <-get_col(df, "m2c"); m2d <-get_col(df, "m2d")
    m2e <-get_col(df, "m2e"); m2f <-get_col(df, "m2f")
    m2g <-get_col(df, "m2g"); m2h <-get_col(df, "m2h")
    m2n <-get_col(df, "m2n")
    
    skilled     <- pmax(m2a, m2b, m2c, m2d, m2e, m2f, na.rm = TRUE)
    traditional <- m2g
    relative    <- m2h
    no_one      <- m2n
    
    # All-zero check across m2a-m2h (skilled + traditional + relative)
    all_provided <- cbind(m2a, m2b, m2c, m2d, m2e, m2f, m2g, m2h)
    all_zero     <- rowSums(all_provided, na.rm = TRUE) == 0 &
      rowSums(!is.na(all_provided)) > 0
    
    case_when(
      !is.na(skilled)     & skilled     == 1 ~ 1,
      !is.na(traditional) & traditional == 1 ~ 2,
      !is.na(relative)    & relative    == 1 ~ 3,
      (!is.na(no_one)     & no_one      == 1) | all_zero ~ 4,
      TRUE ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 27. child_vacc_measles_months
  # Age of child at MCV1 vaccination in months, calculated from the
  # card-recorded vaccination date vs child's birth date (b3 CMC).
  # Formula: 12 * (h9y - 1900) + h9m - b3
  # Country-specific fallbacks use s4mm/s4my, s104m1 date vars, mea2
  # date vars. Set to NA if h9 is 0 (not vaccinated) or 8/9 (DK/NA).
  # Also set to NA if h9y is missing (98/99) or implausible.
  # ------------------------------------------------------------------
  df$child_vacc_measles_months <- {
    h9  <-get_col(df, "h9")
    h9m <-get_col(df, "h9m")
    h9y <-get_col(df, "h9y")
    b3  <-get_col(df, "b3")
    
    # Country-specific fallback date variables
    s4mm  <-get_col(df, "s4mm");  s4my  <-get_col(df, "s4my")
    mea2m <-get_col(df, "mea2m"); mea2y <-get_col(df, "mea2y")
    
    # Calculate age at vaccination using card month/year vs birth CMC
    age_vacc <- ifelse(
      !is.na(h9y) & !(h9y %in% c(98, 99)) &
        !is.na(h9m) & !(h9m %in% c(98, 99)) &
        !is.na(b3),
      12 * (h9y - 1900) + h9m - b3,
      NA_real_
    )
    
    # Fallback: country-specific s4mm / s4my
    age_vacc <- ifelse(
      is.na(age_vacc) &
        !is.na(s4my) & !(s4my %in% c(98, 99)) &
        !is.na(s4mm) & !(s4mm %in% c(98, 99)) &
        !is.na(b3),
      12 * (s4my - 1900) + s4mm - b3,
      age_vacc
    )
    
    # Fallback: mea2m / mea2y (country-specific Measles-2 date vars)
    age_vacc <- ifelse(
      is.na(age_vacc) &
        !is.na(mea2y) & !(mea2y %in% c(98, 99)) &
        !is.na(mea2m) & !(mea2m %in% c(98, 99)) &
        !is.na(b3),
      12 * (mea2y - 1900) + mea2m - b3,
      age_vacc
    )
    
    # Set to NA if child was not vaccinated or h9 is DK/missing
    age_vacc[!is.na(h9) & h9 %in% c(0, 8, 9)] <- NA
    
    # Set to NA if calculated value is implausible (negative or >120 mo)
    age_vacc[!is.na(age_vacc) & (age_vacc < 0 | age_vacc > 120)] <- NA
    
    age_vacc
  }
  
  # ------------------------------------------------------------------
  # 28. child_vacc_measles
  # Whether a child aged 12-23 months received MCV1 and by what source.
  # Uses h9 (MCV1 status) and child_age_months (derived above).
  # Country-specific vars s4m, s104m1, s1508m1, mea2 are checked if h9
  # is absent. Codes:
  # 0 = Not vaccinated (aged 12-23 months)
  # 1 = Card dated (h9 = 1)
  # 2 = Mother's report only (h9 = 2)
  # 3 = Card marked (h9 = 3)
  # 4 = Any source (h9 = 1, 2, or 3)
  # NA = Child not aged 12-23 months or vaccination status unknown
  # ------------------------------------------------------------------
  df$child_vacc_measles <- {
    age <- df$child_age_months
    h9  <-get_col(df, "h9")
    
    # Country-specific fallbacks for h9
    s4m    <-get_col(df, "s4m")
    s104m1 <-get_col(df, "s104m1")
    s1508m1<-get_col(df, "s1508m1")
    mea2   <-get_col(df, "mea2")
    
    # If h9 is missing, try country-specific alternatives
    vacc_source <- h9
    vacc_source[is.na(vacc_source)] <- s4m[is.na(vacc_source)]
    vacc_source[is.na(vacc_source)] <- s104m1[is.na(vacc_source)]
    vacc_source[is.na(vacc_source)] <- s1508m1[is.na(vacc_source)]
    vacc_source[is.na(vacc_source)] <- mea2[is.na(vacc_source)]
    
    in_window <- !is.na(age) & age >= 12 & age <= 23
    
    case_when(
      !in_window                               ~ NA_real_,
      vacc_source %in% c(8, 9) | is.na(vacc_source) ~ NA_real_,
      vacc_source == 0                         ~ 0,
      vacc_source == 1                         ~ 1,
      vacc_source == 2                         ~ 2,
      vacc_source == 3                         ~ 3,
      TRUE                                     ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 29. child_vacc_measles_on_time
  # Whether MCV1 was received on time (at or before 12 months of age),
  # derived from child_vacc_measles_months (age at vaccination).
  # 1 = On time (vaccinated at <= 12 months), 0 = Not on time or NA.
  # Only set if the date is available on card (h9 = 1 or 3).
  # ------------------------------------------------------------------
  df$child_vacc_measles_on_time <- {
    vacc_months <- df$child_vacc_measles_months
    h9          <-get_col(df, "h9")
    
    # Only use card-confirmed dates for timing (h9 = 1 card dated,
    # h9 = 3 card marked with date — both have recorded date)
    card_dated <- !is.na(h9) & h9 %in% c(1, 3)
    
    ifelse(
      card_dated & !is.na(vacc_months) & vacc_months <= 12, 1,
      ifelse(
        card_dated & !is.na(vacc_months) & vacc_months > 12, 0,
        NA_real_
      )
    )
  }
  
  # ------------------------------------------------------------------
  # 30. child_vacc_measles_late
  # Whether MCV1 was received late (after 12 months of age), derived
  # from child_vacc_measles_months.
  # 1 = Late (vaccinated at > 12 months), 0 = On time or NA.
  # ------------------------------------------------------------------
  df$child_vacc_measles_late <- {
    vacc_months <- df$child_vacc_measles_months
    h9          <-get_col(df, "h9")
    
    card_dated <- !is.na(h9) & h9 %in% c(1, 3)
    
    ifelse(
      card_dated & !is.na(vacc_months) & vacc_months > 12, 1,
      ifelse(
        card_dated & !is.na(vacc_months) & vacc_months <= 12, 0,
        NA_real_
      )
    )
  }
  
  # ------------------------------------------------------------------
  # 31. child_health_card
  # Whether the child has a health/vaccination card. Checks h1 (card
  # seen by interviewer = 1, not seen but reported = 2), h1a (yes/no),
  # and country-specific sj466a.
  # 1 = Has card (any source), 0 = No card, NA = unknown
  # ------------------------------------------------------------------
  df$child_health_card <- {
    h1    <-get_col(df, "h1")
    h1a   <-get_col(df, "h1a")
    sj466a<-get_col(df, "sj466a")
    
    has_card <- ((!is.na(h1)     & h1 %in% c(1, 2)) |
                   (!is.na(h1a)    & h1a == 1)         |
                   (!is.na(sj466a) & sj466a == 1))
    
    no_card  <- ((!is.na(h1)  & h1 == 0) |
                   (!is.na(h1a) & h1a == 0))
    
    any_available <- !is.na(h1) | !is.na(h1a) | !is.na(sj466a)
    
    case_when(
      has_card                       ~ 1,
      no_card & !has_card            ~ 0,
      any_available & !has_card      ~ 0,
      TRUE                           ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 32. vacc_campaign
  # Whether the child received any vaccination through a campaign
  # (supplemental immunisation activity). Checks h35 (received vacc in
  # a campaign), and individual campaign variables h36a-h36f.
  # 1 = Received campaign vaccination, 0 = No, NA = not available
  # ------------------------------------------------------------------
  df$vacc_campaign <- {
    h35  <-get_col(df, "h35")
    h36a <-get_col(df, "h36a"); h36b <-get_col(df, "h36b")
    h36c <-get_col(df, "h36c"); h36d <-get_col(df, "h36d")
    h36e <-get_col(df, "h36e"); h36f <-get_col(df, "h36f")
    
    any_campaign <- (!is.na(h35)  & h35  == 1) |
      (!is.na(h36a) & h36a == 1) |
      (!is.na(h36b) & h36b == 1) |
      (!is.na(h36c) & h36c == 1) |
      (!is.na(h36d) & h36d == 1) |
      (!is.na(h36e) & h36e == 1) |
      (!is.na(h36f) & h36f == 1)
    
    any_no_campaign <- (!is.na(h35)  & h35  == 0) |
      (!is.na(h36a) & h36a == 0)
    
    any_available <- !is.na(h35)  | !is.na(h36a) | !is.na(h36b) |
      !is.na(h36c) | !is.na(h36d) | !is.na(h36e) |
      !is.na(h36f)
    
    case_when(
      any_campaign                          ~ 1,
      any_no_campaign & !any_campaign       ~ 0,
      any_available   & !any_campaign       ~ 0,
      TRUE                                  ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 33. child_vitamin_a
  # Whether the child received vitamin A supplementation. Checks h34,
  # h33, h40, h41a, h41b (binary/card-report codes), m54, s463a.
  # Values 1, 2, or 3 in the card/report variables count as received.
  # 8 (DK) is ignored. 1 = Received, 0 = Not received, NA = unknown.
  # ------------------------------------------------------------------
  df$child_vitamin_a <- {
    h34   <-get_col(df, "h34")
    h33   <-get_col(df, "h33")
    h40   <-get_col(df, "h40")
    h41a  <-get_col(df, "h41a")
    h41b  <-get_col(df, "h41b")
    m54   <-get_col(df, "m54")
    s463a <-get_col(df, "s463a")
    
    received <- (!is.na(h34)   & h34   == 1)            |
      (!is.na(h33)   & h33   %in% c(1, 2, 3)) |
      (!is.na(h40)   & h40   %in% c(1, 2, 3)) |
      (!is.na(h41a)  & h41a  %in% c(1, 2, 3)) |
      (!is.na(h41b)  & h41b  %in% c(1, 2, 3)) |
      (!is.na(m54)   & m54   %in% c(1, 2, 3)) |
      (!is.na(s463a) & s463a %in% c(1, 2, 3))
    
    not_received <- (!is.na(h34)   & h34   == 0) |
      (!is.na(h33)   & h33   == 0) |
      (!is.na(h40)   & h40   == 0) |
      (!is.na(h41a)  & h41a  == 0) |
      (!is.na(h41b)  & h41b  == 0) |
      (!is.na(m54)   & m54   == 0) |
      (!is.na(s463a) & s463a == 0)
    
    any_available <- !is.na(h34) | !is.na(h33) | !is.na(h40) |
      !is.na(h41a) | !is.na(h41b) | !is.na(m54) |
      !is.na(s463a)
    
    case_when(
      received                       ~ 1,
      not_received & !received       ~ 0,
      any_available & !received      ~ 0,
      TRUE                           ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 34. household_size
  # Total number of household members (de jure or de facto). Tries v136
  # (IR file total), hv009, hv012, hv013 in order. Values 99 set to NA.
  # RECODE (Utazi et al. Table A):
  # 0 = Large    (>= 9 members)
  # 1 = Medium   (5 - 8 members)
  # 2 = Small    (<= 4 members)
  # Note: paper's polarity has small=2 as the (presumed lowest-risk)
  # reference category in the multi-level non-vaccination model.
  # ------------------------------------------------------------------
  df$household_size <- {
    x <- first_nonNA(
      get_col(df, "v136"),get_col(df, "hv009"),
      get_col(df, "hv012"),get_col(df, "hv013")
    )
    x[x == 99] <- NA
    case_when(
      is.na(x) ~ NA_real_,
      x >= 9   ~ 0,   # large
      x >= 5   ~ 1,   # medium (5-8)
      x <= 4   ~ 2    # small
    )
  }
  
  # ------------------------------------------------------------------
  # 35. hh_head_sex
  # Sex of the household head. Tries v151 (IR file) then hv219 (HR
  # file). Raw DHS codes: 1=Male, 2=Female.
  # RECODE (Utazi et al. Table A): Female -> 0, Male -> 1.
  # All other values set to NA.
  # ------------------------------------------------------------------
  df$hh_head_sex <- {
    x <- first_nonNA(get_col(df, "v151"),get_col(df, "hv219"))
    case_when(
      x == 2 ~ 0,   # Female
      x == 1 ~ 1,   # Male
      TRUE   ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 36. hh_head_age
  # Age of household head in completed years. Tries v152 then hv220.
  # Values 98 (DK) and 99 (missing) set to NA. Grouped into 10-year
  # age bins:
  # 1=10-19, 2=20-29, 3=30-39, 4=40-49, 5=50-59,
  # 6=60-69, 7=70-79, 8=80+
  # ------------------------------------------------------------------
  df$hh_head_age <- {
    x <- first_nonNA(get_col(df, "v152"),get_col(df, "hv220"))
    x[x %in% c(98, 99)] <- NA
    case_when(
      is.na(x) ~ NA_real_,
      x < 20   ~ 1,
      x < 30   ~ 2,
      x < 40   ~ 3,
      x < 50   ~ 4,
      x < 60   ~ 5,
      x < 70   ~ 6,
      x < 80   ~ 7,
      x >= 80  ~ 8
    )
  }
  
  # ------------------------------------------------------------------
  # 37. relationship_to_hh_head
  # Respondent's relationship to the household head (v150). Retained
  # as-is with standard DHS codes. Values 99 set to NA.
  # ------------------------------------------------------------------
  df$relationship_to_hh_head <- {
    x <-get_col(df, "v150")
    x[x == 99] <- NA
    x
  }
  
  # ------------------------------------------------------------------
  # 38. hh_electricity
  # Whether the household has electricity. Checks v119 (IR), hv206 (HR)
  # and country-specific sh121a. 1=Yes, 0=No.
  # ------------------------------------------------------------------
  df$hh_electricity <- {
    v119   <-get_col(df, "v119")
    hv206  <-get_col(df, "hv206")
    sh121a <-get_col(df, "sh121a")
    
    has_it  <- (!is.na(v119)   & v119   == 1) |
      (!is.na(hv206)  & hv206  == 1) |
      (!is.na(sh121a) & sh121a == 1)
    
    any_available <- !is.na(v119) | !is.na(hv206) | !is.na(sh121a)
    
    case_when(
      has_it                          ~ 1,
      any_available & !has_it         ~ 0,
      TRUE                            ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 39. hh_radio
  # Whether the household owns a radio. Checks v120 (IR) and hv207 (HR).
  # 1=Yes, 0=No.
  # ------------------------------------------------------------------
  df$hh_radio <- {
    v120  <-get_col(df, "v120")
    hv207 <-get_col(df, "hv207")
    
    has_it        <- (!is.na(v120)  & v120  == 1) | (!is.na(hv207) & hv207 == 1)
    any_available <- !is.na(v120) | !is.na(hv207)
    
    case_when(
      has_it                          ~ 1,
      any_available & !has_it         ~ 0,
      TRUE                            ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 40. hh_television
  # Whether the household owns a television. Checks v121 (IR) and
  # hv208 (HR). 1=Yes, 0=No.
  # ------------------------------------------------------------------
  df$hh_television <- {
    v121  <-get_col(df, "v121")
    hv208 <-get_col(df, "hv208")
    
    has_it        <- (!is.na(v121)  & v121  == 1) | (!is.na(hv208) & hv208 == 1)
    any_available <- !is.na(v121) | !is.na(hv208)
    
    case_when(
      has_it                          ~ 1,
      any_available & !has_it         ~ 0,
      TRUE                            ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 41. hh_mobile_phone
  # Whether the household owns a mobile phone. Checks hv243a (HR),
  # v169a (IR), country-specific sh122b, sh122c. 1=Yes, 0=No.
  # ------------------------------------------------------------------
  df$hh_mobile_phone <- {
    hv243a <-get_col(df, "hv243a"); v169a  <-get_col(df, "v169a")
    sh122b <-get_col(df, "sh122b"); sh122c <-get_col(df, "sh122c")
    
    has_it <- (!is.na(hv243a) & hv243a == 1) |
      (!is.na(v169a)  & v169a  == 1) |
      (!is.na(sh122b) & sh122b == 1) |
      (!is.na(sh122c) & sh122c == 1)
    
    any_available <- !is.na(hv243a) | !is.na(v169a) |
      !is.na(sh122b) | !is.na(sh122c)
    
    case_when(
      has_it                          ~ 1,
      any_available & !has_it         ~ 0,
      TRUE                            ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 42. hh_bicycle
  # Whether the household owns a bicycle. Checks v123 (IR) and hv210
  # (HR). 1=Yes, 0=No.
  # ------------------------------------------------------------------
  df$hh_bicycle <- {
    v123  <-get_col(df, "v123")
    hv210 <-get_col(df, "hv210")
    
    has_it        <- (!is.na(v123)  & v123  == 1) | (!is.na(hv210) & hv210 == 1)
    any_available <- !is.na(v123) | !is.na(hv210)
    
    case_when(
      has_it                          ~ 1,
      any_available & !has_it         ~ 0,
      TRUE                            ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 43. hh_motorcycle
  # Whether the household owns a motorcycle/scooter. Checks v124 (IR)
  # and hv211 (HR). 1=Yes, 0=No.
  # ------------------------------------------------------------------
  df$hh_motorcycle <- {
    v124  <-get_col(df, "v124")
    hv211 <-get_col(df, "hv211")
    
    has_it        <- (!is.na(v124)  & v124  == 1) | (!is.na(hv211) & hv211 == 1)
    any_available <- !is.na(v124) | !is.na(hv211)
    
    case_when(
      has_it                          ~ 1,
      any_available & !has_it         ~ 0,
      TRUE                            ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 44. hh_car
  # Whether the household owns a car/truck. Checks v125 (IR) and hv212
  # (HR). 1=Yes, 0=No.
  # ------------------------------------------------------------------
  df$hh_car <- {
    v125  <-get_col(df, "v125")
    hv212 <-get_col(df, "hv212")
    
    has_it        <- (!is.na(v125)  & v125  == 1) | (!is.na(hv212) & hv212 == 1)
    any_available <- !is.na(v125) | !is.na(hv212)
    
    case_when(
      has_it                          ~ 1,
      any_available & !has_it         ~ 0,
      TRUE                            ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 45. hh_computer
  # Whether the household owns a computer/laptop. Checks hv243e (HR),
  # and country-specific sh111l, sh110n, sh117h, sh117i (note: typo
  # 'sh1117i' corrected to 'sh117i'). 1=Yes, 0=No.
  # ------------------------------------------------------------------
  df$hh_computer <- {
    hv243e <-get_col(df, "hv243e"); sh111l <-get_col(df, "sh111l")
    sh110n <-get_col(df, "sh110n"); sh117h <-get_col(df, "sh117h")
    sh117i <-get_col(df, "sh117i")
    
    has_it <- (!is.na(hv243e) & hv243e == 1) |
      (!is.na(sh111l) & sh111l == 1) |
      (!is.na(sh110n) & sh110n == 1) |
      (!is.na(sh117h) & sh117h == 1) |
      (!is.na(sh117i) & sh117i == 1)
    
    any_available <- !is.na(hv243e) | !is.na(sh111l) | !is.na(sh110n) |
      !is.na(sh117h) | !is.na(sh117i)
    
    case_when(
      has_it                          ~ 1,
      any_available & !has_it         ~ 0,
      TRUE                            ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 46. internet_use
  # Whether the respondent or household uses the internet. Checks v171a
  # (individual use frequency; 1=not last month, 2=monthly, 3=daily),
  # v171b, and household access variables sh121g, sh132h, sh110o, sh121l.
  # 1 = Uses internet (any frequency), 0 = Does not use, NA = not asked.
  # ------------------------------------------------------------------
  df$internet_use <- {
    v171a  <-get_col(df, "v171a");  v171b  <-get_col(df, "v171b")
    sh121g <-get_col(df, "sh121g"); sh132h <-get_col(df, "sh132h")
    sh110o <-get_col(df, "sh110o"); sh121l <-get_col(df, "sh121l")
    
    uses <- (!is.na(v171a)  & v171a  %in% c(1, 2, 3)) |
      (!is.na(v171b)  & v171b  %in% c(1, 2, 3)) |
      (!is.na(sh121g) & sh121g == 1) |
      (!is.na(sh132h) & sh132h == 1) |
      (!is.na(sh110o) & sh110o == 1) |
      (!is.na(sh121l) & sh121l == 1)
    
    does_not_use <- (!is.na(v171a)  & v171a  == 0) |
      (!is.na(v171b)  & v171b  == 0) |
      (!is.na(sh121g) & sh121g == 0)
    
    any_available <- !is.na(v171a) | !is.na(v171b) | !is.na(sh121g) |
      !is.na(sh132h) | !is.na(sh110o) | !is.na(sh121l)
    
    case_when(
      uses                              ~ 1,
      does_not_use & !uses              ~ 0,
      any_available & !uses             ~ 0,
      TRUE                              ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 47. media_newspaper
  # Frequency of reading newspapers. Harmonises v157 (DHS5+, 4-point
  # scale) with v109 (older phases, binary: read at least weekly).
  # Country-specific s107. Codes:
  # 0 = Not at all, 1 = Less than weekly, 2 = At least weekly,
  # 3 = Almost every day
  # ------------------------------------------------------------------
  df$media_newspaper <- {
    v157 <-get_col(df, "v157"); v109 <-get_col(df, "v109")
    s107 <-get_col(df, "s107")
    
    case_when(
      !is.na(v157) & v157 == 3 ~ 3,
      !is.na(v157) & v157 == 2 ~ 2,
      !is.na(v157) & v157 == 1 ~ 1,
      !is.na(v109) & v109 == 1 ~ 2,   # at least once a week
      !is.na(s107) & s107 == 1 ~ 2,
      !is.na(v157) & v157 == 0 ~ 0,
      !is.na(v109) & v109 == 0 ~ 0,
      !is.na(s107) & s107 == 0 ~ 0,
      TRUE                      ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 48. media_radio
  # Frequency of listening to radio. Harmonises v158 (4-point scale)
  # with v111 and v112 (older binary phases). Country-specific 
  # Codes: 0=Not at all, 1=Less than weekly, 2=At least weekly,
  # 3=Almost every day
  # ------------------------------------------------------------------
  df$media_radio <- {
    v158 <-get_col(df, "v158"); v111 <-get_col(df, "v111")
    v112 <-get_col(df, "v112")
    
    case_when(
      !is.na(v158) & v158 == 3 ~ 3,
      !is.na(v158) & v158 == 2 ~ 2,
      !is.na(v111) & v111 == 1 ~ 2,   # at least weekly (old binary)
      !is.na(v158) & v158 == 1 ~ 1,
      !is.na(v112) & v112 == 1 ~ 1,   # less than weekly (old binary)
      !is.na(v158) & v158 == 0 ~ 0,
      !is.na(v111) & v111 == 0 ~ 0,
      TRUE                      ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 49. media_television
  # Frequency of watching television. Harmonises v159 (4-point scale)
  # with v110 (older binary). Country-specific s117.
  # Codes: 0=Not at all, 1=Less than weekly, 2=At least weekly,
  # 3=Almost every day
  # ------------------------------------------------------------------
  df$media_television <- {
    v159 <-get_col(df, "v159"); v110 <-get_col(df, "v110")
    s117 <-get_col(df, "s117")
    
    case_when(
      !is.na(v159) & v159 == 3 ~ 3,
      !is.na(v159) & v159 == 2 ~ 2,
      !is.na(v110) & v110 == 1 ~ 2,   # at least weekly (old binary)
      !is.na(v159) & v159 == 1 ~ 1,
      !is.na(s117) & s117 == 1 ~ 2,
      !is.na(v159) & v159 == 0 ~ 0,
      !is.na(v110) & v110 == 0 ~ 0,
      !is.na(s117) & s117 == 0 ~ 0,
      TRUE                      ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 50. health_insurance
  # Whether the respondent is covered by any health insurance. Checks
  # v481 (any coverage) and individual type variables v481a-v481d, v481x
  # (mutual/community, employer, social security, private, other).
  # 1 = Covered, 0 = Not covered, NA = not asked.
  # ------------------------------------------------------------------
  df$health_insurance <- {
    v481  <-get_col(df, "v481")
    v481a <-get_col(df, "v481a"); v481b <-get_col(df, "v481b")
    v481c <-get_col(df, "v481c"); v481d <-get_col(df, "v481d")
    v481x <-get_col(df, "v481x")
    
    covered <- (!is.na(v481)  & v481  == 1) |
      (!is.na(v481a) & v481a == 1) |
      (!is.na(v481b) & v481b == 1) |
      (!is.na(v481c) & v481c == 1) |
      (!is.na(v481d) & v481d == 1) |
      (!is.na(v481x) & v481x == 1)
    
    not_covered <- (!is.na(v481) & v481 == 0) |
      (rowSums(cbind(v481a, v481b, v481c, v481d, v481x),
               na.rm = TRUE) == 0 &
         rowSums(!is.na(cbind(v481a, v481b, v481c, v481d, v481x))) > 0)
    
    any_available <- !is.na(v481) | !is.na(v481a) | !is.na(v481b) |
      !is.na(v481c) | !is.na(v481d) | !is.na(v481x)
    
    case_when(
      covered                         ~ 1,
      not_covered & !covered          ~ 0,
      any_available & !covered        ~ 0,
      TRUE                            ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 51. respondent_working
  # Whether the respondent is currently working. Checks v714 (IR
  # standard) and country-specific s707.
  # 1 = Currently working, 0 = Not working, NA = not asked/missing.
  # ------------------------------------------------------------------
  df$respondent_working <- {
    v714 <-get_col(df, "v714"); s707 <-get_col(df, "s707")
    
    currently_working <- (!is.na(v714) & v714 == 1) |
      (!is.na(s707) & s707 == 1)
    
    not_working <- (!is.na(v714) & v714 == 0) |
      (!is.na(s707) & s707 == 0)
    
    any_available <- !is.na(v714) | !is.na(s707)
    
    case_when(
      currently_working               ~ 1,
      not_working & !currently_working~ 0,
      any_available                   ~ 0,
      TRUE                            ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 52. respondent_working_12mo
  # Whether the respondent worked at any point in the past 12 months
  # (even if not currently). Derived from v731.
  # 0 = Did not work in past 12 months, 1 = Worked in past 12 months
  # (including currently working), NA = not asked.
  # ------------------------------------------------------------------
  df$respondent_working_12mo <- {
    v731 <-get_col(df, "v731")
    ifelse(is.na(v731), NA_real_, ifelse(v731 %in% c(1, 2, 3), 1, 0))
  }
  
  # ------------------------------------------------------------------
  # 53. respondent_occupation
  # Respondent's occupation type. Tries v716 (grouped recode) first,
  # then v717 (raw), then sh57. Raw DHS codes:
  # 0=Not working, 1=Professional/technical/managerial, 2=Clerical,
  # 3=Sales, 4=Agric self-employed, 5=Agric employee,
  # 6=Household/domestic, 7=Services, 8=Skilled manual, 9=Unskilled
  # RECODE (Utazi et al. Table A):
  # 0 = Agricultural (raw 4, 5)
  # 1 = Clerical / sales / services (raw 2, 3, 7)
  # 2 = Professional / technical / managerial (raw 1)
  # 3 = Skilled / unskilled manual / other (raw 6, 8, 9)
  # NA = Not working / DK / missing (raw 0, 98, 99)
  # ------------------------------------------------------------------
  df$respondent_occupation <- {
    x <- first_nonNA(
      get_col(df, "v716"),get_col(df, "v717")
    )
    x[x %in% c(98, 99)] <- NA
    case_when(
      is.na(x)        ~ NA_real_,
      x %in% c(4, 5)  ~ 0,    # agricultural
      x %in% c(2, 3, 7) ~ 1,  # clerical/sales/services
      x == 1          ~ 2,    # prof/tech/managerial
      x %in% c(6, 8, 9) ~ 3,  # skilled/unskilled manual/domestic/other
      TRUE            ~ NA_real_   # x == 0 (not working) -> NA
    )
  }
  
  # ------------------------------------------------------------------
  # 54. partner_occupation
  # Partner's/husband's occupation. Tries v704 (grouped) first, then
  # v705 (raw), then country-specific s706, s908, s908a. Values 98
  # (DK) and 99 (missing/NA) set to NA.
  # ------------------------------------------------------------------
  df$partner_occupation <- {
    x <- first_nonNA(
      get_col(df, "v704"),get_col(df, "v705"),
      get_col(df, "s706"),
      get_col(df, "s908a")
    )
    x[x %in% c(98, 99)] <- NA
    x
  }
  
  # ------------------------------------------------------------------
  # 55. womens_autonomy_no_say
  # Count of decision-making domains (v743a-v743f) where the respondent
  # has NO say (code 4=husband alone, 5=someone else). Maximum score=6.
  # NA if none of v743a-v743f are present in the survey.
  # ------------------------------------------------------------------
  df$womens_autonomy_no_say <- {
    vars <- c("v743a","v743b","v743c","v743d","v743e","v743f")
    present <- vars[vars %in% names(df)]
    if (length(present) == 0) {
      rep(NA_real_, nrow(df))
    } else {
      mat <- sapply(present, function(v) {
        x <- get_col(df, v)
        ifelse(x %in% c(4, 5), 1, ifelse(x %in% c(1, 2, 3), 0, NA_real_))
      })
      rowSums(mat, na.rm = FALSE)
    }
  }
  
  # ------------------------------------------------------------------
  # 56. womens_autonomy_has_a_say
  # Count of decision-making domains (v743a-v743f) where the respondent
  # HAS a say (code 1=alone, 2=jointly, 3=with someone else).
  # Maximum score=6. NA if none of v743a-v743f present.
  # ------------------------------------------------------------------
  df$womens_autonomy_has_a_say <- {
    vars <- c("v743a","v743b","v743c","v743d","v743e","v743f")
    present <- vars[vars %in% names(df)]
    if (length(present) == 0) {
      rep(NA_real_, nrow(df))
    } else {
      mat <- sapply(present, function(v) {
        x <- get_col(df, v)
        ifelse(x %in% c(1, 2, 3), 1, ifelse(x %in% c(4, 5), 0, NA_real_))
      })
      rowSums(mat, na.rm = FALSE)
    }
  }
  
  # ------------------------------------------------------------------
  # 57. healthcare_barriers
  # Whether the respondent faces ANY barrier ("big problem") to
  # accessing healthcare. Checks all v467a-v467i (0=not a problem,
  # 1=big problem in DHS raw coding).
  # RECODE (Utazi et al. Table A): the polarity is INVERTED so that
  # the variable reads as "no problem seeking care = 1":
  #   0 = Had problem seeking medical advice/treatment
  #       (>=1 of the v467 items is a big problem)
  #   1 = Did not have a problem (none of v467a-i is a big problem)
  #   NA = no v467 variables present in this survey
  # ------------------------------------------------------------------
  df$healthcare_barriers <- {
    vars <- c("v467a","v467b","v467c","v467d","v467e","v467f","v467g","v467h","v467i")
    present <- vars[vars %in% names(df)]
    if (length(present) == 0) {
      rep(NA_real_, nrow(df))
    } else {
      mat <- sapply(present, function(v) get_col(df, v))
      any_barrier <- rowSums(mat == 1, na.rm = TRUE) > 0
      any_avail   <- rowSums(!is.na(mat)) > 0
      case_when(
        any_barrier              ~ 0,   # had problem
        any_avail & !any_barrier ~ 1,   # no problem
        TRUE                     ~ NA_real_
      )
    }
  }
  
  # ------------------------------------------------------------------
  # 58. healthcare_barriers_cost
  # Whether cost (getting money) is a barrier to healthcare access.
  # Directly from v467c: 0=Not a problem, 1=Big problem, NA=not asked.
  # ------------------------------------------------------------------
  df$healthcare_barriers_cost <- {
    x <-get_col(df, "v467c")
    case_when(
      is.na(x)  ~ NA_real_,
      x == 0    ~ 1, # no problem
      x %in% c(1, 2) ~ 0, # any problem
      TRUE      ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 59. healthcare_barriers_safety
  # Whether safety is a barrier (not wanting to go alone = v467f, or
  # no female provider = v467g). 1=Either is a big problem, 0=Neither.
  # ------------------------------------------------------------------
  
  df$healthcare_barriers_safety <- {
    x <-get_col(df, "v467f")
    case_when(
      is.na(x)  ~ NA_real_,
      x == 0    ~ 1, # no problem
      x %in% c(1, 2) ~ 0, # any problem
      TRUE      ~ NA_real_
    )
  }
  # 
  # df$healthcare_barriers_safety <- {
  #   v467f <-get_col(df, "v467f"); v467g <-get_col(df, "v467g")
  #   
  #   has_barrier <- (!is.na(v467f) & v467f == 1) |
  #     (!is.na(v467g) & v467g == 1)
  #   
  #   any_available <- !is.na(v467f) | !is.na(v467g)
  #   
  #   case_when(
  #     has_barrier                     ~ 1,
  #     any_available & !has_barrier    ~ 0,
  #     TRUE                            ~ NA_real_
  #   )
  # }
  # 
  # ------------------------------------------------------------------
  # 60. healthcare_barriers_transport
  # Whether transport/distance is a barrier (distance to facility =
  # v467d, having to take transport = v467e).
  # 1=Either is a big problem, 0=Neither.
  # NOTE: The spreadsheet listed v467h/v467i for this variable but those
  # are stockout variables; corrected here to use v467d/v467e per DHS
  # recode definitions.
  # ------------------------------------------------------------------
  
  
  
  df$healthcare_barriers_transport <- {
    x <-get_col(df, "v467d")
    case_when(
      is.na(x)  ~ NA_real_,
      x == 0    ~ 1, # no problem
      x %in% c(1, 2) ~ 0, # any problem
      TRUE      ~ NA_real_
    )
  }
  
  # df$healthcare_barriers_transport <- {
  #   v467d <-get_col(df, "v467d"); v467e <-get_col(df, "v467e")
  #   
  #   has_barrier <- (!is.na(v467d) & v467d == 1) |
  #     (!is.na(v467e) & v467e == 1)
  #   
  #   any_available <- !is.na(v467d) | !is.na(v467e)
  #   
  #   case_when(
  #     has_barrier                     ~ 1,
  #     any_available & !has_barrier    ~ 0,
  #     TRUE                            ~ NA_real_
  #   )
  # }
  
  # ------------------------------------------------------------------
  # 61. healthcare_barriers_stockout
  # Whether supply/provider stockout is a barrier (no provider available
  # = v467h, no drugs available = v467i).
  # 1=Either is a big problem, 0=Neither.
  # ------------------------------------------------------------------
  
  
  
  df$healthcare_barriers_stockout <- {
    x <-get_col(df, "v467h")
    case_when(
      is.na(x)  ~ NA_real_,
      x == 0    ~ 1, # no problem
      x %in% c(1, 2) ~ 0, # any problem
      TRUE      ~ NA_real_
    )
  }
  
  # df$healthcare_barriers_stockout <- {
  #   v467h <-get_col(df, "v467h"); v467i <-get_col(df, "v467i")
  #   
  #   has_barrier <- (!is.na(v467h) & v467h == 1) |
  #     (!is.na(v467i) & v467i == 1)
  #   
  #   any_available <- !is.na(v467h) | !is.na(v467i)
  #   
  #   case_when(
  #     has_barrier                     ~ 1,
  #     any_available & !has_barrier    ~ 0,
  #     TRUE                            ~ NA_real_
  #   )
  # }
  
  # ------------------------------------------------------------------
  # 62. religion
  # Respondent's religion (v130). Country-specific numeric codes.
  # Country-specific fallback: s118. Values 98 (DK) and 99 (missing)
  # set to NA. Raw codes retained here.
  # NOTE (Utazi et al. Table A): Final analysis groups religion into
  # categories based on the major religions practised in each country.
  # That collapsing requires per-country reference tables and must be
  # performed downstream (e.g. in the country-specific modelling step),
  # not in this generic harmonisation function.
  # ------------------------------------------------------------------
  df$religion <- {
    x <- first_nonNA(get_col(df, "v130"),get_col(df, "s118"))
    x[x %in% c(98, 99)] <- NA
    x
  }
  
  # ------------------------------------------------------------------
  # 63. ethnicity
  # Respondent's ethnicity (v131). Country-specific numeric codes.
  # Values 96 (Other), 98 (DK), 99 (missing) set to NA.
  # NOTE (Utazi et al. Table A): "Categories dependent on major ethnic
  # groups within each country" - any further collapsing must be done
  # per-country downstream.
  # ------------------------------------------------------------------
  df$ethnicity <- {
    x <-get_col(df, "v131")
    x[x %in% c(96, 98, 99)] <- NA
    x
  }
  
  # ------------------------------------------------------------------
  # 64. postnatal_care_newborn
  # Whether the newborn's health was checked before hospital discharge
  # (m74). 0=No, 1=Yes. NA if m74 not present.
  # ------------------------------------------------------------------
  df$postnatal_care_newborn <- {
    x <-get_col(df, "m74")
    ifelse(x %in% c(0, 1), x, NA_real_)
  }
  
  # ------------------------------------------------------------------
  # 65. postnatal_care_newborn_timing
  # Whether the baby received a postnatal check within 2 months (m70).
  # RECODE (Utazi et al. Table A): No/DK -> 0, Yes -> 1.
  # ------------------------------------------------------------------
  df$postnatal_care_newborn_timing <- {
    x <-get_col(df, "m70")
    case_when(
      x == 1                ~ 1,
      x == 0                ~ 0,
      is.na(x)              ~ 0,   # DK/missing collapsed with No per paper
      TRUE                  ~ 0
    )
  }
  
  # ------------------------------------------------------------------
  # 66. postnatal_care_newborn_person
  # Type of provider who conducted the postnatal check (m76). Codes
  # are country-specific but follow DHS provider type scheme. Retained
  # as-is; 98/99 set to NA.
  # ------------------------------------------------------------------
  df$postnatal_care_newborn_person <- {
    x <-get_col(df, "m76")
    x[x %in% c(98, 99)] <- NA
    x
  }
  
  # ------------------------------------------------------------------
  # 67. postnatal_care_newborn_first_person
  # Place/person of first postnatal contact for the newborn (m73).
  # Codes follow the DHS provider/place scheme (10=home etc.).
  # 98/99 set to NA.
  # NOTE: The spreadsheet description references m73 (place) rather
  # than m76 again. This variable uses m73.
  # ------------------------------------------------------------------
  df$postnatal_care_newborn_person <- {
    m76 <- get_col(df, "m76")
    
    # Handle missing values
    m76[m76 %in% c(98, 99)] <- NA
    
    skilled     <- m76 %in% c(10, 11, 12, 13)   # adjust if needed
    traditional_community <- m76 %in% c(21,22)

    case_when(
      skilled     ~ 1,
      traditional_community ~ 2,
      TRUE ~ NA_real_
    )
  }
  
  
  # ==================================================================
  # NEW VARIABLES (Utazi et al. 2022 framework — Table A)
  # ==================================================================
  
  # ------------------------------------------------------------------
  # 68. skilled_birth_attendance
  # Whether the most recent birth was attended by a skilled provider.
  # m3a doctor, m3b nurse/midwife, m3c auxiliary midwife, m3d/m3e/m3f
  # other trained personnel are considered "skilled". m3g (TBA),
  # m3h (relative/other), m3n (no one) are NOT skilled.
  # RECODE (Utazi et al. Table A): No skilled attendant -> 0,
  #                                 Skilled attendant   -> 1.
  # ------------------------------------------------------------------
  df$skilled_birth_attendance <- {
    m3a <- get_col(df, "m3a"); m3b <- get_col(df, "m3b")
    m3c <- get_col(df, "m3c"); m3d <- get_col(df, "m3d")
    m3e <- get_col(df, "m3e"); m3f <- get_col(df, "m3f")
    m3g <- get_col(df, "m3g"); m3h <- get_col(df, "m3h")
    m3n <- get_col(df, "m3n")
    
    skilled_mat <- cbind(m3a, m3b, m3c, m3d, m3e, m3f)
    unskilled_mat <- cbind(m3g, m3h, m3n)
    
    has_skilled  <- rowSums(skilled_mat == 1, na.rm = TRUE) > 0
    has_unskill  <- rowSums(unskilled_mat == 1, na.rm = TRUE) > 0
    any_avail    <- rowSums(!is.na(cbind(skilled_mat, unskilled_mat))) > 0
    
    case_when(
      has_skilled                ~ 1,
      any_avail & !has_skilled   ~ 0,
      TRUE                       ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 69. birth_quarter
  # Calendar quarter in which the index child was born, derived from
  # b1 (month of birth, 1-12).
  # RECODE (Utazi et al. Table A):
  # 0 = Jan-Mar, 1 = Apr-Jun, 2 = Jul-Sep, 3 = Oct-Dec
  # ------------------------------------------------------------------
  df$birth_quarter <- {
    m <- get_col(df, "b1")
    case_when(
      is.na(m)         ~ NA_real_,
      m %in% 1:3       ~ 0,
      m %in% 4:6       ~ 1,
      m %in% 7:9       ~ 2,
      m %in% 10:12     ~ 3,
      TRUE             ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 70. marital_status
  # Mother's current marital status from v501. Raw DHS codes:
  # 0 = Never in union, 1 = Married, 2 = Living with partner,
  # 3 = Widowed, 4 = Divorced, 5 = Not living together / separated.
  # RECODE (Utazi et al. Table A):
  # 0 = Never in union          (raw 0)
  # 1 = Married / living together (raw 1, 2)
  # 2 = Widowed/divorced/separated (raw 3, 4, 5)
  # ------------------------------------------------------------------
  df$marital_status <- {
    x <- get_col(df, "v501")
    case_when(
      x == 0           ~ 0,
      x %in% c(1, 2)   ~ 1,
      x %in% c(3, 4, 5) ~ 2,
      TRUE             ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 71. mother_owns_land
  # Whether the mother owns land, alone or jointly. From v745a:
  # 0 = does not own, 1 = alone, 2 = jointly only,
  # 3 = both alone and jointly.
  # RECODE (Utazi et al. Table A): Does not own -> 0,
  # Owns alone / jointly / both -> 1.
  # ------------------------------------------------------------------
  df$mother_owns_land_or_house <- {
    v745a <- get_col(df, "v745a")  # house
    v745b <- get_col(df, "v745b")  # land
    
    # Recode missing values if needed
    v745a[v745a %in% c(8, 9)] <- NA
    v745b[v745b %in% c(8, 9)] <- NA
    
    owns_house <- v745a %in% c(1, 2, 3)
    owns_land  <- v745b %in% c(1, 2, 3)
    
    case_when(
      owns_house | owns_land ~ 1,
      v745a == 0 & v745b == 0 ~ 0,
      TRUE ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 72. mother_bank_account
  # Whether the mother has a bank account (v170).
  # RECODE (Utazi et al. Table A): No -> 0, Yes -> 1.
  # ------------------------------------------------------------------
  df$mother_bank_account <- {
    x <- get_col(df, "v170")
    case_when(
      x == 0 ~ 0,
      x == 1 ~ 1,
      TRUE   ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 73. hh_mosquito_net
  # Whether the household owns at least one mosquito bednet.
  # Tries hv227 (have a bednet, DHS-V+) first, then hml1 (number of
  # mosquito nets owned, DHS-VI+, ownership = >=1).
  # RECODE (Utazi et al. Table A): No -> 0, Yes -> 1.
  # ------------------------------------------------------------------
  df$hh_mosquito_net <- {
    hv227 <- get_col(df, "hv227")
    hml1  <- get_col(df, "hml1")
    
    has_net <- (!is.na(hv227) & hv227 == 1) |
      (!is.na(hml1)  & hml1  >= 1)
    no_net  <- (!is.na(hv227) & hv227 == 0) |
      (!is.na(hml1)  & hml1  == 0)
    any_avail <- !is.na(hv227) | !is.na(hml1)
    
    case_when(
      has_net                  ~ 1,
      any_avail & !has_net     ~ 0,
      TRUE                     ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 74. mother_knows_malaria
  # Mother's knowledge of malaria, defined per Utazi et al. as
  # "agreed that malaria can lead to death OR that malaria can be cured".
  # These items live in country-specific malaria modules and use sml/sm
  # prefixes that vary across surveys. We search for any column whose
  # name matches typical malaria-knowledge patterns and code 1 if the
  # respondent gave an affirmative answer to ANY of them, 0 otherwise.
  # RECODE (Utazi et al. Table A): No -> 0, Yes -> 1.
  # NOTE: This is a generic best-effort match. For each country in the
  # final analysis you should verify the exact variables used and may
  # need to override this with country-specific logic.
  # ------------------------------------------------------------------
  df$mother_knows_malaria <- {
    candidate_patterns <- c(
      "malaria_death", "malaria_cure",
      "sm.*death", "sm.*cure",
      "sml.*death", "sml.*cure"
    )
    nm <- names(df)
    matches <- unique(unlist(lapply(candidate_patterns, function(p) {
      grep(p, nm, ignore.case = TRUE, value = TRUE)
    })))
    if (length(matches) == 0) {
      rep(NA_real_, nrow(df))
    } else {
      mat <- sapply(matches, function(v) get_col(df, v))
      knows <- rowSums(mat == 1, na.rm = TRUE) > 0
      any_avail <- rowSums(!is.na(mat)) > 0
      case_when(
        knows                ~ 1,
        any_avail & !knows   ~ 0,
        TRUE                 ~ NA_real_
      )
    }
  }
  
  # ------------------------------------------------------------------
  # 75. length_of_stay
  # Mother's length of continuous residence in current place (v104).
  # Raw DHS codes: years lived (1-94), 95 = always, 96 = visitor,
  # 98 = DK, 99 = missing.
  # RECODE (Utazi et al. Table A):
  # 0 = < 1 year / visitor    (raw 0, 96)
  # 1 = 1 - 3 years            (raw 1-3)
  # 2 = 4 - 5 years            (raw 4, 5)
  # 3 = > 5 years / always     (raw 6-94, 95)
  # ------------------------------------------------------------------
  df$length_of_stay <- {
    x <- get_col(df, "v104")
    x[x %in% c(98, 99)] <- NA
    case_when(
      is.na(x)         ~ NA_real_,
      #x == 96          ~ 0,    # visitor
      x == 0           ~ 0,    # less than 1 year
      x >= 1 & x <= 3  ~ 1,
      x >= 4 & x <= 5  ~ 2,
      x == 95          ~ 3,    # always
      x >= 6 & x <= 94 ~ 3,
      TRUE             ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 76a. hh_slum_dwelling (household-level slum flag)
  # Per the Utazi et al. supplementary methods (Table C and surrounding
  # text), households are scored 0/1 against four UN-Habitat criteria
  # and classified as a slum dwelling if their score is >= 2. The
  # cluster-level urban_slum variable (#76b) is then derived in a
  # second pass after this function returns, since aggregation needs
  # cluster IDs (created earlier in this function).
  #
  # Criteria and DHS variables (HR-file names tried first, IR-file
  # equivalents as fallback):
  # 1. UNimproved water:    hv201 / v113
  # 2. UNimproved sanitation: hv205 / v116, plus shared (hv225 / v160)
  # 3. NON-durable floor:   hv213 / v127
  # 4. Overcrowding:         hv009 (members) and hv216 (sleep rooms);
  #                          IR fallback uses v136 and v161
  #
  # Improved water (per JMP, used by DHS): codes 11-14, 21, 31, 41, 51,
  # 71 -> improved. Anything else -> unimproved.
  # Improved sanitation: 11-15, 21, 22, 41 -> improved AND not shared.
  # Durable floor: first digit 3 (finished, e.g. 30-37) -> durable.
  # Overcrowding: members per sleeping room > 3 -> overcrowded.
  # ------------------------------------------------------------------
  df$hh_slum_dwelling <- {
    # Water
    water <- first_nonNA(get_col(df, "hv201"), get_col(df, "v113"))
    improved_water_codes <- c(10, 11, 12, 13, 14, 21, 31, 41, 51, 71)
    unimproved_water <- !is.na(water) & !(water %in% improved_water_codes)
    
    # Sanitation
    toilet <- first_nonNA(get_col(df, "hv205"), get_col(df, "v116"))
    shared <- first_nonNA(get_col(df, "hv225"), get_col(df, "v160"))
    improved_toilet_codes <- c(10, 11, 12, 13, 14, 15, 21, 22, 41)
    is_improved_toilet <- !is.na(toilet) & toilet %in% improved_toilet_codes
    is_shared_toilet   <- !is.na(shared) & shared == 1
    unimproved_san <- !is.na(toilet) &
      (!is_improved_toilet | (is_improved_toilet & is_shared_toilet))
    
    # Durability of housing (floor material)
    floor <- first_nonNA(get_col(df, "hv213"), get_col(df, "v127"))
    # Finished/durable floors have first digit 3 (codes 30-39)
    nondurable_floor <- !is.na(floor) & !(floor >= 30 & floor < 40)
    
    # Overcrowding: > 3 people per sleeping room
    members <- first_nonNA(get_col(df, "hv009"), get_col(df, "hv012"),
                           get_col(df, "v136"))
    rooms   <- first_nonNA(get_col(df, "hv216"), get_col(df, "v161"))
    rooms[!is.na(rooms) & rooms == 0] <- 1   # avoid div-by-zero
    crowding_ratio <- members / rooms
    overcrowded <- !is.na(crowding_ratio) & crowding_ratio > 3
    
    # Score: how many of the four criteria are met (treats NA as 0
    # for the purposes of scoring, but flags rows with no criterion
    # data at all as NA)
    score_mat <- cbind(unimproved_water, unimproved_san,
                       nondurable_floor, overcrowded)
    score <- rowSums(score_mat, na.rm = TRUE)
    any_data <- !is.na(water) | !is.na(toilet) | !is.na(floor) |
      (!is.na(members) & !is.na(rooms))
    
    case_when(
      !any_data    ~ NA_real_,
      score >= 2   ~ 1,
      TRUE         ~ 0
    )
  }
  
  # ------------------------------------------------------------------
  # 47b. media_exposure  (combined indicator)
  # Whether the mother is exposed to ANY mass media at least weekly.
  # Built from media_radio, media_television, media_newspaper which
  # were defined above in #47-49 (codes 0-3, with >=2 = at least
  # weekly).
  # RECODE (Utazi et al. Table A): No -> 0,
  # Yes (radio OR TV OR newspaper at least weekly) -> 1.
  # ------------------------------------------------------------------
  df$media_exposure <- {
    radio <- df$media_radio
    tv    <- df$media_television
    paper <- df$media_newspaper
    
    weekly <- (!is.na(radio) & radio >= 2) |
      (!is.na(tv)    & tv    >= 2) |
      (!is.na(paper) & paper >= 2)
    any_avail <- !is.na(radio) | !is.na(tv) | !is.na(paper)
    
    case_when(
      weekly                 ~ 1,
      any_avail & !weekly    ~ 0,
      TRUE                   ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # 41b/46b. mobile_internet_use  (combined indicator)
  # Whether the mother uses a mobile phone OR the internet.
  # Built from internet_use (#46) and hh_mobile_phone (#41).
  # RECODE (Utazi et al. Table A): Never used phone or internet -> 0,
  # Used phone OR internet -> 1.
  # ------------------------------------------------------------------
  df$mobile_internet_use <- {
    phone <- df$hh_mobile_phone
    inet  <- df$internet_use
    
    uses <- (!is.na(phone) & phone == 1) |
      (!is.na(inet)  & inet  == 1)
    any_avail <- !is.na(phone) | !is.na(inet)
    
    case_when(
      uses                  ~ 1,
      any_avail & !uses     ~ 0,
      TRUE                  ~ NA_real_
    )
  }
  
  # ------------------------------------------------------------------
  # conflict_area  (admin-2 level, from upstream geospatial merge)
  # The conflict indicator must be attached at adm-2 level by the
  # earlier pipeline step that merged ACLED / similar conflict data
  # against admin-2 polygons. We search for that pre-existing column
  # by common naming patterns and recode whatever 0/1 indicator is
  # present so the polarity matches the paper.
  # RECODE (Utazi et al. Table A): Yes (conflict) -> 0, No -> 1.
  # If no upstream column is found, this is left as NA and you should
  # add the merge in an earlier pipeline stage.
  # ------------------------------------------------------------------
  df$conflict_area <- {
    candidates <- c("recent_violence_adm2")
    found <- candidates[candidates %in% names(df)]
    if (length(found) == 0) {
      rep(NA_real_, nrow(df))
    } else {
      x <- get_col(df, found[1])
      case_when(
        x >= 1 ~ 1,   # conflict present (1 or count > 0)
        x == 0 ~ 0,   # no conflict
        TRUE   ~ NA_real_
      )
    }
  }
  
  return(df)
}



# ================================
# Apply to all surveys in the list
# ================================

final_DHS_data_v2 <- purrr::map(final_DHS_data, function(survey) {
  if (!is.data.frame(survey)) return(survey)
  add_grouped_variables(survey)
}, .progress = TRUE)



# ================================
# Cluster-level urban_slum classification (Utazi et al., Table C)
# ================================
# A cluster is classified as a slum area if:
#   (i)  it is urban (urban_rural == 1),
#   (ii) it contains at least 10 households, and
#   (iii) at least 75% of those households are slum dwellings
#         (hh_slum_dwelling == 1).
# Classification proceeds per (country_phase, cluster_psu) within each
# survey, then is broadcast back to every child row.
# ------------------------------------------------------------------

final_DHS_data_v2 <- purrr::map(final_DHS_data_v2, function(survey) {
  
  if (!is.data.frame(survey)) return(survey)
  if (!all(c("cluster_psu", "urban_rural", "hh_slum_dwelling") %in%
           names(survey))) {
    survey$urban_slum <- NA_real_
    return(survey)
  }
  
  # Identify household ID column - try common DHS HR identifiers
  hh_id_col <- intersect(c("hhid", "hv002", "v002"), names(survey))[1]
  if (is.na(hh_id_col)) {
    # fall back to row-level (less accurate but better than nothing)
    survey$.hh_key <- paste(survey$cluster_psu, seq_len(nrow(survey)))
    hh_id_col <- ".hh_key"
  }
  
  # One row per household per cluster for the aggregation
  hh_lvl <- survey %>%
    dplyr::select(cluster_psu, !!rlang::sym(hh_id_col),
                  urban_rural, hh_slum_dwelling) %>%
    dplyr::distinct()
  
  cluster_summary <- hh_lvl %>%
    dplyr::group_by(cluster_psu) %>%
    dplyr::summarise(
      n_hh           = sum(!is.na(hh_slum_dwelling)),
      n_slum         = sum(hh_slum_dwelling == 1, na.rm = TRUE),
      pct_slum       = ifelse(n_hh > 0, n_slum / n_hh, NA_real_),
      cluster_urban  = any(urban_rural == 1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      urban_slum = dplyr::case_when(
        !cluster_urban             ~ 0,   # rural clusters: not slum
        n_hh < 10                  ~ NA_real_,
        pct_slum >= 0.75           ~ 1,   # slum (paper polarity: yes=0)
        TRUE                       ~ 0    # urban non-slum
      )
    ) %>%
    dplyr::select(cluster_psu, urban_slum)
  
  survey <- dplyr::left_join(survey, cluster_summary, by = "cluster_psu")
  if (".hh_key" %in% names(survey)) survey$.hh_key <- NULL
  survey
}, .progress = TRUE)



identical(sapply(final_DHS_data, nrow), sapply(final_DHS_data_v2, nrow))


# writexl::write_xlsx(head(final_DHS_data_v2[[6]], n=200), "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/final_DHS_data_v2_TEST.xlsx")

# ================================
# Save output
# ================================
saveRDS(
  final_DHS_data_v2,
  "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/final_DHS_data.rds"
)

message("Done. All surveys processed and saved.")



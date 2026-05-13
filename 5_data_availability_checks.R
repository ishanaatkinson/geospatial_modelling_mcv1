

final_DHS_data<- readRDS("C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/final_DHS_data.rds")



vars_to_check <- c(
  # ── Identifiers & survey structure ──────────────────────────────
  "v000", "v001", "v002", "v003", "caseid",
  "SurveyYear", "DHS_CountryCode", "CountryName",
  "cluster", "hv000", "hv001", "hv004", "hv021",
  
  # ── GPS / geography ──────────────────────────────────────────────
  "LATNUM", "LONGNUM", "DHSID", "DHSCC", "DHSYEAR",
  "ADM1NAME", "ADM1DHS", "ADM1FIPS",
  "travel_time_to_city",
  "travel_time_to_HC_motor",
  "travel_time_to_HC_walk",
  
  # ── Child identifiers ────────────────────────────────────────────
  "bidx", "bord", "b3", "b4", "b5",
  
  # ── Vaccination (source variables) ──────────────────────────────
  "h9", "h9d", "h9m", "h9y",
  "h9a", "h9ad", "h9am", "h9ay",
  "s4m", "s4md", "s4mm", "s4my",
  "s104m1", "s1508m1",
  "mea2", "mea2d", "mea2m", "mea2y",
  "h1", "h1a", "sj466a",
  "h35", "h36a", "h36b", "h36c", "h36d", "h36e", "h36f",
  
  # ── Child age (source variables) ─────────────────────────────────
  "b19", "b8", "hw1", "b18", "v008", "v008a",
  
  # ── Birth history ────────────────────────────────────────────────
  "bord", "bord98", "b11", "b12",
  
  # ── Education (source variables) ─────────────────────────────────
  "v106", "v107", "v149", "v133", "v108",
  "s109", "s105", "s114", "s904",
  "v155", "v701", "v715", "v729",
  "s703", "s704", "s805", "s705", "v207",
  
  # ── Wealth (source variables) ────────────────────────────────────
  "v190", "hv270", "v190a", "hv270a",
  "v191", "hv271", "v191a", "hv271a",
  "swlthir", "shwlthir", "swlthiu", "shwlthiu",
  "s051", "sh051", "s052", "sh052",
  
  # ── Geography / residence (source variables) ─────────────────────
  "v025", "v102", "hv025",
  "v101", "v024", "hv024",
  
  # ── ANC (source variables) ───────────────────────────────────────
  "m14", "p409", "m13", "p408", "s306", "s411",
  "m2a", "m2b", "m2c", "m2d", "m2e",
  "m2f", "m2g", "m2h", "m2n",
  
  # ── Delivery (source variables) ──────────────────────────────────
  "m15",
  
  # ── Household assets (source variables) ──────────────────────────
  "v119", "hv206", "sh121a",
  "v120", "hv207",
  "v121", "hv208",
  "v123", "hv210",
  "v124", "hv211",
  "v125", "hv212",
  "hv243a", "v169a", "sh122b", "sh122c",
  "hv243e", "sh111l", "sh110n", "sh117h", "sh117i",
  "v171a", "v171b", "sh121g", "sh132h", "sh110o", "sh121l",
  
  # ── Household composition ────────────────────────────────────────
  "v136", "hv009", "hv012", "hv013",
  "v151", "hv219", "v152", "hv220", "v150",
  "v202", "v203", "v201", "v218",
  
  # ── Media (source variables) ─────────────────────────────────────
  "v157", "v109", "s107",
  "v158", "v111", "v112", "s116",
  "v159", "v110", "s117",
  
  # ── Work & occupation (source variables) ─────────────────────────
  "v714", "s707", "v731",
  "v716", "v717", "sh57",
  "v704", "v705", "s706", "s908", "s908a",
  
  # ── Women's autonomy (source variables) ──────────────────────────
  "v743a", "v743b", "v743c", "v743d", "v743e", "v743f",
  
  # ── Healthcare barriers (source variables) ───────────────────────
  "v467a", "v467b", "v467c", "v467d",
  "v467e", "v467f", "v467g", "v467h", "v467i",
  
  # ── Health insurance (source variables) ──────────────────────────
  "v481", "v481a", "v481b", "v481c", "v481d", "v481x",
  
  # ── Religion & ethnicity (source variables) ──────────────────────
  "v130", "s118", "v131",
  
  # ── Postnatal care (source variables) ────────────────────────────
  "m70", "m73", "m74", "m76",
  
  # ── Vitamin A & campaign (source variables) ──────────────────────
  "h33", "h34", "h40", "h41a", "h41b", "m54", "s463a",
  
  
  # ---- Geospatial covariates of interest --------------
  
    "NIGHTLIGHTS_COMPOSITE",
    "RAINFALL_2000",
    "RAINFALL_2005",
    "RAINFALL_2010",
    "RAINFALL_2015",
    "U5_POPULATION_2000",
    "U5_POPULATION_2005",
    "U5_POPULATION_2010",
    "U5_POPULATION_2015",
    "UN_POPULATION_COUNT_2000",
    "UN_POPULATION_COUNT_2005",
    "UN_POPULATION_COUNT_2010",
    "UN_POPULATION_COUNT_2015",
    "UN_POPULATION_DENSITY_2000",
    "UN_POPULATION_DENSITY_2005",
    "UN_POPULATION_DENSITY_2010",
    "UN_POPULATION_DENSITY_2015",
    "WET_DAYS_2000",
    "WET_DAYS_2005",
    "WET_DAYS_2010",
    "WET_DAYS_2015",
    "MALARIA_INCIDENCE_2000",
    "MALARIA_INCIDENCE_2005",
    "MALARIA_INCIDENCE_2010",
    "MALARIA_INCIDENCE_2015",
    "MALARIA_PREVALENCE_2000",
    "MALARIA_PREVALENCE_2005",
    "MALARIA_PREVALENCE_2010",
    "MALARIA_PREVALENCE_2015",
    "TRAVEL_TIMES_2000",
    "TRAVEL_TIMES_2015",
    "PRECIPITATION_2000",
    "PRECIPITATION_2005",
    "PRECIPITATION_2010",
    "PRECIPITATION_2015",
    "PRECIPITATION_2020",
    "RAINFALL_2020",
    "TRAVEL_TIMES",
    "U5_POPULATION_2020",
    "UN_POPULATION_COUNT_2020",
    "UN_POPULATION_DENSITY_2020",
    "WET_DAYS_2020",
    "MALARIA_INCIDENCE_2020",
    "MALARIA_PREVALENCE_2020",
  
  # ── Derived grouped variables ─────────────────────────────────────
  "respondent_edu_level", "respondent_edu_attainment",
  "respondent_edu_years", "respondent_edu_grade",
  "respondent_literacy",
  "partner_edu_level", "partner_edu_years", "partner_edu_attainment",
  "wealth_index_quintile", "wealth_index_score",
  "urban_rural", "region",
  "respondent_age", "respondent_age_group",
  "cluster_psu", "country_phase",
  "total_children_ever_born", "number_living_children",
  "child_sex", "child_age_months",
  "birth_order", "birth_interval_preceding",
  "anc_visits_number", "anc_first_visit_timing",
  "place_of_delivery", "anc_provider",
  "child_vacc_measles_months", "child_vacc_measles",
  "child_vacc_measles_on_time", "child_vacc_measles_late",
  "child_health_card", "vacc_campaign", "child_vitamin_a",
  "household_size", "hh_head_sex", "hh_head_age",
  "relationship_to_hh_head",
  "hh_electricity", "hh_radio", "hh_television",
  "hh_mobile_phone", "hh_bicycle", "hh_motorcycle",
  "hh_car", "hh_computer", "internet_use",
  "media_newspaper", "media_radio", "media_television",
  "health_insurance", "respondent_working", "respondent_working_12mo",
  "respondent_occupation", "partner_occupation",
  "womens_autonomy_no_say", "womens_autonomy_has_a_say",
  "healthcare_barriers", "healthcare_barriers_cost",
  "healthcare_barriers_safety", "healthcare_barriers_transport",
  "healthcare_barriers_stockout",
  "religion", "ethnicity",
  "postnatal_care_newborn", "postnatal_care_newborn_timing",
  "postnatal_care_newborn_person", "postnatal_care_newborn_first_person"
)

# ── Check availability across all surveys ────────────────────────
# Initialize summary dataframe
KR_var_summary <- lapply(names(final_DHS_data), function(survey_name) {
  
  df <- final_DHS_data[[survey_name]]
  
  # Compute percentage of non-missing for each variable
  pct_avail <- sapply(vars_to_check, function(var) {
    if (var %in% names(df)) {
      round(sum(!is.na(df[[var]])) / nrow(df) * 100, 1)
    } else {
      NA  # Variable not present in this survey
    }
  })
  
  # Helper to safely extract unique scalar values
  safe_unique <- function(x) {
    vals <- unique(x[!is.na(x)])
    if (length(vals) == 1) as.character(vals) else NA_character_
  }
  
  # Return as dataframe row
  data.frame(
    FileName     = survey_name,
    CountryCode  = safe_unique(df$DHS_CountryCode),
    CountryName  = safe_unique(df$CountryName),
    SurveyYear   = safe_unique(df$SurveyYear),
    CountryPhase = safe_unique(df$country_phase),
    NRows        = nrow(df),
    as.list(pct_avail),
    stringsAsFactors = FALSE
  )
  
}) %>% dplyr::bind_rows()

# Preview
dim(KR_var_summary)


writexl::write_xlsx(KR_var_summary, "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/data_availability.xlsx")


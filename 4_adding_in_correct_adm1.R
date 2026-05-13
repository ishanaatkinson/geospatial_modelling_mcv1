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

# ================================================================
# Build per-survey region → GADM mapping key + cluster-level ADM2
# ================================================================
#
# PURPOSE:
# Two outputs from a single pass over the surveys:
#
#   (1) REGION MAPPING KEY (region_gadm_mapping_key.xlsx)
#       One row per unique (region_code × ADM1NAME × GADM ADM1)
#       combination per survey, used to align DHS region naming
#       conventions with the canonical GADM ADM1 names.
#
#   (2) CLUSTER-LEVEL ADM1/ADM2 LOOKUP (cluster_adm1_adm2_lookup.xlsx)
#       One row per cluster with: survey, country, year, coords,
#       DHS region info, GADM ADM1 name+code, GADM ADM2 name+code.
#       This is the file you'd join into model_data to add a clean
#       district variable for sub-national analysis.
#
# Strict point-in-polygon assignment is used (no nearest-feature
# snapping) so that misplaced jittered/coastal points don't pollute
# the canonical mapping. Unmatched clusters are counted and reported.
# ================================================================




library(sf)
library(dplyr)
library(purrr)
library(geodata)

# ══════════════════════════════════════════════════════════════════
# STEP 0: Extract survey_name, CountryName, SurveyYear, adm1, 
#         region, n_records from ALL surveys BEFORE any GADM work
# ══════════════════════════════════════════════════════════════════

adm1_region_names_all_surveys <- purrr::map_dfr(names(final_DHS_data), function(sn) {
  
  df <- final_DHS_data[[sn]]
  
  cname <- unique(df$CountryName)[1]
  syear <- as.integer(unique(df$SurveyYear)[1])
  
  has_region   <- "region"   %in% names(df)
  has_adm1name <- "ADM1NAME" %in% names(df)
  
  if (!has_region && !has_adm1name) return(NULL)
  
  df <- df %>%
    mutate(
      ADM1NAME = if (has_adm1name) ADM1NAME else NA_character_,
      region   = if (has_region)   region   else NA
    )
  
  df %>%
    group_by(ADM1NAME, region, .drop = FALSE) %>%
    summarise(n_records = n(), .groups = "drop") %>%
    mutate(
      survey_name = sn,
      CountryName = cname,
      SurveyYear  = as.integer(syear),
      .before     = 1
    ) %>%
    rename(adm1 = ADM1NAME) %>%
    select(survey_name, CountryName, SurveyYear, adm1, region, n_records)
})

adm1_region_names_all_surveys <- adm1_region_names_all_surveys %>%
  arrange(CountryName, survey_name, region)

# Save it out so you have a complete record
writexl::write_xlsx(
  adm1_region_names_all_surveys,
  "YEAR 1/PhD/MEASLES PROJECT/TABLES/adm1_region_names_all_surveys_v2.xlsx"
)

missing <- adm1_region_names_all_surveys %>% filter(survey_name %in% survey_names)

writexl::write_xlsx(
  missing,
  "YEAR 1/PhD/MEASLES PROJECT/TABLES/missing.xlsx"
)

message("Saved ", nrow(adm1_region_names_all_surveys), " rows across ",
        length(unique(adm1_region_names_all_surveys$survey_name)), " surveys")


# ── Country name → GADM ISO3 lookup ─────────────────────────────
country_iso_lookup <- c(
  "Angola"                       = "AGO",
  "Benin"                        = "BEN",
  "Burkina Faso"                 = "BFA",
  "Burundi"                      = "BDI",
  "Cameroon"                     = "CMR",
  "Central African Republic"     = "CAF",
  "Chad"                         = "TCD",
  "Comoros"                      = "COM",
  "Congo"                        = "COG",
  "Congo Democratic Republic"    = "COD",
  "Cote d'Ivoire"                = "CIV",
  "Eritrea"                      = "ERI",
  "Eswatini"                     = "SWZ",
  "Ethiopia"                     = "ETH",
  "Gabon"                        = "GAB",
  "Gambia"                       = "GMB",
  "Ghana"                        = "GHA",
  "Guinea"                       = "GIN",
  "Kenya"                        = "KEN",
  "Lesotho"                      = "LSO",
  "Liberia"                      = "LBR",
  "Madagascar"                   = "MDG",
  "Malawi"                       = "MWI",
  "Mali"                         = "MLI",
  "Mauritania"                   = "MRT",
  "Mozambique"                   = "MOZ",
  "Namibia"                      = "NAM",
  "Niger"                        = "NER",
  "Nigeria"                      = "NGA",
  "Rwanda"                       = "RWA",
  "Sao Tome and Principe"        = "STP",
  "Senegal"                      = "SEN",
  "Sierra Leone"                 = "SLE",
  "South Africa"                 = "ZAF",
  "Sudan"                        = "SDN",
  "Tanzania"                     = "TZA",
  "Togo"                         = "TGO",
  "Uganda"                       = "UGA",
  "Zambia"                       = "ZMB",
  "Zimbabwe"                     = "ZWE"
)


# ── GADM cache ──────────────────────────────────────────────────
gadm_cache <- list()

get_gadm_adm1 <- function(iso3) {
  if (iso3 %in% names(gadm_cache)) {
    cached <- gadm_cache[[iso3]]
    if (inherits(cached, "sf")) return(cached) else return(NULL)
  }
  message("  Downloading GADM ", iso3, " level 1...")
  v <- tryCatch(
    geodata::gadm(country = iso3, level = 1, path = gadm_cache_path),
    error = function(e) { warning(e$message); NULL }
  )
  if (is.null(v)) {
    gadm_cache[[iso3]] <<- NA
    return(NULL)
  }
  v_sf <- sf::st_make_valid(sf::st_as_sf(v))
  gadm_cache[[iso3]] <<- v_sf
  v_sf
}


gadm_table <- purrr::map_dfr(names(final_DHS_data), function(sn) {
  
  df <- final_DHS_data[[sn]]
  
  cname <- unique(df$CountryName)[1]
  syear <- unique(df$SurveyYear)[1]
  
  if (!("cluster_psu" %in% names(df))) return(NULL)
  
  # Skip if no geo
  if (!all(c("LATNUM", "LONGNUM") %in% names(df))) return(NULL)
  if (is.na(cname) || !(cname %in% names(country_iso_lookup))) return(NULL)
  
  iso3  <- country_iso_lookup[[cname]]
  gadm1 <- get_gadm_adm1(iso3)
  
  if (is.null(gadm1)) return(NULL)
  
  pts <- df %>%
    dplyr::filter(!is.na(LATNUM), !is.na(LONGNUM),
                  LATNUM != 0, LONGNUM != 0) %>%
    dplyr::select(cluster_psu, LATNUM, LONGNUM,
                  ADM1NAME, region)
  
  if (nrow(pts) == 0) return(NULL)
  
  pts_sf <- sf::st_as_sf(
    pts,
    coords = c("LONGNUM", "LATNUM"),
    crs = 4326
  )
  
  joined <- sf::st_join(
    pts_sf,
    gadm1[, "NAME_1"],
    join = sf::st_intersects,
    left = TRUE
  ) %>%
    sf::st_drop_geometry() %>%
    dplyr::rename(gadm_adm1_clean = NAME_1)
  
  # 🔑 Aggregate at SAME LEVEL as your final table
  joined %>%
    dplyr::group_by(ADM1NAME, region) %>%
    dplyr::summarise(
      gadm_adm1_clean = {
        tab <- table(gadm_adm1_clean, useNA = "no")
        if (length(tab) == 0) NA_character_ else names(tab)[which.max(tab)]
      },
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      survey_name = sn,
      CountryName = cname,
      SurveyYear  = syear
    ) %>%
    dplyr::rename(
      adm1   = ADM1NAME,
      region = region
    )
})


#
#  IPUMS-DHS identifies these units as "grouped provinces," 
# so there may be inconsistent boundaries compared to later 
# DHS surveys, meaning they cannot be combined directly with 
# later data using the same regional codes. We use here the closest adm1
# unit name according to GADM to align clusters with the given region to 
# an admin1 unit name 
#
#
#
#
#
#



cleanup_table <- purrr::map_dfr(names(final_DHS_data), function(sn) {
  
  df <- final_DHS_data[[sn]]
  
  cname <- unique(df$CountryName)[1]
  syear <- as.integer(unique(df$SurveyYear)[1])
  
  has_region   <- "region" %in% names(df)
  has_adm1name <- "ADM1NAME"    %in% names(df)
  
  # Keep survey if it has EITHER variable
  if (!has_region && !has_adm1name) return(NULL)
  
  # ── Ensure both columns exist ────────────────────────────────
  df <- df %>%
    mutate(
      ADM1NAME    = if (has_adm1name) ADM1NAME else NA_character_,
      region = if (has_region)   region else NA
    )
  
  # ── Build cluster → GADM ADM1 lookup ────────────────────────
  gadm_lookup <- NULL
  
  if (all(c("LATNUM", "LONGNUM", "cluster_psu") %in% names(df)) &&
      !is.na(cname) && cname %in% names(country_iso_lookup)) {
    
    iso3  <- country_iso_lookup[[cname]]
    gadm1 <- get_gadm_adm1(iso3)
    
    if (!is.null(gadm1)) {
      cluster_pts <- df %>%
        filter(!is.na(LATNUM), !is.na(LONGNUM),
               LATNUM != 0, LONGNUM != 0) %>%
        distinct(cluster_psu, LATNUM, LONGNUM)
      
      if (nrow(cluster_pts) > 0) {
        pts_sf <- sf::st_as_sf(
          cluster_pts,
          coords = c("LONGNUM", "LATNUM"),
          crs    = 4326
        )
        
        gadm_lookup <- sf::st_join(
          pts_sf,
          gadm1[, "NAME_1"],
          join = sf::st_intersects,
          left = TRUE
        ) %>%
          sf::st_drop_geometry() %>%
          distinct(cluster_psu, .keep_all = TRUE) %>%
          dplyr::rename(gadm_adm1_clean = NAME_1)
      }
    }
  }
  
  # ── Attach GADM to individual records ───────────────────────
  grp <- df
  
  if (!is.null(gadm_lookup) && "cluster_psu" %in% names(grp)) {
    grp <- grp %>% left_join(gadm_lookup, by = "cluster_psu")
  } else {
    grp$gadm_adm1_clean <- NA_character_
  }
  
  # ── Grouping: keep ALL combinations (including NA cases) ────
  grp %>%
    group_by(ADM1NAME, region, .drop = FALSE) %>%
    summarise(
      n_records = n(),
      gadm_adm1_clean = {
        tab <- table(gadm_adm1_clean, useNA = "no")
        if (length(tab) == 0) NA_character_ else names(tab)[which.max(tab)]
      },
      .groups = "drop"
    ) %>%
    mutate(
      survey_name = sn,
      CountryName = cname,
      SurveyYear  = as.integer(syear),
      .before     = 1
    ) %>%
    rename(
      adm1   = ADM1NAME,
      region = region
    ) %>%
    select(
      survey_name, CountryName, SurveyYear,
      adm1, region, n_records, gadm_adm1_clean
    )
})

# ── Final sorting ─────────────────────────────────────────────
cleanup_table <- cleanup_table %>%
  arrange(
    CountryName,
    tolower(coalesce(as.character(adm1), "")),
    as.character(region),
    as.integer(SurveyYear)
  )

View(cleanup_table)

#writexl::write_xlsx((cleanup_table %>% filter(survey_name %in% survey_names)), paste0(output_folder_table, "/adm1_lookup_code.xlsx"))

output_folder_table <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/TABLES/mcv1_coverage_dhs/"
output_folder_figures <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/FIGURES/mcv1_coverage_dhs/"

adm1_lookup <- readxl::read_xlsx(paste0(output_folder_table, "ADM1_LOOKUP_FINAL.xlsx"))

for (svy in names(final_DHS_data)) {
  
  df <- final_DHS_data[[svy]]
  
  # Ensure columns exist regardless
  if (!"ADM1NAME_CLEAN" %in% names(df)) df$ADM1NAME_CLEAN <- NA_character_
  if (!"REGION_CLEAN"   %in% names(df)) df$REGION_CLEAN   <- NA_integer_
  
  # Get the subset of the lookup for this survey
  lkp <- adm1_lookup[adm1_lookup$survey_name == svy, ]
  
  if (nrow(lkp) == 0) {
    warning(paste("No lookup entries found for:", svy))
    final_DHS_data[[svy]] <- df
    next
  }
  
  # Build mapping, keeping only rows with valid clean values
  region_map <- lkp[!duplicated(lkp$region), c("region", "clean_country", "clean_adm1", "clean_region")]
  region_map <- region_map %>% filter(!is.na(clean_adm1), !is.na(clean_region))
  
  if (nrow(region_map) == 0) {
    final_DHS_data[[svy]] <- df
    next
  }
  
  # Match — returns NA automatically for unmatched rows
  matched_idx <- match(df$region, region_map$region)
  df$ADM1NAME_CLEAN <- region_map$clean_adm1[matched_idx]
  df$REGION_CLEAN   <- region_map$clean_region[matched_idx]
  
  final_DHS_data[[svy]] <- df
}


# ================================
# Save output
# ================================
saveRDS(
  final_DHS_data,
  "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/final_DHS_data.rds"
)
 
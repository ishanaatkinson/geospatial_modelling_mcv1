# ================================================================
# 12_attach_conflict_data.R
# Attach ACLED recent violence indicators at Admin1/2/3
# recent_violence_admN = 1 if >=1 fatal conflict in the cluster's
# AdminN within YEARS_LOOKBACK years prior to (and including)
# SurveyYear, else 0
# ================================================================

library(dplyr)
library(sf)
library(geodata)
library(purrr)
library(readxl)
library(writexl)
library(lubridate)

sf::sf_use_s2(FALSE)

# ================================================================
# CONFIGURATION
# ================================================================

input_rds   <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/all_KR_surveys_with_gps_HR_IR_maps_GC.rds"
output_rds  <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/all_KR_surveys_with_gps_HR_IR_maps_GC_ACLED.rds"
output_log  <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/KR_ACLED_log.xlsx"

acled_xlsx   <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/ACLED Conflict Data/ACLED Data_2026-04-10.xlsx"
gadm_dir    <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/GADM"

dir.create(gadm_dir, showWarnings = FALSE, recursive = TRUE)

YEARS_LOOKBACK <- 2
ADMIN_LEVELS   <- c(1, 2, 3)

# ================================================================
# 1. Load surveys and ACLED
# ================================================================

message("Loading KR surveys...")
all_KR_surveys <- readRDS(input_rds)
cat("Surveys loaded:", length(all_KR_surveys), "\n")

message("Loading ACLED data...")
acled <- read_excel(acled_xlsx) %>%
  mutate(
    event_date = as.Date(event_date),
    year       = as.integer(year),
    fatalities = as.integer(fatalities),
    latitude   = as.numeric(latitude),
    longitude  = as.numeric(longitude)
  ) %>%
  filter(
    !is.na(latitude), !is.na(longitude),
    !is.na(fatalities), fatalities >= 1
  )

cat("ACLED fatal events loaded:", nrow(acled), "\n")
cat("Year range:", min(acled$year), "-", max(acled$year), "\n")

acled_sf <- st_as_sf(acled, coords = c("longitude", "latitude"),
                     crs = 4326, remove = FALSE)

# ================================================================
# 2. Helper: get GADM at a given level (cached, validated)
# ================================================================

get_gadm <- function(iso3, level) {
  tryCatch({
    g <- geodata::gadm(country = iso3, level = level,
                       path = gadm_dir, version = "latest")
    g_sf <- sf::st_as_sf(g)
    g_sf <- sf::st_make_valid(g_sf)
    g_sf
  }, error = function(e) {
    message("    GADM level ", level, " failed for ", iso3, ": ", e$message)
    NULL
  })
}

dhs_to_iso3 <- c(
  AO="AGO", BJ="BEN", BF="BFA", BU="BDI", CM="CMR", CF="CAF", TD="TCD",
  KM="COM", CG="COG", CD="COD", CI="CIV", ER="ERI", SZ="SWZ", ET="ETH",
  GA="GAB", GM="GMB", GH="GHA", GN="GIN", KE="KEN", LS="LSO", LB="LBR",
  MD="MDG", MW="MWI", ML="MLI", MR="MRT", MZ="MOZ", NM="NAM", NI="NER",
  NG="NGA", RW="RWA", ST="STP", SN="SEN", SL="SLE", ZA="ZAF", SD="SDN",
  TZ="TZA", TG="TGO", UG="UGA", ZM="ZMB", ZW="ZWE", AL="ALB", EG="EGY",
  MA="MAR"
)

# ================================================================
# 3. Initialise log and output container
# ================================================================

acled_log <- tibble(
  survey_name           = character(),
  CountryName           = character(),
  SurveyYear            = integer(),
  iso3                  = character(),
  status                = character(),
  n_children            = integer(),
  n_clusters            = integer(),
  n_clusters_geo        = integer(),
  n_children_violence_adm1 = integer(),
  n_children_violence_adm2 = integer(),
  n_children_violence_adm3 = integer(),
  pct_violence_adm1     = numeric(),
  pct_violence_adm2     = numeric(),
  pct_violence_adm3     = numeric(),
  n_acled_events        = integer(),
  error_message         = character()
)

all_KR_surveys_acled <- list()

# Caches keyed by paste(iso3, level)
gadm_cache <- list()

# ================================================================
# 4. Loop over surveys
# ================================================================

survey_names <- names(all_KR_surveys)

for (i in seq_along(survey_names)) {
  
  sname <- survey_names[i]
  df    <- all_KR_surveys[[sname]]
  message(sprintf("[%d/%d] %s", i, length(survey_names), sname))
  
  n_children <- nrow(df)
  
  if (!all(c("v001", "LATNUM", "LONGNUM", "SurveyYear") %in% names(df))) {
    acled_log <- bind_rows(acled_log, tibble(
      survey_name = sname, status = "MISSING_COLS",
      n_children = n_children,
      error_message = "v001/LATNUM/LONGNUM/SurveyYear missing"
    ))
    df$recent_violence_adm1 <- NA_integer_
    df$recent_violence_adm2 <- NA_integer_
    df$recent_violence_adm3 <- NA_integer_
    all_KR_surveys_acled[[sname]] <- df
    next
  }
  
  dhs_cc <- toupper(substr(sname, 1, 2))
  
  if (sname == "OSKR01FL") {
    dhs_cc <- "NG"
  }
  iso3   <- dhs_to_iso3[dhs_cc]
  if (is.na(iso3)) {
    acled_log <- bind_rows(acled_log, tibble(
      survey_name = sname, status = "NO_ISO3",
      n_children = n_children,
      error_message = paste("No ISO3 mapping for DHS code", dhs_cc)
    ))
    df$recent_violence_adm1 <- NA_integer_
    df$recent_violence_adm2 <- NA_integer_
    df$recent_violence_adm3 <- NA_integer_
    all_KR_surveys_acled[[sname]] <- df
    next
  }
  
  survey_year <- as.integer(unique(df$SurveyYear))[1]
  
  # ── Load GADM levels 1, 2, 3 (cached) ──────────────────────
  gadm_levels <- list()
  for (lvl in ADMIN_LEVELS) {
    key <- paste0(iso3, "_", lvl)
    if (is.null(gadm_cache[[key]])) {
      message("  Downloading GADM level ", lvl, " for ", iso3, "...")
      gadm_cache[[key]] <- get_gadm(iso3, lvl)
    }
    gadm_levels[[as.character(lvl)]] <- gadm_cache[[key]]
  }
  
  if (any(sapply(gadm_levels, is.null))) {
    acled_log <- bind_rows(acled_log, tibble(
      survey_name = sname, SurveyYear = survey_year, iso3 = iso3,
      status = "GADM_FAILED", n_children = n_children,
      error_message = "GADM download failed for at least one level"
    ))
    df$recent_violence_adm1 <- NA_integer_
    df$recent_violence_adm2 <- NA_integer_
    df$recent_violence_adm3 <- NA_integer_
    all_KR_surveys_acled[[sname]] <- df
    next
  }
  
  # ── Build cluster sf, drop missing/zero coords ──────────────
  clusters <- df %>%
    distinct(v001, LATNUM, LONGNUM) %>%
    filter(!is.na(LATNUM), !is.na(LONGNUM),
           !(LATNUM == 0 & LONGNUM == 0))
  
  n_clusters     <- n_distinct(df$v001)
  n_clusters_geo <- nrow(clusters)
  
  if (n_clusters_geo == 0) {
    acled_log <- bind_rows(acled_log, tibble(
      survey_name = sname, SurveyYear = survey_year, iso3 = iso3,
      status = "NO_GPS", n_children = n_children,
      n_clusters = n_clusters, n_clusters_geo = 0L
    ))
    df$recent_violence_adm1 <- NA_integer_
    df$recent_violence_adm2 <- NA_integer_
    df$recent_violence_adm3 <- NA_integer_
    all_KR_surveys_acled[[sname]] <- df
    next
  }
  
  clusters_sf <- st_as_sf(clusters,
                          coords = c("LONGNUM", "LATNUM"),
                          crs = 4326)
  
  # ── Filter ACLED to country bbox + time window ─────────────
  yr_min <- survey_year - YEARS_LOOKBACK
  yr_max <- survey_year
  
  country_bbox <- st_as_sfc(st_bbox(gadm_levels[["1"]]))
  acled_window <- acled_sf %>%
    filter(year >= yr_min, year <= yr_max)
  acled_window <- acled_window[st_intersects(acled_window,
                                             country_bbox,
                                             sparse = FALSE)[, 1], ]
  
  # ── Per-level spatial joins ────────────────────────────────
  cluster_flags <- tibble(v001 = clusters$v001)
  log_extra <- list()
  
  for (lvl in ADMIN_LEVELS) {
    
    gid_col   <- paste0("GID_", lvl)
    flag_col  <- paste0("recent_violence_adm", lvl)
    gadm_lvl  <- gadm_levels[[as.character(lvl)]]
    
    if (!gid_col %in% names(gadm_lvl)) {
      cluster_flags[[flag_col]] <- NA_integer_
      log_extra[[paste0("n_children_violence_adm", lvl)]] <- NA_integer_
      log_extra[[paste0("pct_violence_adm",        lvl)]] <- NA_real_
      next
    }
    
    # Cluster -> AdminN
    cl_admin <- st_join(clusters_sf,
                        gadm_lvl[, gid_col],
                        join = st_intersects,
                        left = TRUE) %>%
      st_drop_geometry() %>%
      as_tibble() %>%
      select(v001, !!gid_col)
    
    if (nrow(acled_window) == 0) {
      violent_ids <- character(0)
      n_events_lvl <- 0L
    } else {
      acled_admin <- st_join(acled_window,
                             gadm_lvl[, gid_col],
                             join = st_intersects,
                             left = FALSE) %>%
        st_drop_geometry() %>%
        as_tibble()
      violent_ids <- unique(acled_admin[[gid_col]])
      violent_ids <- violent_ids[!is.na(violent_ids)]
      n_events_lvl <- nrow(acled_admin)
    }
    
    cl_admin[[flag_col]] <- if_else(
      cl_admin[[gid_col]] %in% violent_ids, 1L, 0L,
      missing = NA_integer_
    )
    
    cluster_flags <- cluster_flags %>%
      left_join(cl_admin %>% select(v001, !!flag_col), by = "v001")
  }
  
  # ── Join flags back to child records ───────────────────────
  df <- df %>%
    mutate(v001 = as.numeric(v001)) %>%
    left_join(cluster_flags %>% mutate(v001 = as.numeric(v001)),
              by = "v001")
  
  n_v_adm1 <- sum(df$recent_violence_adm1 == 1, na.rm = TRUE)
  n_v_adm2 <- sum(df$recent_violence_adm2 == 1, na.rm = TRUE)
  n_v_adm3 <- sum(df$recent_violence_adm3 == 1, na.rm = TRUE)
  
  acled_log <- bind_rows(acled_log, tibble(
    survey_name = sname, SurveyYear = survey_year, iso3 = iso3,
    status = "OK", n_children = n_children,
    n_clusters = n_clusters, n_clusters_geo = n_clusters_geo,
    n_children_violence_adm1 = n_v_adm1,
    n_children_violence_adm2 = n_v_adm2,
    n_children_violence_adm3 = n_v_adm3,
    pct_violence_adm1 = round(100 * n_v_adm1 / n_children, 1),
    pct_violence_adm2 = round(100 * n_v_adm2 / n_children, 1),
    pct_violence_adm3 = round(100 * n_v_adm3 / n_children, 1),
    n_acled_events = nrow(acled_window)
  ))
  
  all_KR_surveys_acled[[sname]] <- df
  message(sprintf("  ✅ violence %% adm1/2/3: %.1f / %.1f / %.1f | %d events",
                  100 * n_v_adm1 / n_children,
                  100 * n_v_adm2 / n_children,
                  100 * n_v_adm3 / n_children,
                  nrow(acled_window)))
}

# ================================================================
# 5. Save
# ================================================================

cat("\n=== ACLED ATTACHMENT SUMMARY ===\n")
print(table(acled_log$status, useNA = "ifany"))
cat("\nMean % children with recent violence (OK surveys):\n")
cat("  Admin1:", round(mean(acled_log$pct_violence_adm1[acled_log$status == "OK"], na.rm = TRUE), 1), "%\n")
cat("  Admin2:", round(mean(acled_log$pct_violence_adm2[acled_log$status == "OK"], na.rm = TRUE), 1), "%\n")
cat("  Admin3:", round(mean(acled_log$pct_violence_adm3[acled_log$status == "OK"], na.rm = TRUE), 1), "%\n")

saveRDS(all_KR_surveys_acled, output_rds)
writexl::write_xlsx(list(acled_log = acled_log), output_log)
message("\n✅ Saved: ", output_rds)
message("✅ Log:   ", output_log)

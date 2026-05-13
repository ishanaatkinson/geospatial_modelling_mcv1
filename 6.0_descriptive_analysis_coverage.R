# ================================================================
# 10_descriptive_inequality_analysis.R
# ----------------------------------------------------------------
# Systematic descriptive analysis of MCV1 coverage inequalities:
#   - between countries
#   - within every country, at adm1 / adm2 / cluster level
#   - overall or disaggregated by any population subgroup
#
# Single entry point:  run_inequality_analysis()
#
# Output: per-call Excel file written to `output_folder`, containing:
#   1. coverage          weighted coverage for every (geo unit x subgroup)
#   2. inequality        per-country / per-survey inequality summary
#                        (range, IQR, SD, 90-10 gap, highest/lowest,
#                         ratio, absolute gap, concentration index where
#                         applicable)
#   3. subgroup_summary  per-subgroup coverage with gap vs overall mean
#   4. metadata          parameters of the run
#
# Dependencies: dplyr, purrr, survey, haven, writexl, tidyr, rlang, sf
#   sf + geodata only needed if you want adm2 spatial joins (and you
#   have already produced a cluster -> adm2 lookup upstream).
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(rlang)
  library(survey)
  library(haven)
  library(writexl)
})

options(survey.lonely.psu = "adjust")


# ================================================================
# 1. Subgroup catalogue
# ----------------------------------------------------------------
# Central place to define how each subgroup variable is re-binned
# for disaggregation. Add new subgroups here - the rest of the
# pipeline picks them up automatically.
# ================================================================
subgroup_definitions <- list(
  
  wealth_index_quintile = list(
    var    = "wealth_index_quintile",
    groups = list("Poorest/Poorer" = 0, "Middle" = 1, "Richer/Richest" = 2),
    ordered = TRUE,
    low_group  = "Poorest/Poorer",
    high_group = "Richer/Richest"
  ),
  
  respondent_edu_level = list(
    var    = "respondent_edu_level",
    groups = list("No education" = 0, "Primary" = 1,
                  "Secondary" = 2, "Higher" = 3),
    ordered = TRUE,
    low_group  = "No education",
    high_group = "Higher"
  ),
  
  urban_rural = list(
    var    = "urban_rural",
    groups = list("Rural" = 0, "Urban" = 1),
    ordered = FALSE,
    low_group  = "Rural",
    high_group = "Urban"
  ),
  
  respondent_age_group = list(
    var    = "respondent_age_group",
    groups = list("15-19" = 1, "20-24" = 2, "25-29" = 3,
                  "30-34" = 4, "35-39" = 5, "40-44" = 6, "45-49" = 7),
    ordered = TRUE,
    low_group  = "15-19",
    high_group = "30-34"
  ),
  
  birth_order = list(
    var    = "birth_order",
    groups = list("1st-2nd" = 0, "3rd-5th" = 1, "6th+" = 2),
    ordered = TRUE,
    low_group  = "6th+",
    high_group = "1st-2nd"
  ),
  
  child_sex = list(
    var    = "child_sex",
    groups = list("Male" = 0, "Female" = 1),
    ordered = FALSE,
    low_group  = "Female",
    high_group = "Male"
  ),
  
  anc_visits_number = list(
    var    = "anc_visits_number",
    groups = list("0 visits" = 0, "1-3 visits" = 1, "4+ visits" = 2),
    ordered = TRUE,
    low_group  = "0 visits",
    high_group = "4+ visits"
  ),
  
  place_of_delivery = list(
    var    = "place_of_delivery",
    groups = list("Home" = 1, "Govt facility" = 2, "Private facility" = 3),
    ordered = FALSE,
    low_group  = "Home",
    high_group = "Private facility"
  ),
  
  child_health_card = list(
    var    = "child_health_card",
    groups = list("No card" = 0, "Has card" = 1),
    ordered = FALSE,
    low_group  = "No card",
    high_group = "Has card"
  ),
  
  healthcare_barriers = list(
    var    = "healthcare_barriers",
    groups = list("Had barrier" = 0, "No barrier" = 1),
    ordered = FALSE,
    low_group  = "Had barrier",
    high_group = "No barrier"
  ),
  
  healthcare_barriers_cost = list(
    var    = "healthcare_barriers_cost",
    groups = list("Cost barrier" = 0, "No cost barrier" = 1),
    ordered = FALSE,
    low_group  = "Cost barrier",
    high_group = "No cost barrier"
  ),
  
  healthcare_barriers_transport = list(
    var    = "healthcare_barriers_transport",
    groups = list("Transport barrier" = 0, "No transport barrier" = 1),
    ordered = FALSE,
    low_group  = "Transport barrier",
    high_group = "No transport barrier"
  ),
  
  hh_slum_dwelling = list(
    var    = "hh_slum_dwelling",
    groups = list("Not slum" = 0, "Slum dwelling" = 1),
    ordered = FALSE,
    low_group  = "Slum dwelling",
    high_group = "Not slum"
  ),
  
  media_exposure = list(
    var    = "media_exposure",
    groups = list("No media" = 0, "Any media weekly" = 1),
    ordered = FALSE,
    low_group  = "No media",
    high_group = "Any media weekly"
  ),
  
  mobile_internet_use = list(
    var    = "mobile_internet_use",
    groups = list("No phone/internet" = 0, "Phone or internet" = 1),
    ordered = FALSE,
    low_group  = "No phone/internet",
    high_group = "Phone or internet"
  ),
  
  conflict_area = list(
    var    = "conflict_area",
    groups = list("Conflict" = 0, "No conflict" = 1),
    ordered = FALSE,
    low_group  = "Conflict",
    high_group = "No conflict"
  )
)


# ================================================================
# 2. Helpers
# ================================================================

# Which ID columns define the geographical unit at each level
.geo_vars_for <- function(geo_level) {
  switch(geo_level,
         "country" = c("CountryName", "SurveyYear", "country_phase"),
         "adm1"    = c("CountryName", "SurveyYear", "country_phase",
                       "region", "ADM1NAME"),
         "adm2"    = c("CountryName", "SurveyYear", "country_phase",
                       "region", "ADM1NAME", "adm2_name"),
         "cluster" = c("CountryName", "SurveyYear", "country_phase",
                       "region", "ADM1NAME", "cluster_psu"),
         stop("geo_level must be one of: country, adm1, adm2, cluster")
  )
}

# The "country identity" columns - used when we summarise across
# geographical units within a country (for within-country inequality).
.country_id_vars <- c("CountryName", "SurveyYear", "country_phase")


# ---- weighted quantile (survey weights) -------------------------
.w_quantile <- function(x, w, probs = c(0.1, 0.25, 0.5, 0.75, 0.9)) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (sum(ok) == 0) return(setNames(rep(NA_real_, length(probs)),
                                    paste0("p", probs * 100)))
  x <- x[ok]; w <- w[ok]
  ord <- order(x)
  x <- x[ord]; w <- w[ord]
  cw <- cumsum(w) / sum(w)
  setNames(
    sapply(probs, function(p) x[which(cw >= p)[1]]),
    paste0("p", probs * 100)
  )
}

# ---- weighted mean ----------------------------------------------
.w_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (sum(ok) == 0) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

# ---- concentration index (Wagstaff style) -----------------------
# rank_var: continuous ranking variable (e.g. wealth_index_score);
# outcome:  proportion (coverage); w: survey weights
.concentration_index <- function(outcome, rank_var, w) {
  ok <- is.finite(outcome) & is.finite(rank_var) & is.finite(w) & w > 0
  if (sum(ok) < 30) return(NA_real_)
  o <- outcome[ok]; r <- rank_var[ok]; w <- w[ok]
  ord <- order(r)
  o <- o[ord]; w <- w[ord]; r <- r[ord]
  cw <- cumsum(w) / sum(w)
  frac_rank <- cw - w / (2 * sum(w))         # fractional rank
  mu <- sum(o * w) / sum(w)
  if (!is.finite(mu) || mu == 0) return(NA_real_)
  2 / mu * sum(w * (o - mu) * (frac_rank - 0.5)) / sum(w)
}


# ================================================================
# 3. Per-survey MCV1 coverage with arbitrary geo x disagg grouping
# ----------------------------------------------------------------
# Returns one row per (geo unit x subgroup category) with survey-
# weighted MCV1 coverage, SE, and unweighted n.
# ================================================================
.coverage_one_survey <- function(df, survey_name, geo_vars, subgroup_def) {
  
  if (nrow(df) == 0) return(NULL)
  
  # Zap labels
  df <- df %>% mutate(across(where(haven::is.labelled), haven::zap_labels))
  
  # Core IDs
  if (!"CountryName" %in% names(df)) df$CountryName <- NA_character_
  if (!"SurveyYear"  %in% names(df)) df$SurveyYear  <- NA_character_
  
  df <- df %>%
    mutate(
      CountryName   = as.character(CountryName),
      SurveyYear    = as.character(SurveyYear),
      country_phase = if ("country_phase" %in% names(.)) as.character(country_phase) else NA_character_,
      region        = if ("region"        %in% names(.)) as.character(region)        else NA_character_,
      ADM1NAME      = if ("ADM1NAME"      %in% names(.)) as.character(ADM1NAME)      else NA_character_
    )
  
  # Require child_age_months
  if (!"child_age_months" %in% names(df)) {
    message("  ", survey_name, ": child_age_months missing, skipping")
    return(NULL)
  }
  df <- df %>%
    mutate(child_age_months = as.numeric(child_age_months)) %>%
    filter(!is.na(child_age_months),
           child_age_months >= 12,
           child_age_months <= 23)
  if (nrow(df) == 0) return(NULL)
  
  # MCV1 binary (any source)
  if ("child_vacc_measles" %in% names(df)) {
    df <- df %>% mutate(vacc = as.numeric(child_vacc_measles))
  } else if ("h9" %in% names(df)) {
    df <- df %>% mutate(vacc = as.numeric(h9),
                        vacc = ifelse(vacc %in% c(8, 9), NA, vacc))
  } else {
    message("  ", survey_name, ": no MCV1 variable, skipping")
    return(NULL)
  }
  df <- df %>%
    mutate(vacc      = ifelse(vacc %in% c(0, 1, 2, 3), vacc, NA_real_),
           mcv1_any  = as.integer(vacc %in% c(1, 2, 3)))
  df$mcv1_any[is.na(df$vacc)] <- NA_integer_
  
  # Disaggregation category
  if (is.null(subgroup_def)) {
    df$disagg <- "All"
  } else {
    v <- subgroup_def$var
    if (!v %in% names(df)) {
      df$disagg <- NA_character_
    } else {
      raw <- suppressWarnings(as.numeric(df[[v]]))
      labels <- names(subgroup_def$groups)
      values <- subgroup_def$groups
      cat_vec <- rep(NA_character_, length(raw))
      for (lbl in labels) {
        cat_vec[raw %in% values[[lbl]]] <- lbl
      }
      df$disagg <- cat_vec
    }
  }
  df <- df %>% filter(!is.na(disagg))
  if (nrow(df) == 0) return(NULL)
  
  # Ensure all geo_vars exist (for adm2 we expect upstream join)
  for (gv in geo_vars) {
    if (!gv %in% names(df)) df[[gv]] <- NA_character_
  }
  # ADM1NAME fallback to region code for surveys without GPS
  if ("ADM1NAME" %in% geo_vars && all(is.na(df$ADM1NAME))) {
    df$ADM1NAME <- df$region
  }
  df <- df %>% mutate(across(all_of(geo_vars), as.character))
  
  # Survey design variables
  if (!"v005" %in% names(df)) {
    message("  ", survey_name, ": v005 missing, skipping")
    return(NULL)
  }
  df$weight <- as.numeric(df$v005) / 1e6
  
  use_nest <- TRUE
  if ("v023" %in% names(df) && !all(is.na(df$v023))) {
    df$strata_var <- as.numeric(df$v023)
  } else if (all(c("v024", "v025") %in% names(df))) {
    df$strata_var <- as.integer(interaction(df$v024, df$v025, drop = TRUE))
  } else {
    df$strata_var <- rep(1L, nrow(df)); use_nest <- FALSE
  }
  df$psu_var <- if ("v021" %in% names(df) && !all(is.na(df$v021))) {
    as.numeric(df$v021)
  } else if ("cluster_psu" %in% names(df) && !all(is.na(df$cluster_psu))) {
    as.numeric(df$cluster_psu)
  } else if ("v001" %in% names(df)) {
    as.numeric(df$v001)
  } else seq_len(nrow(df))
  
  keep_cols <- unique(c(geo_vars, "disagg", "mcv1_any",
                        "weight", "strata_var", "psu_var"))
  df_design <- df %>% select(all_of(keep_cols)) %>% filter(!is.na(mcv1_any))
  if (nrow(df_design) == 0) return(NULL)
  
  # Survey design
  design <- tryCatch(
    survey::svydesign(
      ids     = ~psu_var,
      strata  = ~strata_var,
      weights = ~weight,
      data    = df_design,
      nest    = use_nest
    ),
    error = function(e) { message("  ", survey_name, ": ", e$message); NULL }
  )
  if (is.null(design)) return(NULL)
  
  # Weighted coverage
  group_formula <- as.formula(paste("~", paste(c(geo_vars, "disagg"),
                                               collapse = " + ")))
  coverage <- tryCatch(
    survey::svyby(~mcv1_any, group_formula, design = design,
                  FUN = survey::svymean, na.rm = TRUE, keep.var = TRUE),
    error = function(e) { message("  ", survey_name, ": svyby ", e$message); NULL }
  )
  if (is.null(coverage)) return(NULL)
  
  # Unweighted n (total and vaccinated)
  denom <- df_design %>%
    group_by(across(all_of(c(geo_vars, "disagg")))) %>%
    summarise(
      n_total_12_23 = n(),
      n_vacc        = sum(mcv1_any == 1, na.rm = TRUE),
      .groups = "drop"
    )
  
  coverage %>%
    as_tibble() %>%
    rename(coverage_prop = mcv1_any, se_prop = se) %>%
    mutate(coverage_pct = round(coverage_prop * 100, 2),
           se_pct       = round(se_prop       * 100, 2)) %>%
    left_join(denom, by = c(geo_vars, "disagg")) %>%
    mutate(survey_name = survey_name, .before = 1)
}


# ================================================================
# 4. Inequality summary per grouping unit
# ----------------------------------------------------------------
# For within-country geographical dispersion we summarise across
# ADM1 / ADM2 / cluster units *within each country-survey*.
# For between-country dispersion we summarise across countries
# *within each survey era* (or, simpler, treat each country-survey
# as its own row and dispersion is just the distribution of those).
# ================================================================

# (a) Geographical dispersion summary
.geo_dispersion <- function(coverage_df, geo_level, scope) {
  
  # coverage_df here is the overall coverage per geo unit (no subgroup),
  # i.e. disagg == "All".
  coverage_df <- coverage_df %>% filter(disagg == "All")
  
  if (scope == "between_countries") {
    # One row per country-survey; dispersion is across country means.
    out <- coverage_df %>%
      group_by(SurveyYear) %>%
      summarise(
        n_units       = n(),
        mean_cov      = mean(coverage_pct, na.rm = TRUE),
        median_cov    = median(coverage_pct, na.rm = TRUE),
        sd_cov        = sd(coverage_pct, na.rm = TRUE),
        min_cov       = min(coverage_pct, na.rm = TRUE),
        max_cov       = max(coverage_pct, na.rm = TRUE),
        p10_cov       = quantile(coverage_pct, 0.1, na.rm = TRUE),
        p90_cov       = quantile(coverage_pct, 0.9, na.rm = TRUE),
        range_cov     = max_cov - min_cov,
        p90_p10_gap   = p90_cov - p10_cov,
        iqr_cov       = IQR(coverage_pct, na.rm = TRUE),
        prop_below_50 = mean(coverage_pct < 50,  na.rm = TRUE),
        prop_below_80 = mean(coverage_pct < 80,  na.rm = TRUE),
        prop_above_95 = mean(coverage_pct >= 95, na.rm = TRUE),
        .groups = "drop"
      )
    return(out)
  }
  
  # within_country: dispersion across adm1/adm2/cluster within each
  # (country, survey)
  out <- coverage_df %>%
    group_by(across(all_of(.country_id_vars))) %>%
    summarise(
      geo_level     = geo_level,
      n_units       = n(),
      mean_cov      = mean(coverage_pct, na.rm = TRUE),
      median_cov    = median(coverage_pct, na.rm = TRUE),
      sd_cov        = sd(coverage_pct, na.rm = TRUE),
      min_cov       = min(coverage_pct, na.rm = TRUE),
      max_cov       = max(coverage_pct, na.rm = TRUE),
      p10_cov       = quantile(coverage_pct, 0.1, na.rm = TRUE),
      p90_cov       = quantile(coverage_pct, 0.9, na.rm = TRUE),
      range_cov     = max_cov - min_cov,
      p90_p10_gap   = p90_cov - p10_cov,
      iqr_cov       = IQR(coverage_pct, na.rm = TRUE),
      prop_below_50 = mean(coverage_pct < 50,  na.rm = TRUE),
      prop_below_80 = mean(coverage_pct < 80,  na.rm = TRUE),
      prop_above_95 = mean(coverage_pct >= 95, na.rm = TRUE),
      .groups = "drop"
    )
  out
}

# (b) Subgroup gap summary per grouping unit
#
# For each (geo unit), compute:
#   - coverage in low_group, high_group, overall
#   - absolute gap (high - low)
#   - relative ratio (high / low)
#   - max - min across any pair of subgroup categories
.subgroup_gap_summary <- function(coverage_df, geo_vars, subgroup_def) {
  
  if (is.null(subgroup_def)) return(NULL)
  if (!"disagg" %in% names(coverage_df)) return(NULL)
  if (all(coverage_df$disagg == "All")) return(NULL)
  
  low  <- subgroup_def$low_group
  high <- subgroup_def$high_group
  
  wide <- coverage_df %>%
    select(all_of(geo_vars), disagg, coverage_pct) %>%
    pivot_wider(names_from  = disagg,
                values_from = coverage_pct,
                names_prefix = "cov_")
  
  # Range across all present subgroup categories (robust to
  # country-surveys that only have a subset of categories)
  cat_cols <- setdiff(names(wide), geo_vars)
  wide$max_subgroup_cov  <- apply(wide[cat_cols], 1,
                                  function(x) suppressWarnings(max(x, na.rm = TRUE)))
  wide$min_subgroup_cov  <- apply(wide[cat_cols], 1,
                                  function(x) suppressWarnings(min(x, na.rm = TRUE)))
  wide$subgroup_range    <- wide$max_subgroup_cov - wide$min_subgroup_cov
  
  # High/low anchored gap (where both anchors are present)
  low_col  <- paste0("cov_", low)
  high_col <- paste0("cov_", high)
  if (low_col %in% names(wide) && high_col %in% names(wide)) {
    wide$absolute_gap_high_minus_low <- wide[[high_col]] - wide[[low_col]]
    wide$relative_ratio_high_over_low <- ifelse(
      !is.na(wide[[low_col]]) & wide[[low_col]] > 0,
      wide[[high_col]] / wide[[low_col]], NA_real_)
  } else {
    wide$absolute_gap_high_minus_low <- NA_real_
    wide$relative_ratio_high_over_low <- NA_real_
  }
  
  wide %>%
    mutate(
      low_anchor  = low,
      high_anchor = high,
      max_subgroup_cov = ifelse(is.infinite(max_subgroup_cov), NA_real_,
                                max_subgroup_cov),
      min_subgroup_cov = ifelse(is.infinite(min_subgroup_cov), NA_real_,
                                min_subgroup_cov)
    )
}


# ================================================================
# 5. Main entry point
# ----------------------------------------------------------------
# run_inequality_analysis(
#   surveys        = <named list of DHS data.frames>,
#   scope          = "between_countries" or "within_country",
#   geo_level      = "country" | "adm1" | "adm2" | "cluster",
#   subgroup       = NULL or a name in subgroup_definitions,
#   countries      = NULL (all) or character vector of CountryName,
#   output_folder  = directory to write .xlsx into,
#   file_prefix    = prefix for the filename
# )
# ----------------------------------------------------------------
# Rules:
#   scope = "between_countries" forces geo_level = "country"
#   scope = "within_country"    requires geo_level in {adm1, adm2, cluster}
#   subgroup = NULL gives overall coverage + geo dispersion only
#   subgroup = <name> adds subgroup disaggregation and gap tables
# ================================================================
run_inequality_analysis <- function(
    surveys,
    scope          = c("between_countries", "within_country"),
    geo_level      = c("country", "adm1", "adm2", "cluster"),
    subgroup       = NULL,
    countries      = NULL,
    output_folder  = ".",
    file_prefix    = "mcv1_inequality"
) {
  
  scope     <- match.arg(scope)
  geo_level <- match.arg(geo_level)
  
  # Enforce scope / geo_level coherence
  if (scope == "between_countries" && geo_level != "country") {
    message("scope=between_countries forces geo_level=country; overriding.")
    geo_level <- "country"
  }
  if (scope == "within_country" && geo_level == "country") {
    stop("scope=within_country requires geo_level in {adm1, adm2, cluster}.")
  }
  
  # Subgroup spec
  sg_def <- NULL
  if (!is.null(subgroup)) {
    if (!subgroup %in% names(subgroup_definitions)) {
      stop("Unknown subgroup '", subgroup,
           "'. Options: ", paste(names(subgroup_definitions), collapse = ", "))
    }
    sg_def <- subgroup_definitions[[subgroup]]
  }
  
  # Country filter
  if (!is.null(countries)) {
    surveys <- keep(surveys, function(s) {
      is.data.frame(s) &&
        "CountryName" %in% names(s) &&
        any(as.character(s$CountryName) %in% countries)
    })
  }
  
  geo_vars <- .geo_vars_for(geo_level)
  
  message("Running: scope=", scope, " | geo_level=", geo_level,
          " | subgroup=", ifelse(is.null(subgroup), "none", subgroup),
          " | n_surveys=", length(surveys))
  
  # --- 1. Overall coverage at the requested geo_level (no subgroup)
  coverage_overall <- purrr::imap_dfr(surveys, function(df, nm) {
    if (!is.data.frame(df)) return(NULL)
    .coverage_one_survey(df, nm, geo_vars, subgroup_def = NULL)
  })
  
  # --- 2. Subgroup-disaggregated coverage (if requested)
  coverage_by_subgroup <- NULL
  if (!is.null(sg_def)) {
    coverage_by_subgroup <- purrr::imap_dfr(surveys, function(df, nm) {
      if (!is.data.frame(df)) return(NULL)
      .coverage_one_survey(df, nm, geo_vars, subgroup_def = sg_def)
    })
  }
  
  # --- 3. Geographical dispersion summary
  geo_dispersion <- .geo_dispersion(coverage_overall, geo_level, scope)
  
  # --- 4. Subgroup gap summary (if subgroup requested)
  subgroup_gaps <- NULL
  subgroup_means <- NULL
  if (!is.null(coverage_by_subgroup)) {
    
    subgroup_gaps <- .subgroup_gap_summary(coverage_by_subgroup,
                                           geo_vars, sg_def)
    
    # Per-subgroup summary across units (eg pooled across ADM1 within country)
    group_cols <- if (scope == "between_countries") {
      c("SurveyYear", "disagg")
    } else {
      c(.country_id_vars, "disagg")
    }
    overall_per_unit <- coverage_overall %>%
      select(all_of(c(geo_vars, "coverage_pct"))) %>%
      rename(overall_cov = coverage_pct)
    
    subgroup_means <- coverage_by_subgroup %>%
      left_join(overall_per_unit, by = geo_vars) %>%
      mutate(gap_vs_overall = coverage_pct - overall_cov) %>%
      group_by(across(all_of(group_cols))) %>%
      summarise(
        n_units            = n(),
        mean_subgroup_cov  = mean(coverage_pct,    na.rm = TRUE),
        median_subgroup    = median(coverage_pct,  na.rm = TRUE),
        mean_overall_cov   = mean(overall_cov,     na.rm = TRUE),
        mean_gap_vs_overall = mean(gap_vs_overall, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  # --- 5. Metadata
  metadata <- tibble::tibble(
    parameter = c("scope", "geo_level", "subgroup",
                  "n_surveys_input", "n_surveys_used",
                  "countries_filter", "timestamp",
                  "low_anchor", "high_anchor"),
    value = c(
      scope, geo_level,
      ifelse(is.null(subgroup), "none", subgroup),
      as.character(length(surveys)),
      as.character(length(unique(coverage_overall$survey_name))),
      ifelse(is.null(countries), "all",
             paste(countries, collapse = "; ")),
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      ifelse(is.null(sg_def), NA, sg_def$low_group),
      ifelse(is.null(sg_def), NA, sg_def$high_group)
    )
  )
  
  # --- 6. Assemble sheets and save
  sheets <- list(
    metadata        = metadata,
    coverage        = coverage_overall,
    geo_dispersion  = geo_dispersion
  )
  if (!is.null(coverage_by_subgroup)) {
    sheets$coverage_by_subgroup <- coverage_by_subgroup
    sheets$subgroup_gaps        <- subgroup_gaps
    sheets$subgroup_summary     <- subgroup_means
  }
  
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  sub_tag <- ifelse(is.null(subgroup), "none", subgroup)
  out_file <- file.path(
    output_folder,
    sprintf("%s_%s_%s_by_%s.xlsx",
            file_prefix, scope, geo_level, sub_tag)
  )
  writexl::write_xlsx(sheets, out_file)
  message("Saved: ", out_file)
  
  invisible(sheets)
}


# ================================================================
# 6. Batch runner - run ALL combinations in one go
# ----------------------------------------------------------------
# Use this for the full descriptive chapter: it loops over every
# scope x geo_level x subgroup combination and writes one xlsx per
# combination, plus a single combined summary.
# ================================================================
run_all_inequality_analyses <- function(
    surveys,
    output_folder = ".",
    subgroups     = names(subgroup_definitions),
    include_no_subgroup = TRUE
) {
  
  specs <- list()
  
  # Between-countries (geo_level fixed to country)
  if (include_no_subgroup) {
    specs[[length(specs) + 1]] <- list(scope = "between_countries",
                                       geo_level = "country",
                                       subgroup = NULL)
  }
  for (sg in subgroups) {
    specs[[length(specs) + 1]] <- list(scope = "between_countries",
                                       geo_level = "country",
                                       subgroup = sg)
  }
  
  # Within-country at adm1 / adm2 / cluster
  for (gl in c("adm1", "adm2", "cluster")) {
    if (include_no_subgroup) {
      specs[[length(specs) + 1]] <- list(scope = "within_country",
                                         geo_level = gl,
                                         subgroup = NULL)
    }
    for (sg in subgroups) {
      specs[[length(specs) + 1]] <- list(scope = "within_country",
                                         geo_level = gl,
                                         subgroup = sg)
    }
  }
  
  message("Batch run: ", length(specs), " combinations to process.")
  
  results <- purrr::imap(specs, function(s, i) {
    message(sprintf("[%d/%d] scope=%s geo=%s sub=%s",
                    i, length(specs), s$scope, s$geo_level,
                    ifelse(is.null(s$subgroup), "none", s$subgroup)))
    tryCatch(
      run_inequality_analysis(
        surveys       = surveys,
        scope         = s$scope,
        geo_level     = s$geo_level,
        subgroup      = s$subgroup,
        output_folder = output_folder
      ),
      error = function(e) {
        message("    FAILED: ", e$message); NULL
      }
    )
  })
  
  invisible(results)
}


# ================================================================
# 7. Example usage (commented out - uncomment to run)
# ================================================================
#
# final_DHS_data <- readRDS(
#   "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/final_DHS_data.rds"
# )
#
# output_folder <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/TABLES/mcv1_inequality"
#
# # Example 1: between-country dispersion, overall (no subgroup)
# run_inequality_analysis(
#   surveys       = final_DHS_data,
#   scope         = "between_countries",
#   geo_level     = "country",
#   subgroup      = NULL,
#   output_folder = output_folder
# )
#
# # Example 2: within-country ADM1 inequality, by wealth quintile
# run_inequality_analysis(
#   surveys       = final_DHS_data,
#   scope         = "within_country",
#   geo_level     = "adm1",
#   subgroup      = "wealth_index_quintile",
#   output_folder = output_folder
# )
#
# # Example 3: Nigeria only, cluster-level, by maternal education
# run_inequality_analysis(
#   surveys       = final_DHS_data,
#   scope         = "within_country",
#   geo_level     = "cluster",
#   subgroup      = "respondent_edu_level",
#   countries     = "Nigeria",
#   output_folder = output_folder
# )
#
# # Example 4: run everything in one batch
# run_all_inequality_analyses(
#   surveys       = final_DHS_data,
#   output_folder = output_folder
# )
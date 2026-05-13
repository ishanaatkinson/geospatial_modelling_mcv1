# ================================================================
# 11_temporal_inequality_analysis.R
# ----------------------------------------------------------------
# Systematic temporal analysis of MCV1 coverage inequalities:
# how have coverage gaps changed over time, at every geographical
# level and for every population subgroup.
#
# Builds on 10_descriptive_inequality_analysis.R - reuses the same
# subgroup catalogue and per-survey coverage helpers.
#
# Single entry point:  run_temporal_inequality_analysis()
#
# Output: per-call Excel file with:
#   1. metadata           run parameters
#   2. coverage_ts        per-unit coverage at every survey time point
#   3. country_trends     per-country coverage trajectory (first,
#                         latest, change, annualised slope)
#   4. geo_dispersion_ts  per-country geographical dispersion at
#                         every time point + change in dispersion
#   5. subgroup_trends    per-country per-subgroup trajectory + change
#                         in absolute gap and relative ratio
#   6. convergence_summary  which subgroups/geos widened vs narrowed
#                           across the sample of countries
# ================================================================

# Requires the descriptive script to have been sourced first, so
# that subgroup_definitions, .coverage_one_survey, .geo_vars_for and
# the other helpers are in scope.
#
# source("10_descriptive_inequality_analysis.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(rlang)
  library(writexl)
})


# ================================================================
# 1. Helpers for temporal slopes and change classification
# ================================================================

# Classify a change using a simple "reliable change" heuristic:
# non-overlapping 95% CIs (where available) OR |change| >= threshold pp.
.classify_change <- function(change_pp, se_first = NA, se_latest = NA,
                             threshold_pp = 3) {
  if (is.na(change_pp)) return(NA_character_)
  # If we have SEs for both endpoints, test overlap of 95% CIs
  if (!is.na(se_first) && !is.na(se_latest)) {
    ci_width <- 1.96 * sqrt(se_first^2 + se_latest^2)
    if (abs(change_pp) > ci_width) {
      return(if (change_pp > 0) "increase" else "decline")
    } else {
      return("no reliable change")
    }
  }
  # Otherwise fall back to absolute-threshold rule
  if (abs(change_pp) < threshold_pp) "no reliable change"
  else if (change_pp > 0) "increase"
  else "decline"
}

# Fit a weighted linear slope (pp per year) across survey years,
# weighting by inverse variance if SEs are available.
.fit_annual_slope <- function(years, values, ses = NULL) {
  ok <- is.finite(years) & is.finite(values)
  if (sum(ok) < 2) return(c(slope = NA_real_, slope_se = NA_real_,
                            n_points = sum(ok)))
  years <- years[ok]; values <- values[ok]
  if (!is.null(ses)) ses <- ses[ok]
  w <- if (!is.null(ses) && all(is.finite(ses)) && all(ses > 0)) {
    1 / ses^2
  } else {
    rep(1, length(years))
  }
  fit <- tryCatch(stats::lm(values ~ years, weights = w),
                  error = function(e) NULL)
  if (is.null(fit)) return(c(slope = NA_real_, slope_se = NA_real_,
                             n_points = length(years)))
  co <- summary(fit)$coefficients
  c(slope    = unname(co["years", "Estimate"]),
    slope_se = unname(co["years", "Std. Error"]),
    n_points = length(years))
}


# ================================================================
# 2. Country-level coverage trajectory
# ----------------------------------------------------------------
# Input: overall coverage_ts from descriptive pipeline (one row per
# survey, disagg == "All"). Output: one row per country with first,
# latest, change, annualised slope, and change classification.
# ================================================================
.country_trajectory <- function(coverage_ts) {
  
  ts <- coverage_ts %>%
    filter(disagg == "All") %>%
    mutate(year_num = suppressWarnings(as.numeric(SurveyYear))) %>%
    filter(!is.na(year_num))
  
  if (nrow(ts) == 0) return(NULL)
  
  ts %>%
    group_by(CountryName) %>%
    arrange(year_num, .by_group = TRUE) %>%
    summarise(
      n_surveys        = n(),
      first_year       = first(year_num),
      latest_year      = last(year_num),
      years_span       = latest_year - first_year,
      first_coverage   = first(coverage_pct),
      latest_coverage  = last(coverage_pct),
      first_se         = first(se_pct),
      latest_se        = last(se_pct),
      abs_change_pp    = latest_coverage - first_coverage,
      rel_change_pct   = ifelse(first_coverage > 0,
                                (latest_coverage - first_coverage) /
                                  first_coverage * 100, NA_real_),
      ann_slope_pp_yr  = .fit_annual_slope(year_num, coverage_pct,
                                           se_pct)["slope"],
      ann_slope_se     = .fit_annual_slope(year_num, coverage_pct,
                                           se_pct)["slope_se"],
      .groups = "drop"
    ) %>%
    mutate(
      change_class = mapply(.classify_change,
                            abs_change_pp, first_se, latest_se)
    )
}


# ================================================================
# 3. Geographical dispersion trajectory (within-country, over time)
# ----------------------------------------------------------------
# For each country, at each survey year, compute the sub-national
# dispersion in coverage across the chosen geo_level (adm1, adm2 or
# cluster). Then report the change in dispersion from first to latest.
# ================================================================
.geo_dispersion_trajectory <- function(coverage_ts) {
  
  ts <- coverage_ts %>%
    filter(disagg == "All") %>%
    mutate(year_num = suppressWarnings(as.numeric(SurveyYear))) %>%
    filter(!is.na(year_num))
  
  if (nrow(ts) == 0) return(list(ts = NULL, change = NULL))
  
  # Per country x survey, dispersion across units
  per_survey <- ts %>%
    group_by(CountryName, SurveyYear, year_num, country_phase) %>%
    summarise(
      n_units        = n(),
      mean_cov       = mean(coverage_pct, na.rm = TRUE),
      median_cov     = median(coverage_pct, na.rm = TRUE),
      sd_cov         = sd(coverage_pct, na.rm = TRUE),
      min_cov        = min(coverage_pct, na.rm = TRUE),
      max_cov        = max(coverage_pct, na.rm = TRUE),
      p10_cov        = quantile(coverage_pct, 0.1, na.rm = TRUE),
      p90_cov        = quantile(coverage_pct, 0.9, na.rm = TRUE),
      iqr_cov        = IQR(coverage_pct, na.rm = TRUE),
      range_cov      = max_cov - min_cov,
      p90_p10_gap    = p90_cov - p10_cov,
      prop_below_50  = mean(coverage_pct < 50,  na.rm = TRUE),
      prop_below_80  = mean(coverage_pct < 80,  na.rm = TRUE),
      prop_above_95  = mean(coverage_pct >= 95, na.rm = TRUE),
      .groups = "drop"
    )
  
  # First vs latest change per country
  change <- per_survey %>%
    group_by(CountryName) %>%
    arrange(year_num, .by_group = TRUE) %>%
    summarise(
      n_surveys          = n(),
      first_year         = first(year_num),
      latest_year        = last(year_num),
      first_range        = first(range_cov),
      latest_range       = last(range_cov),
      range_change_pp    = latest_range - first_range,
      first_iqr          = first(iqr_cov),
      latest_iqr         = last(iqr_cov),
      iqr_change_pp      = latest_iqr - first_iqr,
      first_p90_p10      = first(p90_p10_gap),
      latest_p90_p10     = last(p90_p10_gap),
      p90_p10_change_pp  = latest_p90_p10 - first_p90_p10,
      first_sd           = first(sd_cov),
      latest_sd          = last(sd_cov),
      sd_change_pp       = latest_sd - first_sd,
      first_min          = first(min_cov),
      latest_min         = last(min_cov),
      min_change_pp      = latest_min - first_min,
      first_max          = first(max_cov),
      latest_max         = last(max_cov),
      max_change_pp      = latest_max - first_max,
      .groups = "drop"
    ) %>%
    mutate(
      dispersion_direction = case_when(
        is.na(range_change_pp)        ~ NA_character_,
        range_change_pp < -3          ~ "narrowing",
        range_change_pp >  3          ~ "widening",
        TRUE                          ~ "stable"
      ),
      worst_unit_catchup = case_when(
        is.na(min_change_pp) | is.na(max_change_pp) ~ NA_character_,
        min_change_pp > max_change_pp + 3           ~ "worst catching up",
        max_change_pp > min_change_pp + 3           ~ "best pulling ahead",
        TRUE                                        ~ "parallel change"
      )
    )
  
  list(ts = per_survey, change = change)
}


# ================================================================
# 4. Subgroup gap trajectory
# ----------------------------------------------------------------
# For each country x survey year, get coverage per subgroup category.
# Compute absolute gap (high - low) and relative ratio (high / low)
# using the subgroup_def anchors. Then report change in both metrics
# from first to latest survey.
# ================================================================
.subgroup_trajectory <- function(coverage_ts_subgroup, subgroup_def) {
  
  if (is.null(coverage_ts_subgroup)) return(list(ts = NULL, change = NULL))
  if (!"disagg" %in% names(coverage_ts_subgroup)) {
    return(list(ts = NULL, change = NULL))
  }
  if (all(coverage_ts_subgroup$disagg == "All")) {
    return(list(ts = NULL, change = NULL))
  }
  
  low  <- subgroup_def$low_group
  high <- subgroup_def$high_group
  
  ts <- coverage_ts_subgroup %>%
    mutate(year_num = suppressWarnings(as.numeric(SurveyYear))) %>%
    filter(!is.na(year_num))
  
  # Per country x survey x disagg, pool across subnational units
  # (i.e. for geo_level=country this is one row per disagg; for adm1
  # etc. this pools the ADM1 unit values up to a country mean within
  # each subgroup, weighted by unweighted sample size).
  per_survey_sub <- ts %>%
    group_by(CountryName, SurveyYear, year_num, country_phase, disagg) %>%
    summarise(
      n_units       = n(),
      n_total       = sum(n_total_12_23, na.rm = TRUE),
      pooled_cov    = weighted.mean(coverage_pct,
                                    w = n_total_12_23, na.rm = TRUE),
      mean_cov      = mean(coverage_pct, na.rm = TRUE),
      median_cov    = median(coverage_pct, na.rm = TRUE),
      sd_cov        = sd(coverage_pct, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Pivot wide to compute gaps per country x survey
  gaps_wide <- per_survey_sub %>%
    select(CountryName, SurveyYear, year_num, country_phase,
           disagg, pooled_cov) %>%
    pivot_wider(names_from  = disagg,
                values_from = pooled_cov,
                names_prefix = "cov_")
  
  cat_cols <- setdiff(names(gaps_wide),
                      c("CountryName", "SurveyYear", "year_num",
                        "country_phase"))
  
  gaps_wide$max_subgroup_cov <- apply(gaps_wide[cat_cols], 1, function(x)
    suppressWarnings(max(x, na.rm = TRUE)))
  gaps_wide$min_subgroup_cov <- apply(gaps_wide[cat_cols], 1, function(x)
    suppressWarnings(min(x, na.rm = TRUE)))
  gaps_wide$max_subgroup_cov <- ifelse(is.infinite(gaps_wide$max_subgroup_cov),
                                       NA_real_, gaps_wide$max_subgroup_cov)
  gaps_wide$min_subgroup_cov <- ifelse(is.infinite(gaps_wide$min_subgroup_cov),
                                       NA_real_, gaps_wide$min_subgroup_cov)
  gaps_wide$subgroup_range <- gaps_wide$max_subgroup_cov -
    gaps_wide$min_subgroup_cov
  
  low_col  <- paste0("cov_", low)
  high_col <- paste0("cov_", high)
  if (low_col %in% names(gaps_wide) && high_col %in% names(gaps_wide)) {
    gaps_wide$absolute_gap   <- gaps_wide[[high_col]] - gaps_wide[[low_col]]
    gaps_wide$relative_ratio <- ifelse(
      !is.na(gaps_wide[[low_col]]) & gaps_wide[[low_col]] > 0,
      gaps_wide[[high_col]] / gaps_wide[[low_col]], NA_real_)
  } else {
    gaps_wide$absolute_gap   <- NA_real_
    gaps_wide$relative_ratio <- NA_real_
  }
  gaps_wide$low_anchor  <- low
  gaps_wide$high_anchor <- high
  
  # Per-country change from first to latest
  change <- gaps_wide %>%
    group_by(CountryName) %>%
    arrange(year_num, .by_group = TRUE) %>%
    summarise(
      n_surveys            = n(),
      first_year           = first(year_num),
      latest_year          = last(year_num),
      
      first_low_cov        = first(.data[[low_col]]),
      latest_low_cov       = last(.data[[low_col]]),
      low_change_pp        = latest_low_cov - first_low_cov,
      
      first_high_cov       = first(.data[[high_col]]),
      latest_high_cov      = last(.data[[high_col]]),
      high_change_pp       = latest_high_cov - first_high_cov,
      
      first_abs_gap        = first(absolute_gap),
      latest_abs_gap       = last(absolute_gap),
      abs_gap_change_pp    = latest_abs_gap - first_abs_gap,
      
      first_rel_ratio      = first(relative_ratio),
      latest_rel_ratio     = last(relative_ratio),
      rel_ratio_change     = latest_rel_ratio - first_rel_ratio,
      
      first_subgroup_range = first(subgroup_range),
      latest_subgroup_range = last(subgroup_range),
      range_change_pp      = latest_subgroup_range - first_subgroup_range,
      
      ann_slope_abs_gap    = .fit_annual_slope(year_num, absolute_gap)["slope"],
      .groups = "drop"
    ) %>%
    mutate(
      gap_direction = case_when(
        is.na(abs_gap_change_pp)        ~ NA_character_,
        abs_gap_change_pp < -3          ~ "narrowing",
        abs_gap_change_pp >  3          ~ "widening",
        TRUE                            ~ "stable"
      ),
      pro_poor_progress = case_when(
        is.na(low_change_pp) | is.na(high_change_pp) ~ NA_character_,
        low_change_pp > high_change_pp + 3           ~ "low-anchor catching up",
        high_change_pp > low_change_pp + 3           ~ "high-anchor pulling ahead",
        TRUE                                         ~ "parallel change"
      ),
      low_anchor  = low,
      high_anchor = high
    )
  
  list(ts = gaps_wide, change = change)
}


# ================================================================
# 5. Cross-country convergence summary
# ----------------------------------------------------------------
# Pools across all countries: of those with >= 2 surveys, how many
# narrowed / widened / were stable on each metric? Gives the headline
# equity-trajectory numbers for the results chapter.
# ================================================================
.convergence_summary <- function(geo_change = NULL,
                                 subgroup_change = NULL) {
  
  out <- list()
  
  if (!is.null(geo_change) && "dispersion_direction" %in% names(geo_change)) {
    out$geo_dispersion <- geo_change %>%
      count(dispersion_direction, name = "n_countries") %>%
      mutate(pct_countries = round(100 * n_countries /
                                     sum(n_countries), 1),
             metric = "geographical dispersion (range across units)",
             .before = 1)
    
    out$geo_worst_vs_best <- geo_change %>%
      count(worst_unit_catchup, name = "n_countries") %>%
      mutate(pct_countries = round(100 * n_countries /
                                     sum(n_countries), 1),
             metric = "worst-vs-best unit catch-up",
             .before = 1)
  }
  
  if (!is.null(subgroup_change) && "gap_direction" %in% names(subgroup_change)) {
    out$subgroup_gaps <- subgroup_change %>%
      count(gap_direction, name = "n_countries") %>%
      mutate(pct_countries = round(100 * n_countries /
                                     sum(n_countries), 1),
             metric = "subgroup absolute gap",
             .before = 1)
    
    out$subgroup_pro_poor <- subgroup_change %>%
      count(pro_poor_progress, name = "n_countries") %>%
      mutate(pct_countries = round(100 * n_countries /
                                     sum(n_countries), 1),
             metric = "low-anchor vs high-anchor change",
             .before = 1)
    
    # Headline numbers
    out$subgroup_headline <- subgroup_change %>%
      summarise(
        n_countries_total  = n(),
        median_first_gap   = median(first_abs_gap,      na.rm = TRUE),
        median_latest_gap  = median(latest_abs_gap,     na.rm = TRUE),
        median_gap_change  = median(abs_gap_change_pp,  na.rm = TRUE),
        mean_gap_change    = mean(abs_gap_change_pp,    na.rm = TRUE),
        n_narrowing        = sum(gap_direction == "narrowing",
                                 na.rm = TRUE),
        n_widening         = sum(gap_direction == "widening",
                                 na.rm = TRUE),
        n_stable           = sum(gap_direction == "stable",
                                 na.rm = TRUE)
      )
  }
  
  out
}


# ================================================================
# 6. Main entry point - temporal analysis
# ----------------------------------------------------------------
# run_temporal_inequality_analysis(
#   surveys        = named list of DHS data.frames,
#   scope          = "between_countries" | "within_country",
#   geo_level      = "country" | "adm1" | "adm2" | "cluster",
#   subgroup       = NULL or name in subgroup_definitions,
#   min_surveys    = 2  (drop countries with <N surveys),
#   countries      = NULL or character vector of CountryName,
#   output_folder  = directory to write .xlsx into,
#   file_prefix    = prefix for the filename
# )
# ================================================================
run_temporal_inequality_analysis <- function(
    surveys,
    scope          = c("between_countries", "within_country"),
    geo_level      = c("country", "adm1", "adm2", "cluster"),
    subgroup       = NULL,
    min_surveys    = 2,
    countries      = NULL,
    output_folder  = ".",
    file_prefix    = "mcv1_temporal"
) {
  
  scope     <- match.arg(scope)
  geo_level <- match.arg(geo_level)
  
  # Coherence rules
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
      stop("Unknown subgroup '", subgroup, "'. Options: ",
           paste(names(subgroup_definitions), collapse = ", "))
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
  
  # Restrict to countries with >= min_surveys
  country_counts <- purrr::imap_dfr(surveys, function(df, nm) {
    if (!is.data.frame(df) || !"CountryName" %in% names(df)) return(NULL)
    tibble(survey_name = nm,
           CountryName = as.character(df$CountryName[1]))
  }) %>%
    count(CountryName, name = "n_surveys") %>%
    filter(n_surveys >= min_surveys)
  
  surveys <- keep(surveys, function(s) {
    is.data.frame(s) &&
      as.character(s$CountryName[1]) %in% country_counts$CountryName
  })
  
  message("Temporal run: scope=", scope, " | geo_level=", geo_level,
          " | subgroup=", ifelse(is.null(subgroup), "none", subgroup),
          " | countries with >=", min_surveys, " surveys: ",
          nrow(country_counts))
  
  if (length(surveys) == 0) {
    message("  No countries with enough surveys; nothing to do.")
    return(invisible(NULL))
  }
  
  # --- 1. Per-survey coverage at requested geo_level (overall)
  coverage_overall <- purrr::imap_dfr(surveys, function(df, nm) {
    if (!is.data.frame(df)) return(NULL)
    .coverage_one_survey(df, nm, geo_vars, subgroup_def = NULL)
  })
  
  # --- 2. Per-survey coverage disaggregated by subgroup (if given)
  coverage_sub <- NULL
  if (!is.null(sg_def)) {
    coverage_sub <- purrr::imap_dfr(surveys, function(df, nm) {
      if (!is.data.frame(df)) return(NULL)
      .coverage_one_survey(df, nm, geo_vars, subgroup_def = sg_def)
    })
  }
  
  # --- 3. Country-level coverage trajectory (national mean trend)
  #
  # For scope = between_countries this uses coverage_overall directly
  # (one row per country-survey). For scope = within_country we first
  # pool across sub-national units to a country mean per survey.
  if (scope == "between_countries") {
    country_mean_ts <- coverage_overall
  } else {
    country_mean_ts <- coverage_overall %>%
      group_by(CountryName, SurveyYear, country_phase, disagg) %>%
      summarise(
        coverage_pct   = weighted.mean(coverage_pct, w = n_total_12_23,
                                       na.rm = TRUE),
        se_pct         = sqrt(mean(se_pct^2, na.rm = TRUE)),  # rough pool
        n_total_12_23  = sum(n_total_12_23, na.rm = TRUE),
        n_vacc         = sum(n_vacc,        na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  country_trends <- .country_trajectory(country_mean_ts)
  
  # --- 4. Geographical dispersion trajectory (only for within-country)
  geo_disp_ts      <- NULL
  geo_disp_change  <- NULL
  if (scope == "within_country") {
    gdt <- .geo_dispersion_trajectory(coverage_overall)
    geo_disp_ts     <- gdt$ts
    geo_disp_change <- gdt$change
  }
  
  # --- 5. Subgroup trajectory (if subgroup requested)
  subgroup_ts     <- NULL
  subgroup_change <- NULL
  if (!is.null(coverage_sub)) {
    sgt <- .subgroup_trajectory(coverage_sub, sg_def)
    subgroup_ts     <- sgt$ts
    subgroup_change <- sgt$change
  }
  
  # --- 6. Cross-country convergence summary
  convergence <- .convergence_summary(
    geo_change      = geo_disp_change,
    subgroup_change = subgroup_change
  )
  
  # --- 7. Metadata
  metadata <- tibble::tibble(
    parameter = c("scope", "geo_level", "subgroup",
                  "min_surveys_per_country",
                  "n_countries_retained", "n_surveys_used",
                  "countries_filter", "timestamp",
                  "low_anchor", "high_anchor"),
    value = c(
      scope, geo_level,
      ifelse(is.null(subgroup), "none", subgroup),
      as.character(min_surveys),
      as.character(nrow(country_counts)),
      as.character(length(unique(coverage_overall$survey_name))),
      ifelse(is.null(countries), "all",
             paste(countries, collapse = "; ")),
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      ifelse(is.null(sg_def), NA, sg_def$low_group),
      ifelse(is.null(sg_def), NA, sg_def$high_group)
    )
  )
  
  # --- 8. Assemble sheets and save
  sheets <- list(metadata = metadata)
  sheets$coverage_ts     <- coverage_overall
  sheets$country_trends  <- country_trends
  if (!is.null(geo_disp_ts))     sheets$geo_dispersion_ts     <- geo_disp_ts
  if (!is.null(geo_disp_change)) sheets$geo_dispersion_change <- geo_disp_change
  if (!is.null(coverage_sub))    sheets$coverage_ts_subgroup  <- coverage_sub
  if (!is.null(subgroup_ts))     sheets$subgroup_trends_ts    <- subgroup_ts
  if (!is.null(subgroup_change)) sheets$subgroup_change       <- subgroup_change
  
  # Convergence summary tables (one sheet each where present)
  if (!is.null(convergence$geo_dispersion))
    sheets$convergence_geo         <- convergence$geo_dispersion
  if (!is.null(convergence$geo_worst_vs_best))
    sheets$convergence_worst_best  <- convergence$geo_worst_vs_best
  if (!is.null(convergence$subgroup_gaps))
    sheets$convergence_subgroup    <- convergence$subgroup_gaps
  if (!is.null(convergence$subgroup_pro_poor))
    sheets$convergence_pro_poor    <- convergence$subgroup_pro_poor
  if (!is.null(convergence$subgroup_headline))
    sheets$convergence_headline    <- convergence$subgroup_headline
  
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
# 7. Batch runner - all combinations
# ================================================================
run_all_temporal_analyses <- function(
    surveys,
    output_folder        = ".",
    subgroups            = names(subgroup_definitions),
    include_no_subgroup  = TRUE,
    min_surveys          = 2
) {
  
  specs <- list()
  
  # Between-countries trend (geo = country)
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
  
  message("Temporal batch run: ", length(specs), " combinations.")
  
  results <- purrr::imap(specs, function(s, i) {
    message(sprintf("[%d/%d] scope=%s geo=%s sub=%s",
                    i, length(specs), s$scope, s$geo_level,
                    ifelse(is.null(s$subgroup), "none", s$subgroup)))
    tryCatch(
      run_temporal_inequality_analysis(
        surveys       = surveys,
        scope         = s$scope,
        geo_level     = s$geo_level,
        subgroup      = s$subgroup,
        min_surveys   = min_surveys,
        output_folder = output_folder
      ),
      error = function(e) { message("    FAILED: ", e$message); NULL }
    )
  })
  
  invisible(results)
}


# ================================================================
# 8. Example usage (commented out)
# ================================================================
#
# source("10_descriptive_inequality_analysis.R")   # needed for helpers
#
# final_DHS_data <- readRDS("…/final_DHS_data.rds")
#
# output_folder <- "…/TABLES/mcv1_temporal"
#
# # (a) Between-country national trend, overall
# run_temporal_inequality_analysis(
#   surveys       = final_DHS_data,
#   scope         = "between_countries",
#   geo_level     = "country",
#   subgroup      = NULL,
#   output_folder = output_folder
# )
#
# # (b) Within-country ADM1 dispersion trajectory, by wealth quintile
# run_temporal_inequality_analysis(
#   surveys       = final_DHS_data,
#   scope         = "within_country",
#   geo_level     = "adm1",
#   subgroup      = "wealth_index_quintile",
#   output_folder = output_folder
# )
#
# # (c) Nigeria only, cluster-level, by maternal education
# run_temporal_inequality_analysis(
#   surveys       = final_DHS_data,
#   scope         = "within_country",
#   geo_level     = "cluster",
#   subgroup      = "respondent_edu_level",
#   countries     = "Nigeria",
#   output_folder = output_folder
# )
#
# # (d) Everything
# run_all_temporal_analyses(
#   surveys       = final_DHS_data,
#   output_folder = output_folder
# )
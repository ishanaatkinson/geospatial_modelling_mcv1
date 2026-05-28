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

# RUN DATA AND FUNCTIONS FOR INITIAL COVERAGE ANALYSIS . R

run_type <- "full"

source("~/YEAR 1/PhD/MEASLES PROJECT/CODE/Vaccination coverage/10.0_data_and_functions_for_initial_coverage_analysis.R", local = TRUE)

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
  
  
  
  
  # Posterior sampling failed: Package 'sn' is required to proceed, but is not installed. Please install. 
  # Posterior sampling not available — falling back to point estimate.
  
  
  
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
  
  # ═══════════════════════════════════════════════════════════════════════
  # ADDITION 1
  # PLACE THIS IMMEDIATELY AFTER:
  #   inla_sp <- NULL; inla_fu <- NULL; inla_bb <- NULL
  # ═══════════════════════════════════════════════════════════════════════
  
  m5_result <- NULL
  m6_result <- NULL
  boundary_probs <- NULL
  pm_m5 <- list(r = NA, r2 = NA, mae = NA, rmse = NA)
  pm_m6 <- list(r = NA, r2 = NA, mae = NA, rmse = NA)
  cv_result <- NULL
  cv_summary <- NULL
  cv_predictions <- NULL
  
  
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
      
      # version <- "_(v2)"
      # 
      # output_folder <- file.path(base_path, paste0("TABLES/Statical_Model_Results", version))
      # if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
      # 
      # rg <- max(ext / 3, 20)
      # max_edge <- c(rg / 5, rg)
      # offset_value <- c(rg / 2, rg * 1.5)
      # cutoff_value <- max(rg / 20, 5)
      # prior_on_range <- c(rg/2, 0.5)
      # prior_on_variance <- c(2, 0.05)
      # nugget_effect <- "c(2, 0.05)"
      
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
      # rg <- max(ext / 5, 10)
      # max_edge <- c(rg / 8, rg/2)
      # offset_value <- c(rg / 3, rg)
      # cutoff_value <- max(rg / 30, 2)
      # prior_on_range <- c(rg, 0.5)
      # prior_on_variance <- c(1, 0.01)
      # nugget_effect <- "c(1, 0.01)"
      
      # V6 (varying prediction of country surfaces for mcv1 and covariates)
      
      version <- "_(v6)"
      
      
      # RE RUN WITH THESE PARAMETERS
      
      # uses posterior sample method for mcv1 and covariate surface mapping
      
      output_folder <- file.path(base_path, paste0("TABLES/Statical_Model_Results", version))
      if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
      
      rg <- max(ext / 3, 20)
      max_edge <- c(rg / 5, rg)
      offset_value <- c(rg / 2, rg * 1.5)
      cutoff_value <- max(rg / 20, 5)
      prior_on_range <- c(rg, 0.5)
      prior_on_variance <- c(1, 0.01)
      nugget_effect <- "c(1, 0.01)"
      
      
      # version <- "_(v7)"
      # 
      # # uses prediction stack inla method for mcv1 and covariate surface mapping
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
      
      
      
      
      
      
      prior_info <- list(
        version = version,
        
        range = list(
          rg = rg,
          prior = prior_on_range,
          statement = paste0("P(range < ", prior_on_range[1], ") = ", prior_on_range[2])
        ),
        
        sigma = list(
          prior = prior_on_variance,
          statement = paste0("P(σ > ", prior_on_variance[1], ") = ", prior_on_variance[2])
        ),
        
        nugget = list(
          prior = nugget_effect,
          statement = paste0("PC prior: ", nugget_effect)
        )
      )
      
      
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
        
        
        # ═══════════════════════════════════════════════════════════════════════
        # ADDITION 2
        # PLACE THIS INSIDE THE MESH BLOCK, IMMEDIATELY AFTER MODEL 4 ENDS:
        #
        #   if (!is.null(inla_bb)) {
        #     ...
        #     sub_m$resid_bb <- sub_m$coverage - sub_m$pred_bb
        #   }
        #
        # and BEFORE:
        #   } # end if mesh
        # ═══════════════════════════════════════════════════════════════════════
        
        # MODEL 5-6: Zero-inflated binomial models
        cat("\n  Models 5-6: Zero-inflated binomial models...\n")
        
        boundary_probs <- estimate_excess_zero_prob(sub_m, "pred_glm")
        
        cat("    Observed zeros:", round(boundary_probs$observed_zero * 100, 1), "%\n")
        cat("    Expected zeros under binomial:", round(boundary_probs$expected_zero * 100, 1), "%\n")
        cat("    Excess zero inflation pi0:", round(boundary_probs$pi0_excess * 100, 2), "%\n")
        cat("    Observed ones:", round(boundary_probs$observed_one * 100, 1), "%\n")
        cat("    Note: ones remain in the binomial process; no one-inflation is modelled.\n")
        
        m5_result <- fit_model5_zib_fixed(
          sub_m           = sub_m,
          vars_final      = vars_final,
          mesh            = mesh,
          spde            = spde,
          A               = A,
          stk_template    = stk,
          iid_spec        = iid_spec,
          cov_terms       = cov_terms,
          boundary_probs  = boundary_probs,
          run_inla_fn     = run_inla
        )
        
        if (!is.null(m5_result)) {
          sub_m$pred_m5  <- m5_result$predictions
          sub_m$resid_m5 <- sub_m$coverage - sub_m$pred_m5
        }
        
        m6_result <- fit_model6_zib_estimated(
          sub_m           = sub_m,
          vars_final      = vars_final,
          mesh            = mesh,
          spde            = spde,
          A               = A,
          stk_template    = stk,
          iid_spec        = iid_spec,
          cov_terms       = cov_terms,
          boundary_probs  = boundary_probs,
          run_inla_fn     = run_inla
        )
        
        if (!is.null(m6_result)) {
          sub_m$pred_m6  <- m6_result$predictions
          sub_m$resid_m6 <- sub_m$coverage - sub_m$pred_m6
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
  
  # ═══════════════════════════════════════════════════════════════════════
  # ADDITION 3
  # PLACE THIS IMMEDIATELY AFTER YOUR EXISTING MODEL COMPARISON LINE:
  #   cat("\n"); print(t(mc), quote = FALSE)
  #
  # This updates pred_table, mc, best_model, and runs CV before plotting.
  # ═══════════════════════════════════════════════════════════════════════
  
  pm_m5 <- if ("pred_m5" %in% names(sub_m)) calc_pred_metrics(sub_m$coverage, sub_m$pred_m5)
  else list(r = NA, r2 = NA, mae = NA, rmse = NA)
  
  pm_m6 <- if ("pred_m6" %in% names(sub_m)) calc_pred_metrics(sub_m$coverage, sub_m$pred_m6)
  else list(r = NA, r2 = NA, mae = NA, rmse = NA)
  
  pred_table <- bind_rows(
    pred_table,
    tibble(
      Model = c("ZIB-Fixed", "ZIB-Estimated"),
      r     = c(pm_m5$r, pm_m6$r),
      R2    = c(pm_m5$r2, pm_m6$r2),
      MAE   = c(pm_m5$mae, pm_m6$mae)
    )
  )
  
  mc <- mc %>%
    mutate(
      m5_dic   = if (!is.null(m5_result)) round(m5_result$dic, 1) else NA,
      m5_waic  = if (!is.null(m5_result)) round(m5_result$waic, 1) else NA,
      m5_r     = pm_m5$r,
      m5_r2    = pm_m5$r2,
      m5_mae   = pm_m5$mae,
      m5_pi0   = if (!is.null(m5_result)) round(m5_result$pi0, 4) else NA,
      m5_pi1   = 0,
      m5_range = if (!is.null(m5_result)) round(m5_result$spatial_range, 1) else NA,
      m5_var   = if (!is.null(m5_result)) round(m5_result$spatial_var, 4) else NA,
      
      m6_dic       = if (!is.null(m6_result)) round(m6_result$dic, 1) else NA,
      m6_waic      = if (!is.null(m6_result)) round(m6_result$waic, 1) else NA,
      m6_r         = pm_m6$r,
      m6_r2        = pm_m6$r2,
      m6_mae       = pm_m6$mae,
      m6_mean_pi0  = if (!is.null(m6_result)) round(mean(m6_result$pred_pi0), 4) else NA,
      m6_mean_pi1  = 0,
      m6_range     = if (!is.null(m6_result)) round(m6_result$spatial_range, 1) else NA,
      m6_var       = if (!is.null(m6_result)) round(m6_result$spatial_var, 4) else NA
    )
  
  all_waics <- c(
    Full          = mc$full_waic,
    BetaBinom     = mc$bb_waic,
    ZIB_Fixed     = mc$m5_waic,
    ZIB_Estimated = mc$m6_waic
  )
  all_waics <- all_waics[!is.na(all_waics)]
  
  mc$best_model <- if (length(all_waics) > 0) {
    best_name <- names(which.min(all_waics))
    switch(
      best_name,
      Full          = "Full INLA",
      BetaBinom     = "Beta-binomial",
      ZIB_Fixed     = "ZIB-Fixed",
      ZIB_Estimated = "ZIB-Estimated",
      "Full INLA"
    )
  } else if (!is.null(fit_glm)) {
    "GLM"
  } else {
    "None"
  }
  
  all_comp[[cname]] <- mc
  
  cat("\n  Updated model comparison including Models 5 and 6:\n")
  print(t(mc), quote = FALSE)
  
  if (!is.null(mesh) && !is.null(spde) && inla_available && nrow(sub_m) >= 50) {
    
    cat("\n--- Cross-validation: Models 2-6 ---\n")
    
    cv_result <- tryCatch(
      run_cv_models_2_to_6_utazi_style(
        sub_m = sub_m,
        vars_final = vars_final,
        k = 10,
        schemes = c("random", "spatial"),
        seed = 123,
        boundary_probs = boundary_probs,
        n_crps_samples = 500
      ),
      error = function(e) {
        cat("  Cross-validation failed:", e$message, "\n")
        NULL
      }
    )
    
    if (!is.null(cv_result)) {
      cv_summary <- cv_result$summary
      cv_predictions <- cv_result$predictions
      
      cat("\n  Cross-validation summary:\n")
      print(cv_summary %>% mutate(across(where(is.numeric), ~ round(.x, 4))),
            n = Inf)
      
      get_cv_val <- function(model_name, scheme_name, metric_name) {
        vv <- cv_summary %>%
          dplyr::filter(model == model_name, scheme == scheme_name) %>%
          dplyr::pull({{ metric_name }})
        if (length(vv) == 0) NA_real_ else vv[1]
      }
      
      mc <- mc %>%
        mutate(
          cv_random_sp_mae    = get_cv_val("Model 2: Spatial-only", "random", mae),
          cv_random_full_mae  = get_cv_val("Model 3: Full INLA", "random", mae),
          cv_random_bb_mae    = get_cv_val("Model 4: Beta-binomial", "random", mae),
          cv_random_m5_mae    = get_cv_val("Model 5: ZIB-Fixed", "random", mae),
          cv_random_m6_mae    = get_cv_val("Model 6: ZIB-Estimated", "random", mae),
          cv_spatial_sp_mae   = get_cv_val("Model 2: Spatial-only", "spatial", mae),
          cv_spatial_full_mae = get_cv_val("Model 3: Full INLA", "spatial", mae),
          cv_spatial_bb_mae   = get_cv_val("Model 4: Beta-binomial", "spatial", mae),
          cv_spatial_m5_mae   = get_cv_val("Model 5: ZIB-Fixed", "spatial", mae),
          cv_spatial_m6_mae   = get_cv_val("Model 6: ZIB-Estimated", "spatial", mae)
        )
      
      all_comp[[cname]] <- mc
    }
  }
  
  
  
  
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
  
  
  # ═══════════════════════════════════════════════════════════════════════
  # ADDITION 4
  # PLACE THIS IMMEDIATELY AFTER YOUR EXISTING BEST-PREDICTION BLOCK:
  #
  #   if ("pred_full" %in% names(sub_m)) { ... }
  #   else if ("pred_glm" %in% names(sub_m)) { ... }
  #   else { ... }
  #
  # and BEFORE:
  #   sub_sf <- st_as_sf(...)
  # ═══════════════════════════════════════════════════════════════════════
  
  if (mc$best_model == "ZIB-Estimated" && "pred_m6" %in% names(sub_m)) {
    sub_m$pred_best  <- sub_m$pred_m6
    sub_m$resid_best <- sub_m$resid_m6
    bl <- "ZIB-Estimated (Zero-inflated binomial + Covariates + Matern GP + Nugget)"
  } else if (mc$best_model == "ZIB-Fixed" && "pred_m5" %in% names(sub_m)) {
    sub_m$pred_best  <- sub_m$pred_m5
    sub_m$resid_best <- sub_m$resid_m5
    bl <- "ZIB-Fixed (Zero-inflated binomial + Covariates + Matern GP + Nugget)"
  } else if (mc$best_model == "Beta-binomial" && "pred_bb" %in% names(sub_m)) {
    sub_m$pred_best  <- sub_m$pred_bb
    sub_m$resid_best <- sub_m$resid_bb
    bl <- "Beta-binomial INLA (Covariates + Matern GP + Nugget)"
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
  
  country_dir <- file.path(output_folder, country_slug)
  
  dir.create(country_dir, recursive = TRUE, showWarnings = FALSE)
  
  
  pdf_path <- file.path(country_dir, paste0(country_slug, "_geostatistical_models_unweighted", version, ".pdf"))
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
  
  # ═══════════════════════════════════════════════════════════════════════
  # ADDITION 5
  # PLACE THIS INSIDE PAGE 2 SUMMARY, AFTER THE EXISTING BETA-BINOMIAL
  # MODEL FIT summary_lines block and BEFORE:
  #   # Prediction metrics
  # ═══════════════════════════════════════════════════════════════════════
  
  if (!is.null(m5_result))
    summary_lines <- c(
      summary_lines,
      paste0(
        "ZIB-Fixed: DIC=", mc$m5_dic,
        " | WAIC=", mc$m5_waic,
        " | pi0=", mc$m5_pi0,
        " | Range=", mc$m5_range, " km",
        " | Variance=", mc$m5_var
      )
    )
  
  if (!is.null(m6_result))
    summary_lines <- c(
      summary_lines,
      paste0(
        "ZIB-Estimated: DIC=", mc$m6_dic,
        " | WAIC=", mc$m6_waic,
        " | mean pi0=", mc$m6_mean_pi0,
        " | Range=", mc$m6_range, " km",
        " | Variance=", mc$m6_var
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
  
  
  # ═══════════════════════════════════════════════════════════════════════
  # ADDITION 6
  # PLACE THIS AFTER THE EXISTING PREDICTION METRICS summary_lines BLOCK
  # and BEFORE:
  #   # Model fit assessment
  # ═══════════════════════════════════════════════════════════════════════
  
  summary_lines <- c(
    summary_lines,
    
    sprintf(
      "  %-30s  r=%s  R²=%s  MAE=%s",
      "ZIB-Fixed:",
      ifelse(is.na(pm_m5$r), " N/A", sprintf("%6.3f", pm_m5$r)),
      ifelse(is.na(pm_m5$r2), " N/A", sprintf("%6.3f", pm_m5$r2)),
      ifelse(is.na(pm_m5$mae), "N/A", sprintf("%5.3f", pm_m5$mae))
    ),
    
    sprintf(
      "  %-30s  r=%s  R²=%s  MAE=%s",
      "ZIB-Estimated:",
      ifelse(is.na(pm_m6$r), " N/A", sprintf("%6.3f", pm_m6$r)),
      ifelse(is.na(pm_m6$r2), " N/A", sprintf("%6.3f", pm_m6$r2)),
      ifelse(is.na(pm_m6$mae), "N/A", sprintf("%5.3f", pm_m6$mae))
    )
  )
  
  if (!is.null(cv_summary)) {
    
    cv_best_random <- cv_summary %>%
      filter(scheme == "random") %>%
      arrange(mae) %>%
      slice(1)
    
    cv_best_spatial <- cv_summary %>%
      filter(scheme == "spatial") %>%
      arrange(mae) %>%
      slice(1)
    
    summary_lines <- c(
      summary_lines,
      "",
      "----- CROSS-VALIDATION: MODELS 2-6 -----",
      "  Lower MAE/RMSE/CRPS is better; higher R² is better."
    )
    
    if (nrow(cv_best_random) == 1) {
      summary_lines <- c(
        summary_lines,
        paste0(
          "Best random CV MAE: ", cv_best_random$model,
          " | MAE=", round(cv_best_random$mae, 3),
          " | R²=", round(cv_best_random$r2, 3),
          " | CRPS=", round(cv_best_random$crps, 3)
        )
      )
    }
    
    if (nrow(cv_best_spatial) == 1) {
      summary_lines <- c(
        summary_lines,
        paste0(
          "Best spatial CV MAE: ", cv_best_spatial$model,
          " | MAE=", round(cv_best_spatial$mae, 3),
          " | R²=", round(cv_best_spatial$r2, 3),
          " | CRPS=", round(cv_best_spatial$crps, 3)
        )
      )
    }
  }
  
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
  
  # After the forest plot section, add:
  if (!is.null(feff)) {
    effect_table <- feff %>%
      filter(term != "intercept") %>%
      mutate(
        OR_formatted = sprintf("%.3f (%.3f – %.3f)", OR, OR_lo, OR_hi),
        beta_formatted = sprintf("%.3f (%.3f)", mean, sd),
        sig = ifelse(sig == "*", "Yes", "No")
      ) %>%
      dplyr::select(
        Covariate = term,
        `β (SD)` = beta_formatted,
        `OR (95% CrI)` = OR_formatted,
        Significant = sig
      )
    
    tg <- gridExtra::tableGrob(
      effect_table, rows = NULL,
      theme = gridExtra::ttheme_minimal(
        core    = list(fg_params = list(fontsize = 8)),
        colhead = list(fg_params = list(fontsize = 9, fontface = "bold"))
      )
    )
    
    grid::grid.newpage()
    grid::grid.text(
      paste(cname, "— Fixed Effect Sizes (Full INLA Model)"),
      x = 0.5, y = 0.97, gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
    grid::grid.draw(tg)
  }
  
  
  # ----- PAGE 4: Plot prior and posterior densities  --------------------------------------------------
  
  # The PC priors for the SPDE are defined on a transformed scale internally, 
  # and the exact prior density on the "natural" (range, sigma) scale involves 
  # a Jacobian. The code above gives the approximate exponential form which is 
  # the standard way PC priors are visualised in the literature
  
  if (!is.null(inla_fu)) {
    
    # --- Spatial range ---
    spr <- INLA::inla.spde2.result(inla_fu, "spatial", spde, do.transf = TRUE)
    
    # Posterior marginal for range
    range_post <- spr$marginals.range.nominal[[1]]
    
    # PC prior for range: exponential on 1/range
    # P(range < rg) = 0.5 → lambda = -log(0.5) / rg
    lambda_range <- -log(0.5) / rg
    range_x <- seq(0.01, max(range_post[, 1]) * 1.2, length.out = 500)
    range_prior_dens <- lambda_range * exp(-lambda_range * range_x)
    
    range_df <- rbind(
      data.frame(x = range_post[, 1], y = range_post[, 2], 
                 Distribution = "Posterior"),
      data.frame(x = range_x, y = range_prior_dens, 
                 Distribution = "Prior")
    )
    
    get_post_mean <- function(marg) {
      INLA::inla.emarginal(function(x) x, marg)
    }
    
    
    range_post_mean <- get_post_mean(range_post)
    
    p_range <- ggplot(range_df, aes(x = x, y = y, colour = Distribution)) +
      geom_line(linewidth = 1) +
      
      geom_vline(xintercept = range_post_mean,
                 linetype = "dashed",
                 linewidth = 0.8,
                 colour = "steelblue") +
      
      scale_colour_manual(values = c("Prior" = "grey50", "Posterior" = "steelblue")) +
      
      labs(
        title = paste(cname, "— Spatial range (PC prior, version", version, ")"),
        subtitle = paste0(
          prior_info$range$statement,
          " | Posterior mean = ",
          round(range_post_mean, 1),
          " km"
        ),
        x = "Range (km)",
        y = "Density"
      ) +
      theme_minimal(base_size = 10)
    
    # --- Spatial sigma (standard deviation) ---
    sigma_post <- spr$marginals.variance.nominal[[1]]
    # Transform to SD
    sigma_sd_post <- INLA::inla.tmarginal(sqrt, sigma_post)
    
    # PC prior for sigma: P(sigma > 1) = 0.01
    # → lambda_sigma = -log(0.01) / 1
    lambda_sigma <- -log(prior_on_variance[2]) / prior_on_variance[1]
    sigma_x <- seq(0.01, max(sigma_sd_post[, 1]) * 1.2, length.out = 500)
    sigma_prior_dens <- lambda_sigma * exp(-lambda_sigma * sigma_x)
    
    sigma_df <- rbind(
      data.frame(x = sigma_sd_post[, 1], y = sigma_sd_post[, 2],
                 Distribution = "Posterior"),
      data.frame(x = sigma_x, y = sigma_prior_dens,
                 Distribution = "Prior")
    )
    
    sigma_post_mean <- get_post_mean(sigma_sd_post)
    
    p_sigma <- ggplot(sigma_df, aes(x = x, y = y, colour = Distribution)) +
      geom_line(linewidth = 1) +
      
      geom_vline(xintercept = sigma_post_mean,
                 linetype = "dashed",
                 linewidth = 0.8,
                 colour = "steelblue") +
      
      scale_colour_manual(values = c("Prior" = "grey50", "Posterior" = "steelblue")) +
      
      labs(
        title = paste(cname, "— Spatial marginal SD (version", version, ")"),
        subtitle = paste0(
          prior_info$sigma$statement,
          " | Posterior mean = ",
          round(sigma_post_mean, 3)
        ),
        x = "σ",
        y = "Density"
      ) +
      theme_minimal(base_size = 10)
    
    # --- Nugget SD ---
    # The nugget is f(cluster_id, model='iid') with PC prior on precision
    nugget_prec_post <- inla_fu$marginals.hyperpar[["Precision for cluster_id"]]
    # Transform precision → SD: sd = 1/sqrt(prec)
    nugget_sd_post <- INLA::inla.tmarginal(function(x) 1/sqrt(x), nugget_prec_post)
    
    # PC prior for nugget: param = c(u, alpha) means P(sd > u) = alpha
    # Parse from nugget_effect string or use the values directly
    nugget_u <- 1  # from c(1, 0.01)
    nugget_alpha <- 0.01
    lambda_nugget <- -log(nugget_alpha) / nugget_u
    nugget_x <- seq(0.01, max(nugget_sd_post[, 1]) * 1.2, length.out = 500)
    nugget_prior_dens <- lambda_nugget * exp(-lambda_nugget * nugget_x)
    
    nugget_df <- rbind(
      data.frame(x = nugget_sd_post[, 1], y = nugget_sd_post[, 2],
                 Distribution = "Posterior"),
      data.frame(x = nugget_x, y = nugget_prior_dens,
                 Distribution = "Prior")
    )
    
    nugget_post_mean <- get_post_mean(nugget_sd_post)
    
    p_nugget <- ggplot(nugget_df, aes(x = x, y = y, colour = Distribution)) +
      geom_line(linewidth = 1) +
      
      geom_vline(xintercept = nugget_post_mean,
                 linetype = "dashed",
                 linewidth = 0.8,
                 colour = "steelblue") +
      
      scale_colour_manual(values = c("Prior" = "grey50", "Posterior" = "steelblue")) +
      
      labs(
        title = paste(cname, "— Nugget SD (iid, version", version, ")"),
        subtitle = paste0(
          prior_info$nugget$statement,
          " | Posterior mean = ",
          round(nugget_post_mean, 3)
        ),
        x = "σ_nugget",
        y = "Density"
      ) +
      theme_minimal(base_size = 10)
    # Print as a combined page
    print(p_range / p_sigma / p_nugget +
            plot_annotation(title = paste(cname, "— Prior vs Posterior distributions")))
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
  r_val <- cor(sub_m$pred_best, sub_m$coverage, use = "complete.obs")
  
  p4 <- ggplot(sub_m, aes(x = pred_best, y = coverage)) +
    
    geom_point(aes(size = n_children),
               alpha = 0.3,
               colour = "steelblue") +
    
    # identity line (perfect prediction)
    geom_abline(intercept = 0, slope = 1,
                linetype = "dashed",
                colour = "red") +
    
    # fitted gradient line (smooth relationship)
    geom_smooth(method = "lm",
                se = FALSE,
                colour = "#7a78a9",
                linewidth = 1) +
    
    scale_size_continuous(name = "N children", range = c(0.5, 4)) +
    
    annotate(
      "text",
      x = 0.05, y = 0.95,
      label = paste0("r = ", round(r_val, 3)),
      hjust = 0,
      vjust = 1,
      size = 4
    ) +
    
    labs(
      x = "Predicted coverage",
      y = "Observed coverage",
      title = paste(cname, "Observed vs Predicted")
    ) +
    
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
  
  
  # ═══════════════════════════════════════════════════════════════════════
  # ADDITION 7
  # PLACE THIS AFTER THE EXISTING 6-PANEL DIAGNOSTIC PAGE IS PRINTED:
  #   print(combined_page)
  # and BEFORE:
  #   # PAGE 6: Full-country prediction...
  # ═══════════════════════════════════════════════════════════════════════
  
  model_pred_long <- sub_m %>%
    dplyr::select(
      coverage,
      n_children,
      any_of(c("pred_glm", "pred_sp", "pred_full", "pred_bb", "pred_m5", "pred_m6"))
    ) %>%
    tidyr::pivot_longer(
      cols = starts_with("pred_"),
      names_to = "model",
      values_to = "pred"
    ) %>%
    mutate(
      model = dplyr::recode(
        model,
        pred_glm = "Model 1: GLM",
        pred_sp = "Model 2: Spatial-only",
        pred_full = "Model 3: Full INLA",
        pred_bb = "Model 4: Beta-binomial",
        pred_m5 = "Model 5: ZIB-Fixed",
        pred_m6 = "Model 6: ZIB-Estimated"
      )
    )
  
  if (nrow(model_pred_long) > 0) {
    p_all_scatter <- ggplot(model_pred_long, aes(x = pred, y = coverage)) +
      geom_point(aes(size = n_children), alpha = 0.25, colour = "steelblue") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red") +
      geom_smooth(method = "lm", se = FALSE, colour = "#7a78a9", linewidth = 0.8) +
      facet_wrap(~ model, ncol = 2) +
      scale_size_continuous(name = "N children", range = c(0.4, 2.5)) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(
        title = paste(cname, "- Observed vs Predicted Coverage for Models 1-6"),
        x = "Predicted coverage",
        y = "Observed coverage"
      ) +
      theme_minimal(base_size = 10) +
      theme(legend.position = "bottom")
    print(p_all_scatter)
  }
  
  model_resid_long <- sub_m %>%
    dplyr::select(
      any_of(c("resid_glm", "resid_sp", "resid_full", "resid_bb", "resid_m5", "resid_m6"))
    ) %>%
    tidyr::pivot_longer(
      cols = everything(),
      names_to = "model",
      values_to = "resid"
    ) %>%
    mutate(
      model = dplyr::recode(
        model,
        resid_glm = "Model 1: GLM",
        resid_sp = "Model 2: Spatial-only",
        resid_full = "Model 3: Full INLA",
        resid_bb = "Model 4: Beta-binomial",
        resid_m5 = "Model 5: ZIB-Fixed",
        resid_m6 = "Model 6: ZIB-Estimated"
      )
    )
  
  if (nrow(model_resid_long) > 0) {
    p_all_resid <- ggplot(model_resid_long, aes(x = resid)) +
      geom_histogram(bins = 40, fill = "grey60", colour = "white") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
      facet_wrap(~ model, ncol = 2) +
      labs(
        title = paste(cname, "- Residual Distributions for Models 1-6"),
        x = "Residual (observed - predicted)",
        y = "Clusters"
      ) +
      theme_minimal(base_size = 10)
    print(p_all_resid)
  }
  
  if (!is.null(cv_result)) {
    if (!is.null(cv_result$plot)) {
      print(cv_result$plot +
              labs(
                title = paste(cname, "- Cross-validated MAE across Models 2-6"),
                subtitle = "Lower is better"
              ))
    }
    
    if (!is.null(cv_result$r2_plot)) {
      print(cv_result$r2_plot +
              labs(
                title = paste(cname, "- Cross-validated R² across Models 2-6"),
                subtitle = "Higher is better"
              ))
    }
    
    cv_rmse_plot <- ggplot(
      cv_summary,
      aes(x = reorder(model, rmse), y = rmse, fill = scheme)
    ) +
      geom_col(position = "dodge") +
      coord_flip() +
      labs(
        title = paste(cname, "- Cross-validated RMSE across Models 2-6"),
        subtitle = "Lower is better",
        x = NULL,
        y = "RMSE"
      ) +
      theme_minimal(base_size = 11)
    print(cv_rmse_plot)
  }
  # PAGE 6: Full-country prediction, spatial field, and covariate surfaces
  
  
  # NEW (for posterior_sample method — DEFAULT, no extra args needed):
  surface_diag <- make_country_surface_diagnostics(
    country_name = cname, adm1_sf = adm1_sf, sub_m = sub_m,
    vars_final = vars_final, inla_fu = inla_fu, mesh = mesh,
    spde = spde, lcrs = lcrs,
    method = "posterior_sample",
    n_posterior_samples = 500
  )
  
  # NEW (for prediction_stack method — must also pass stk, f3, iid_spec):
  # surface_diag <- make_country_surface_diagnostics(
  #   country_name = cname, adm1_sf = adm1_sf, sub_m = sub_m,
  #   vars_final = vars_final, inla_fu = inla_fu, mesh = mesh,
  #   spde = spde, lcrs = lcrs,
  #   stk_est = stk, f3 = f3, iid_spec = iid_spec,
  #   method = "prediction_stack"
  # )
  
  
  if (!is.null(surface_diag)) {
    
    # PAGE 6a: Predicted coverage + spatial field surfaces
    print(surface_diag$surface_plot)
    
    # PAGE 6b (NEW): Uncertainty maps
    if (length(surface_diag$uncertainty_plots) > 0) {
      
      # MCV1 uncertainty: SD and CrI width on one page
      mcv1_unc <- surface_diag$uncertainty_plots[
        names(surface_diag$uncertainty_plots) %in% c("mcv1_sd", "mcv1_ci_width")
      ]
      if (length(mcv1_unc) > 0) {
        print(
          wrap_plots(mcv1_unc, ncol = 2) +
            plot_annotation(
              title = paste(cname, "- MCV1 coverage prediction uncertainty")
            )
        )
      }
      
      # Covariate kriging uncertainty maps (6 per page)
      cov_unc <- surface_diag$uncertainty_plots[
        grep("^cov_", names(surface_diag$uncertainty_plots))
      ]
      if (length(cov_unc) > 0) {
        for (ii in seq(1, length(cov_unc), by = 6)) {
          print(
            wrap_plots(
              cov_unc[ii:min(ii + 5, length(cov_unc))],
              ncol = 3
            ) +
              plot_annotation(
                title = paste(cname, "- covariate interpolation uncertainty (kriging SD)")
              )
          )
        }
      }
    }
    
    # PAGE 6c: Covariate surfaces (now smooth from kriging)
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
    
    # PAGE 6d: Covariate histograms
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
    
    # PAGE 6e: Partial residual plots
    partial_residual_plots <- make_partial_residual_plots(
      sub_m = sub_m, vars_final = vars_final,
      feff = feff, glm_co = glm_co
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
    
  }  # end if (!is.null(surface_diag))
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
  
  
  
  
  
  
  
  
  
  # ═══════════════════════════════════════════════════════════════════════════════
  # SAVE ALL RAW DATA OUTPUTS — MCV1 GEOSTATISTICAL PIPELINE
  # ═══════════════════════════════════════════════════════════════════════════════
  #
  # PURPOSE
  #   Saves every piece of raw data produced by the pipeline so that maps,
  #   graphs, and cross-country analyses can be fully reproduced later without
  #   re-running any models.
  #
  # WHAT IS SAVED (no images, no plots — only data)
  #   ├── csv/                  — all tabular data as flat CSVs
  #   ├── rds/                  — R objects (models, lists, sf objects)
  #   ├── gpkg/                 — spatial vector data (GeoPackage)
  #   └── tif/                  — spatial raster data (GeoTIFF)
  #       ├── prediction/
  #       ├── uncertainty/
  #       ├── covariates/
  #       ├── spatial_field/
  #       └── residuals/
  #
  # REQUIRES (must exist in the environment before sourcing this file)
  #   cname, country_dir, country_slug, version
  #   sub_m, sub_sf, surface_diag, adm1_sf
  #   fit_glm, inla_sp, inla_fu, inla_bb
  #   m5_result, m6_result (can be NULL)
  #   pred_table, mc, biv, vars_final
  #   feff, glm_co, var_audit (can be NULL)
  #   cv_result, cv_summary, cv_predictions (can be NULL)
  #   mesh, spde, stk (can be NULL)
  #   lcrs  — local projected CRS string / EPSG code for raster output
  # ═══════════════════════════════════════════════════════════════════════════════
  
  
  # ═══════════════════════════════════════════════════════════════════════
  # 12.  EXCEL WORKBOOK — one definitive build of every sheet
  # ═══════════════════════════════════════════════════════════════════════
  
  sheets <- list()
  
  # --- Model summary --------------------------------------------------
  sheets[["Model_Summary"]] <- mc
  
  # --- Variables used -------------------------------------------------
  sheets[["Variables_Used"]] <- tibble(
    variable       = vars_final,
    bivariate_p    = biv$p_value[match(vars_final, biv$variable)],
    bivariate_OR   = biv$OR[match(vars_final, biv$variable)],
    bivariate_beta = biv$estimate[match(vars_final, biv$variable)]
  ) %>% mutate(across(where(is.numeric), ~ round(.x, 4)))
  
  # --- GLM coefficients -----------------------------------------------
  if (!is.null(glm_co))
    sheets[["GLM_Coefficients"]] <- glm_co %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))
  
  # --- INLA fixed effects ---------------------------------------------
  if (exists("feff") && !is.null(feff))
    sheets[["INLA_Fixed_Effects"]] <- feff %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))
  
  # --- Spatial parameters ---------------------------------------------
  if (!is.na(rp))
    sheets[["Spatial_Parameters"]] <- tibble(
      parameter = c("Spatial range (km)", "Spatial marginal variance",
                    "DIC (spatial-only)", "WAIC (spatial-only)",
                    "DIC (full model)",   "WAIC (full model)",
                    "DIC (beta-binomial)","WAIC (beta-binomial)"),
      value = c(round(rp, 1), round(sp, 4),
                mc$sp_dic,   mc$sp_waic,
                mc$full_dic, mc$full_waic,
                mc$bb_dic,   mc$bb_waic)
    )
  
  # --- Bivariate screening -------------------------------------------
  sheets[["Bivariate_Screening"]] <- biv %>%
    mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
    arrange(p_value)
  
  # --- Prediction metrics (all 6 models) ------------------------------
  #     Safely handle missing pm_m5 / pm_m6
  has_m5 <- exists("pm_m5") && !is.null(pm_m5)
  has_m6 <- exists("pm_m6") && !is.null(pm_m6)
  
  sheets[["Prediction_Metrics"]] <- tibble(
    model         = c("Model 1: GLM", "Model 2: Spatial-only",
                      "Model 3: Full INLA", "Model 4: Beta-binomial",
                      "Model 5: ZIB-Fixed", "Model 6: ZIB-Estimated"),
    correlation_r = c(pm_glm$r,   pm_sp$r,   pm_full$r,   pm_bb$r,
                      if (has_m5) pm_m5$r   else NA,
                      if (has_m6) pm_m6$r   else NA),
    R_squared     = c(pm_glm$r2,  pm_sp$r2,  pm_full$r2,  pm_bb$r2,
                      if (has_m5) pm_m5$r2  else NA,
                      if (has_m6) pm_m6$r2  else NA),
    MAE           = c(pm_glm$mae, pm_sp$mae, pm_full$mae, pm_bb$mae,
                      if (has_m5) pm_m5$mae else NA,
                      if (has_m6) pm_m6$mae else NA),
    DIC           = c(NA, mc$sp_dic, mc$full_dic, mc$bb_dic,
                      mc$m5_dic, mc$m6_dic),
    WAIC          = c(NA, mc$sp_waic, mc$full_waic, mc$bb_waic,
                      mc$m5_waic, mc$m6_waic),
    AIC           = c(mc$glm_aic, NA, NA, NA, NA, NA),
    AUC_binary    = c(mc$glm_auc, NA, NA, NA, NA, NA)
  )
  
  # --- Model 5 detail sheets -----------------------------------------
  if (!is.null(m5_result)) {
    sheets[["M5_ZIB_Fixed"]] <- tibble(
      parameter = c("pi0", "pi1", "binomial_weight", "spatial_range_km",
                    "spatial_variance", "DIC", "WAIC", "r", "R2", "MAE"),
      value = c(
        round(m5_result$pi0, 6),
        0,
        round(1 - m5_result$pi0, 6),
        round(m5_result$spatial_range, 1),
        round(m5_result$spatial_var, 4),
        round(m5_result$dic, 1),
        round(m5_result$waic, 1),
        if (has_m5) pm_m5$r  else NA,
        if (has_m5) pm_m5$r2 else NA,
        if (has_m5) pm_m5$mae else NA
      )
    )
    if (!is.null(m5_result$fixed_effects))
      sheets[["M5_Fixed_Effects"]] <- m5_result$fixed_effects %>%
        mutate(across(where(is.numeric), ~ round(.x, 4)))
  }
  
  # --- Model 6 detail sheets -----------------------------------------
  if (!is.null(m6_result)) {
    sheets[["M6_ZIB_Estimated"]] <- tibble(
      parameter = c("mean_pi0", "pi1", "spatial_range_km",
                    "spatial_variance", "DIC", "WAIC", "r", "R2", "MAE"),
      value = c(
        round(mean(m6_result$pred_pi0), 6),
        0,
        round(m6_result$spatial_range, 1),
        round(m6_result$spatial_var, 4),
        round(m6_result$dic, 1),
        round(m6_result$waic, 1),
        if (has_m6) pm_m6$r  else NA,
        if (has_m6) pm_m6$r2 else NA,
        if (has_m6) pm_m6$mae else NA
      )
    )
    if (!is.null(m6_result$fixed_effects))
      sheets[["M6_Fixed_Effects"]] <- m6_result$fixed_effects %>%
        mutate(across(where(is.numeric), ~ round(.x, 4)))
  }
  
  # --- Cross-validation sheets ----------------------------------------
  if (exists("cv_summary") && !is.null(cv_summary))
    sheets[["Cross_Validation_Summary"]] <- cv_summary %>%
    mutate(across(where(is.numeric), ~ round(.x, 6)))
  
  if (exists("cv_predictions") && !is.null(cv_predictions))
    sheets[["Cross_Validation_Predictions"]] <- cv_predictions %>%
    mutate(across(where(is.numeric), ~ round(.x, 6)))
  
  # --- Cluster predictions (all models) --------------------------------
  sheets[["Cluster_Predictions"]] <- sub_m %>%
    dplyr::select(LONGNUM, LATNUM, n_children, n_vaccinated, coverage,
                  any_of(c("pred_glm",  "resid_glm",
                           "pred_sp",   "resid_sp",
                           "pred_full", "resid_full",
                           "pred_bb",   "resid_bb",
                           "pred_m5",   "resid_m5",
                           "pred_m6",   "resid_m6",
                           "pred_best", "resid_best"))) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))
  
  # --- Variable audit -------------------------------------------------
  if (exists("var_audit") && !is.null(var_audit))
    sheets[["Variable_Audit"]] <- var_audit %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))
  
  # --- Correlation / VIF drop logs ------------------------------------
  if (exists("cor_drop_log") && nrow(cor_drop_log) > 0)
    sheets[["Correlation_Drops"]] <- cor_drop_log
  if (exists("vif_drop_log") && nrow(vif_drop_log) > 0)
    sheets[["VIF_Drops"]] <- vif_drop_log
  
  # --- Coverage by sample size ----------------------------------------
  if (exists("cov_by_size_country") && !is.null(cov_by_size_country))
    sheets[["Coverage_By_Size"]] <- cov_by_size_country
  
  # --- Model equation -------------------------------------------------
  sheets[["Model_Equation"]] <- tibble(
    item = c("Model type", "Formula (R syntax)",
             "Equation (fitted)", "Intercept",
             paste0("Covariate: ", vars_final)),
    value = c(
      if (exists("eq_label")) eq_label else "see best_model in Model_Summary",
      if (exists("glm_form")) paste(deparse(glm_form), collapse = " ") else "N/A",
      if (exists("eq_full"))  eq_full  else "N/A",
      if (!is.null(feff)) as.character(round(feff$mean[feff$term == "intercept"], 4))
      else if (!is.null(glm_co)) as.character(round(glm_co$Estimate[glm_co$term == "(Intercept)"], 4))
      else "N/A",
      if (!is.null(feff)) {
        sapply(vars_final, function(v) {
          row <- feff %>% dplyr::filter(term == v)
          if (nrow(row) == 1)
            paste0("\u03B2=", round(row$mean, 4), " (OR=", round(row$OR, 4),
                   ", 95% CrI: ", round(row$OR_lo, 4), "-", round(row$OR_hi, 4),
                   ")", ifelse(row$sig == "*", " *", ""))
          else "not in INLA model"
        })
      } else if (!is.null(glm_co)) {
        sapply(vars_final, function(v) {
          row <- glm_co %>% dplyr::filter(term == v)
          if (nrow(row) == 1)
            paste0("\u03B2=", round(row$Estimate, 4), " (OR=", round(row$OR, 4),
                   ", 95% CI: ", round(row$OR_lo, 4), "-", round(row$OR_hi, 4),
                   ") ", row$sig)
          else "not in GLM"
        })
      } else rep("N/A", length(vars_final))
    )
  )
  
  # --- Spatial surface grid -------------------------------------------
  if (!is.null(surface_diag) && "grid" %in% names(surface_diag))
    sheets[["Spatial_Surface_Grid"]] <- surface_diag$grid %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))
  
  # --- Write the workbook ---------------------------------------------
  country_xl_path <- file.path(
    country_dir,
    paste0(country_slug, "_geostatistical_results_unweighted", version, ".xlsx")
  )
  writexl::write_xlsx(sheets, country_xl_path)
  all_xl_paths[[cname]] <- country_xl_path
  cat("  Excel workbook saved to:", country_xl_path, "\n")
  
  
  
  # ═══════════════════════════════════════════════════════════════════════════════
  # MASTER BUNDLE — ONE RDS TO RULE THEM ALL
  # ═══════════════════════════════════════════════════════════════════════════════
  #
  # PURPOSE
  #   Save a SINGLE RDS file per country that contains every piece of data
  #   needed to recreate every page of the PDF, every sheet of the Excel
  #   workbook, and every intermediate result — without re-running any models.
  #
  # PLACEMENT
  #   Drop this block AFTER the Excel workbook is written and BEFORE
  #   the final "Country complete" message at the bottom of the country loop.
  #
  # USAGE LATER
  #   bundle <- readRDS("path/to/country_master_bundle.rds")
  #   # Then: bundle$sub_m, bundle$feff, bundle$surface_grid, etc.
  #
  # ═══════════════════════════════════════════════════════════════════════════════
  
  cat("\n╔══════════════════════════════════════════════════════════╗\n")
  cat("║  SAVING MASTER BUNDLE FOR:", cname, "\n")
  cat("╚══════════════════════════════════════════════════════════╝\n\n")
  
  master_bundle <- list(
    
    # ─────────────────────────────────────────────────────────────────────
    
    # A.  IDENTIFIERS & METADATA
    # ─────────────────────────────────────────────────────────────────────
    country_name  = cname,
    country_slug  = country_slug,
    country_iso   = if (exists("country_iso")) country_iso else NA_character_,
    version       = version,
    output_folder = output_folder,
    country_dir   = country_dir,
    pdf_path      = pdf_path,
    xl_path       = country_xl_path,
    timestamp     = Sys.time(),
    
    # ─────────────────────────────────────────────────────────────────────
    # B.  RAW / FILTERED DATA
    #     sub   = full country subset (before complete-case drop)
    #     sub_m = complete-case modelling dataset (all pred/resid columns)
    #     sub_sf = sf version of sub_m with geometry
    # ─────────────────────────────────────────────────────────────────────
    sub       = if (exists("sub"))    sub    else NULL,
    sub_m     = sub_m,
    sub_sf    = sub_sf,
    
    # ─────────────────────────────────────────────────────────────────────
    # C.  ADMIN BOUNDARIES
    # ─────────────────────────────────────────────────────────────────────
    adm1_sf   = adm1_sf,
    
    # ─────────────────────────────────────────────────────────────────────
    # D.  VARIABLE SELECTION PIPELINE
    # ─────────────────────────────────────────────────────────────────────
    #  D1. All candidate variables and what passed each step
    avail_vars     = if (exists("avail_vars"))   avail_vars   else NULL,
    vars_pass      = if (exists("vars_pass"))    vars_pass    else NULL,
    biv_pass       = if (exists("biv_pass"))     biv_pass     else NULL,
    vars_after_cor = if (exists("vars_after_cor")) vars_after_cor else NULL,
    vars_final     = vars_final,
    
    #  D2. Screening tables
    biv            = biv,
    var_audit      = if (exists("var_audit")) var_audit else NULL,
    cor_drop_log   = if (exists("cor_drop_log")) cor_drop_log else NULL,
    vif_drop_log   = if (exists("vif_drop_log")) vif_drop_log else NULL,
    
    #  D3. Missingness & variance vectors
    miss_pct       = if (exists("miss_pct")) miss_pct else NULL,
    var_sd         = if (exists("var_sd"))   var_sd   else NULL,
    
    #  D4. Coverage by cluster size
    cov_by_size_country = if (exists("cov_by_size_country")) cov_by_size_country else NULL,
    
    # ─────────────────────────────────────────────────────────────────────
    # E.  MODEL 1: GLM
    # ─────────────────────────────────────────────────────────────────────
    fit_glm   = fit_glm,
    glm_co    = glm_co,
    glm_form  = if (exists("glm_form")) glm_form else NULL,
    auc_g     = auc_g,
    
    # ─────────────────────────────────────────────────────────────────────
    # F.  MODELS 2–4: INLA-SPDE
    # ─────────────────────────────────────────────────────────────────────
    #  F1. Model objects
    inla_sp   = inla_sp,
    inla_fu   = inla_fu,
    inla_bb   = inla_bb,
    
    #  F2. Fixed effects from full model
    feff      = if (exists("feff")) feff else NULL,
    
    #  F3. Spatial parameter estimates
    rp        = rp,   # range (full model)
    sp        = sp,   # variance (full model)
    sp_range  = if (exists("sp_range")) sp_range else NA,  # range (spatial-only)
    sp_var    = if (exists("sp_var"))   sp_var   else NA,   # variance (spatial-only)
    bb_range  = if (exists("bb_range")) bb_range else NA,   # range (beta-binomial)
    bb_var    = if (exists("bb_var"))   bb_var   else NA,    # variance (beta-binomial)
    
    #  F4. Mesh & SPDE infrastructure
    mesh      = mesh,
    spde      = spde,
    stk       = if (exists("stk")) stk else NULL,
    A         = if (exists("A"))   A   else NULL,
    
    #  F5. Coordinates & projection
    coords    = if (exists("coords")) coords else NULL,
    cp        = if (exists("cp"))     cp     else NULL,  # projected coordinates (km)
    cent      = if (exists("cent"))   cent   else NULL,  # centroid used for projection
    lcrs      = if (exists("lcrs"))   lcrs   else NULL,  # local CRS string
    
    #  F6. Mesh construction parameters
    rg             = if (exists("rg"))             rg             else NA,
    max_edge       = if (exists("max_edge"))       max_edge       else NULL,
    offset_value   = if (exists("offset_value"))   offset_value   else NULL,
    cutoff_value   = if (exists("cutoff_value"))   cutoff_value   else NA,
    
    #  F7. Prior specification (critical for sensitivity analysis)
    prior_info         = if (exists("prior_info"))         prior_info         else NULL,
    prior_on_range     = if (exists("prior_on_range"))     prior_on_range     else NULL,
    prior_on_variance  = if (exists("prior_on_variance"))  prior_on_variance  else NULL,
    nugget_effect      = if (exists("nugget_effect"))      nugget_effect      else NULL,
    
    #  F8. INLA formulas (as R formula objects AND deparsed text)
    f2             = if (exists("f2")) f2 else NULL,
    f3             = if (exists("f3")) f3 else NULL,
    f2_text        = if (exists("f2")) deparse(f2) else NULL,
    f3_text        = if (exists("f3")) deparse(f3) else NULL,
    iid_spec       = if (exists("iid_spec"))   iid_spec   else NULL,
    cov_terms      = if (exists("cov_terms"))  cov_terms  else NULL,
    
    #  F9. INLA stack indices (for mapping fitted values back to data)
    idx_sp  = if (exists("idx_sp"))  idx_sp  else NULL,
    idx     = if (exists("idx"))     idx     else NULL,  # full model index
    idx_bb  = if (exists("idx_bb"))  idx_bb  else NULL,
    
    # ─────────────────────────────────────────────────────────────────────
    # G.  MODELS 5–6: ZERO-INFLATED BINOMIAL
    # ─────────────────────────────────────────────────────────────────────
    m5_result      = m5_result,
    m6_result      = m6_result,
    boundary_probs = boundary_probs,
    
    # ─────────────────────────────────────────────────────────────────────
    # H.  MODEL COMPARISON & PREDICTION METRICS
    # ─────────────────────────────────────────────────────────────────────
    mc         = mc,
    pred_table = pred_table,
    pm_glm     = pm_glm,
    pm_sp      = pm_sp,
    pm_full    = pm_full,
    pm_bb      = pm_bb,
    pm_m5      = pm_m5,
    pm_m6      = pm_m6,
    
    # ─────────────────────────────────────────────────────────────────────
    # I.  CROSS-VALIDATION
    # ─────────────────────────────────────────────────────────────────────
    cv_result      = if (exists("cv_result"))      cv_result      else NULL,
    cv_summary     = if (exists("cv_summary"))     cv_summary     else NULL,
    cv_predictions = if (exists("cv_predictions")) cv_predictions else NULL,
    
    # ─────────────────────────────────────────────────────────────────────
    # J.  SURFACE DIAGNOSTICS (PREDICTION GRID + PLOTS DATA)
    # ─────────────────────────────────────────────────────────────────────
    #  J1. The full surface_diag object (contains grid, plot objects, etc.)
    surface_diag = surface_diag,
    
    #  J2. Grid extracted separately for convenience
    surface_grid = if (!is.null(surface_diag) && "grid" %in% names(surface_diag))
      surface_diag$grid else NULL,
    
    # ─────────────────────────────────────────────────────────────────────
    # K.  EQUATION STRINGS (for the summary text page)
    # ─────────────────────────────────────────────────────────────────────
    eq_full    = if (exists("eq_full"))  eq_full  else NULL,
    eq_label   = if (exists("eq_label")) eq_label else NULL,
    eq_fixed   = if (exists("eq_fixed")) eq_fixed else NULL,
    bl         = if (exists("bl")) bl else NULL,  # best model label string
    
    # ─────────────────────────────────────────────────────────────────────
    # L.  EXCEL SHEETS (the exact list that was written to xlsx)
    # ─────────────────────────────────────────────────────────────────────
    excel_sheets = sheets,
    
    # ─────────────────────────────────────────────────────────────────────
    # M.  SESSION & REPRODUCIBILITY
    # ─────────────────────────────────────────────────────────────────────
    session_info = sessionInfo(),
    R_version    = R.version.string
  )
  
  
  # ─────────────────────────────────────────────────────────────────────────────
  # SAVE
  # ─────────────────────────────────────────────────────────────────────────────
  
  bundle_path <- file.path(
    country_dir,
    paste0(country_slug, "_master_bundle", version, ".rds")
  )
  
  saveRDS(master_bundle, bundle_path)
  
  bundle_size <- file.size(bundle_path)
  bundle_size_mb <- round(bundle_size / 1024^2, 1)
  
  cat("\n  ✓ Master bundle saved:", bundle_path, "\n")
  cat("    Size:", bundle_size_mb, "MB\n")
  cat("    Contents:", length(master_bundle), "top-level objects\n")
  cat("    Sections: A (metadata), B (data), C (boundaries), D (variable selection),\n")
  cat("              E (GLM), F (INLA spatial), G (ZIB), H (comparison),\n")
  cat("              I (cross-validation), J (surfaces), K (equations),\n")
  cat("              L (excel sheets), M (session info)\n\n")
  
  # ─────────────────────────────────────────────────────────────────────────────
  # ALSO ACCUMULATE INTO CROSS-COUNTRY MASTER LIST
  # ─────────────────────────────────────────────────────────────────────────────
  
  if (!exists("all_master_bundles")) all_master_bundles <- list()
  all_master_bundles[[cname]] <- bundle_path
  
  # -------------------------------------------------------------------
  # 13.  Store lightweight results for cross-country summaries
  # -------------------------------------------------------------------
  
  all_sheets[[cname]] <- sheets
  
  all_results[[cname]] <- list(
    vars_used   = vars_final,
    n_clusters  = nrow(sub_m),
    glm_coefs   = glm_co,
    inla_fixed  = if (exists("feff")) feff else NULL,
    range_km    = rp,
    spatial_var = sp,
    best_model  = if (exists("bl")) bl else mc$best_model
  )
  
  cat("\n Country", cname, "complete.\n")
  cat("    Output directory:", country_dir, "\n")
  
  
  
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










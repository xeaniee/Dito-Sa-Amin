# thesis_master_final.R
# Final master script: RQ1 -> RQ8 with assumption checks and visualizations
# Dataframe: DSA_df

# ---------------------------
# 0) Packages
# ---------------------------
pkgs <- c(
  "tidyverse", "broom", "MASS", "car", "lmtest", "sandwich", "pscl",
  "brant", "psych", "showtext", "scales", "dotwhisker", "ggeffects",
  "forcats", "DescTools"
)
inst <- installed.packages()[, "Package"]
for(p in pkgs) if(!(p %in% inst)) install.packages(p, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---------------------------
# 1) Global theme, font, palette & folders
# ---------------------------
palette_vec <- c("#08326e", "#1a63a8", "#3d8cc1", "#79b3d6", "#aacee5")
scale_fill_cmu <- function(...) scale_fill_manual(values = palette_vec, ...)
scale_color_cmu <- function(...) scale_color_manual(values = palette_vec, ...)

# font (CMU Serif fallback to system serif)
font_family <- "CMU Serif"
showtext_auto()
if(!(font_family %in% font_families())) {
  message("CMU Serif not found — falling back to system 'serif'. To use CMU Serif, call showtext::font_add() with the font file path.")
  font_family <- "serif"
}
theme_set(theme_minimal(base_family = font_family))
ggplot2::theme_update(text = element_text(family = font_family))

dir.create("thesis_outputs", showWarnings = FALSE)
dir.create("thesis_outputs/plots", showWarnings = FALSE)
dir.create("thesis_outputs/tables", showWarnings = FALSE)

# ---------------------------
# 2) Basic dataset & variable checks
# ---------------------------
# Ensure DSA_df exists
if(!exists("DSA_df")) stop("DSA_df not found. Load your dataset into an object named DSA_df and re-run.")

# Check & label community_type
if("community_type" %in% names(DSA_df)) {
  if(is.numeric(DSA_df$community_type) && all(unique(na.omit(DSA_df$community_type)) %in% c(0,1))) {
    DSA_df <- DSA_df %>% mutate(community_type = factor(community_type, levels = c(0,1), labels = c("Non-BHS","BHS")))
  } else {
    DSA_df <- DSA_df %>% mutate(community_type = as.factor(community_type))
  }
} else stop("community_type variable missing in DSA_df.")

# Utility: write safe CSV wrapper
safe_write_csv <- function(df, path) {
  tryCatch(write.csv(df, path, row.names = FALSE), error = function(e) message("Write failed: ", path))
}

# ---------------------------
# Helper functions: assumption checks & reporting
# ---------------------------
# 1) Chi-square expected counts check (returns "chi" or "fisher" and test object)
chi_or_fisher <- function(x, y, data) {
  tb <- table(data[[x]], data[[y]])
  tb <- tb[rowSums(tb) > 0, colSums(tb) > 0, drop = FALSE]
  res <- NULL
  # try chi-square; if expected < 5 in many cells, fallback to fisher
  ch <- tryCatch(chisq.test(tb), warning = function(w) w, error = function(e) e)
  if(inherits(ch, "warning") || inherits(ch, "error")) {
    res <- fisher.test(tb)
    method <- "Fisher"
  } else {
    expc <- ch$expected
    if(any(expc < 5)) {
      res <- fisher.test(tb)
      method <- "Fisher"
    } else {
      res <- ch
      method <- "Chi-square"
    }
  }
  list(method = method, table = tb, test = res)
}

# 2) MLR assumption checks
mlr_assumptions <- function(model) {
  # model: lm object
  res <- list()
  # Linearity: component + residual plots (car::crPlots) — returns invisible plot; we provide VIF and stats
  res$vif <- tryCatch(car::vif(model), error = function(e) NA)
  # Normality of residuals
  resid <- residuals(model)
  res$shapiro_p <- tryCatch(shapiro.test(resid)$p.value, error = function(e) NA)
  # Breusch-Pagan for heteroscedasticity
  res$bptest_p <- tryCatch(lmtest::bptest(model)$p.value, error = function(e) NA)
  # Durbin-Watson for independence of errors (if needed)
  res$dw_p <- tryCatch(lmtest::dwtest(model)$p.value, error = function(e) NA)
  res
}

# 3) OLR proportional odds (Brant)
olr_assumptions <- function(polr_model, data) {
  res <- list()
  # brant test
  br <- tryCatch(brant::brant(polr_model), error = function(e) NULL)
  res$brant <- br
  # VIF: need to fit a linear model on same predictors for VIF
  terms <- names(coef(polr_model))
  # compute VIF on underlying numeric data if possible
  res$note <- "Check VIF on numeric predictors with car::vif on an equivalent lm model if needed."
  res
}

# 4) BLR: linearity of logit for continuous predictors (Box-Tidwell)
blr_assumptions <- function(model, data) {
  res <- list()
  # VIF
  res$vif <- tryCatch(car::vif(model), error = function(e) NA)
  res
}

# 5) Separation check for logistic models
check_separation <- function(glm_model) {
  # detect separation via check of large coefficients / warning; simple heuristic:
  s <- summary(glm_model)
  if(any(is.infinite(coef(glm_model)))) return(TRUE)
  # Check if any predictor perfectly predicts outcome: table cross-tabs
  sep_vars <- c()
  for(term in names(coef(glm_model))[-1]) {
    if(term %in% names(glm_model$model)) {
      tab <- table(glm_model$model[[term]], glm_model$model[[1]])
      if(any(rowSums(tab)==0) || any(colSums(tab)==0)) sep_vars <- c(sep_vars, term)
    }
  }
  list(separation = length(sep_vars) > 0, vars = sep_vars)
}


# ---------------------------
# ---- RQ1: Descriptive Statistics ----
# Means, SDs for continuous; Freqs for categorical, by community_type
# ---------------------------
cat("\n---- RQ1: Descriptive Statistics ----\n")

# Identify variable types (adjust based on your actual data structure if needed)
cont_vars_guess <- names(DSA_df)[sapply(DSA_df, function(x) is.numeric(x) && length(unique(na.omit(x))) > 10)]
cat_vars_guess <- setdiff(names(DSA_df), c(cont_vars_guess, "community_type"))
cat_vars_guess <- cat_vars_guess[sapply(DSA_df[cat_vars_guess], function(x) !is.data.frame(x) && !is.list(x))]

# --- Continuous Variable Descriptives ---
desc_cont <- DSA_df %>%
  group_by(community_type) %>%
  summarise(across(all_of(cont_vars_guess),
                   list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE)),
                   .names = "{.col}_{.fn}")) %>%
  pivot_longer(cols = -community_type, names_to = "variable", values_to = "value") %>%
  separate(variable, into = c("variable", "stat"), sep = "_(mean|sd)$") %>%
  pivot_wider(names_from = stat, values_from = value)
safe_write_csv(desc_cont, "thesis_outputs/tables/RQ1_Descriptives_Continuous.csv")

# --- Categorical Variable Descriptives (Including first_point_of_contact) ---
if("first_point_of_contact" %in% names(DSA_df)) {
  # Ensure first_point_of_contact is treated as a factor for correct display
  DSA_df$first_point_of_contact_fct <- factor(DSA_df$first_point_of_contact,
                                              levels = 0:5,
                                              labels = c("Hindi nagpapatingin", "BHS", "RHU", "TAM", "Distant public facility", "Distant private facility"))
  
  # Add the new variable to the list for frequency count
  cat_vars_for_freq <- unique(c("first_point_of_contact_fct", cat_vars_guess))
  
  desc_cat <- list()
  for (v in cat_vars_for_freq) {
    if (v %in% names(DSA_df)) {
      ftab <- DSA_df %>%
        filter(!is.na(.data[[v]])) %>%
        group_by(community_type, .data[[v]]) %>%
        summarise(n = n(), .groups = "drop_last") %>%
        mutate(percent = n / sum(n) * 100) %>%
        ungroup() %>%
        rename(level = v) %>%
        mutate(variable = v)
      desc_cat[[v]] <- ftab
    }
  }
  desc_cat_df <- bind_rows(desc_cat)
  safe_write_csv(desc_cat_df, "thesis_outputs/tables/RQ1_Descriptives_Categorical.csv")
} else {
  message("RQ1: 'first_point_of_contact' variable missing. Categorical descriptives may be incomplete.")
}


# ---------------------------
# ---- RQ2: Chi-square tests (Aggregate) ----
# personal_illness_consult and fam_mem_consult ~ community_type
# ---------------------------
# RQ2a: personal illness consultation
if("personal_illness_consult" %in% names(DSA_df)) {
  cat("\n---- RQ2a: personal_illness_consult x community_type ----\n")
  r_q2a <- chi_or_fisher("personal_illness_consult", "community_type", DSA_df)
  results_rq2a <- list(method = r_q2a$method, table = r_q2a$table, test = r_q2a$test)
  safe_write_csv(as.data.frame(r_q2a$table), "thesis_outputs/tables/RQ2_personal_table.csv")
  # Plot stacked proportions
  p_rq2a <- DSA_df %>%
    filter(!is.na(personal_illness_consult)) %>%
    group_by(community_type, personal_illness_consult) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(community_type) %>% mutate(prop = n/sum(n)) %>%
    ggplot(aes(x = community_type, y = prop, fill = as.factor(personal_illness_consult))) +
    geom_col(position = "fill") + scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_cmu(name = "Response") + labs(title = "Personal illness consultation by community type", x = "", y = "Percent")
  ggsave("thesis_outputs/plots/RQ2_personal_stack.png", p_rq2a, width = 7, height = 4, dpi = 300)
}

# RQ2b: family illness consultation
if("fam_mem_consult" %in% names(DSA_df)) {
  cat("\n---- RQ2b: fam_mem_consult x community_type ----\n")
  r_q2b <- chi_or_fisher("fam_mem_consult", "community_type", DSA_df)
  results_rq2b <- list(method = r_q2b$method, table = r_q2b$table, test = r_q2b$test)
  safe_write_csv(as.data.frame(r_q2b$table), "thesis_outputs/tables/RQ2_family_table.csv")
  p_rq2b <- DSA_df %>%
    filter(!is.na(fam_mem_consult)) %>%
    group_by(community_type, fam_mem_consult) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(community_type) %>% mutate(prop = n/sum(n)) %>%
    ggplot(aes(x = community_type, y = prop, fill = as.factor(fam_mem_consult))) +
    geom_col(position = "fill") + scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_cmu(name = "Response") + labs(title = "Family illness consultation by community type", x = "", y = "Percent")
  ggsave("thesis_outputs/plots/RQ2_family_stack.png", p_rq2b, width = 7, height = 4, dpi = 300)
}

# ---------------------------
# ---- RQ3: Multiple Linear Regression ----
# ---------------------------
base_ivars <- c("distance", "waiting_time", "awareness", "affordability", "available_services", "neighbordhood_perception")
facilities <- c("BHS", "RHU", "TAM")

for(ct in levels(DSA_df$community_type)) {
  cat("\n==== RQ3 models for community_type:", ct, "====\n")
  subset_df <- DSA_df %>% filter(community_type == ct)
  for(fac in facilities) {
    fac_ivars <- paste0(fac, "_", base_ivars)
    dv <- paste0(fac, "_utilization_rate")
    missing_vars <- setdiff(c(dv, fac_ivars), names(subset_df))
    if(length(missing_vars) > 0) {
      message("RQ3 skipped for ", ct, "-", fac, " (missing: ", paste(missing_vars, collapse = ", "), ")")
      next
    }
    formula_str <- paste0(dv, " ~ ", paste(fac_ivars, collapse = " + "))
    model <- lm(as.formula(formula_str), data = subset_df)
    # save model and tidy coefs with robust SE
    tidy_coefs <- broom::tidy(coeftest(model, vcov = sandwich))
    safe_write_csv(tidy_coefs, paste0("thesis_outputs/tables/RQ3_MLR_", ct, "_", fac, "_coefs.csv"))
    saveRDS(model, paste0("thesis_outputs/RQ3_MLR_", ct, "_", fac, "_model.rds"))
    # assumptions
    asum <- mlr_assumptions(model)
    safe_write_csv(as.data.frame(t(asum)), paste0("thesis_outputs/tables/RQ3_MLR_", ct, "_", fac, "_assumptions.csv"))
    # plots: dotwhisker (coeffs) + partial scatter for each predictor
    p_dw <- dotwhisker::dwplot(model) + labs(title = paste0("MLR: ", fac, " utilization — ", ct)) + scale_color_cmu()
    ggsave(paste0("thesis_outputs/plots/RQ3_dw_", ct, "_", fac, ".png"), p_dw, width = 7, height = 4, dpi = 300)
    for(pred in fac_ivars) {
      p_sc <- ggplot(subset_df, aes_string(x = pred, y = dv)) +
        geom_point(alpha = 0.4) + geom_smooth(method = "lm", se = TRUE) +
        labs(title = paste(fac, "utilization vs", pred, "(", ct, ")"), x = pred, y = dv)
      ggsave(paste0("thesis_outputs/plots/RQ3_scatter_", ct, "_", fac, "_", pred, ".png"), p_sc, width = 6, height = 4, dpi = 300)
    }
  }
}

# ---------------------------
# ---- RQ4: Ordinal Logistic Regression ----
# ---------------------------
for(ct in levels(DSA_df$community_type)) {
  cat("\n==== RQ4 models for community_type:", ct, "====\n")
  subset_df <- DSA_df %>% filter(community_type == ct)
  for(fac in facilities) {
    dv <- paste0(fac, "_satisfaction_level")
    fac_ivars <- paste0(fac, "_", base_ivars)
    missing_vars <- setdiff(c(dv, fac_ivars), names(subset_df))
    if(length(missing_vars) > 0) {
      message("RQ4 skipped for ", ct, "-", fac, " (missing: ", paste(missing_vars, collapse = ", "), ")")
      next
    }
    # ensure ordered factor
    subset_df[[dv]] <- ordered(subset_df[[dv]])
    formula_str <- paste0(dv, " ~ ", paste(fac_ivars, collapse = " + "))
    polr_model <- MASS::polr(as.formula(formula_str), data = subset_df, Hess = TRUE, na.action = na.omit)
    saveRDS(polr_model, paste0("thesis_outputs/RQ4_OLR_", ct, "_", fac, "_model.rds"))
    ctab <- coef(summary(polr_model))
    pvals <- pnorm(abs(ctab[, "t value"]), lower.tail = FALSE) * 2
    ctab <- cbind(ctab, "p.value" = pvals)
    write.csv(as.data.frame(ctab), paste0("thesis_outputs/tables/RQ4_OLR_", ct, "_", fac, "_coefs.csv"))
    # proportional odds test (Brant)
    br <- tryCatch(brant::brant(polr_model), error = function(e) { message("Brant test error for ", ct, "-", fac); NULL })
    if(!is.null(br)) saveRDS(br, paste0("thesis_outputs/tables/RQ4_OLR_", ct, "_", fac, "_brant.rds"))
    # Predicted probabilities example using the facility-specific distance if available
    pred_var <- paste0(fac, "_distance")
    if(pred_var %in% names(subset_df)) {
      preds <- tryCatch(ggeffects::ggpredict(polr_model, terms = paste0(pred_var, " [all]")), error = function(e) NULL)
      if(!is.null(preds)) {
        pprob <- ggplot(preds, aes(x = x, y = predicted, color = group)) +
          geom_line() + geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.12, color = NA) +
          labs(title = paste0("Predicted probabilities of satisfaction by ", pred_var, " — ", fac, " (", ct, ")"), x = pred_var, y = "Predicted probability") +
          scale_color_cmu() + scale_fill_cmu()
        ggsave(paste0("thesis_outputs/plots/RQ4_predprob_", ct, "_", fac, ".png"), pprob, width = 7, height = 4, dpi = 300)
      }
    }
  }
}

# ---------------------------
# ---- RQ5: Frequency counts - Coping Mechanism ----
# ---------------------------
cat("\n---- RQ5: Frequency counts - Coping Mechanisms ----\n")
# New list of coping mechanism variables (binary 0/1)
coping_vars <- c("social_support", "religiosity", "forbearance", "emotional_release",
                 "overactivity", "relaxation_recreation", "substance_use",
                 "cognitive_reappraisal", "problem_solving")
missing_coping <- setdiff(coping_vars, names(DSA_df))

if (length(missing_coping) > 0) {
  message("RQ5 skipped: Missing coping variables: ", paste(missing_coping, collapse = ", "))
} else {
  # Calculate frequency (sum of 'yes' = 1) and percentage across the entire sample (or subset if needed)
  coping_df <- DSA_df %>%
    summarise(across(all_of(coping_vars), ~sum(. == 1, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "mechanism", values_to = "count") %>%
    # Use total non-missing sample size for percentage calculation
    mutate(total_n = sum(!is.na(DSA_df[[coping_vars[1]]])), # Use one var as proxy for total n
           percent = count / total_n * 100) %>%
    arrange(desc(count))
  
  safe_write_csv(coping_df, "thesis_outputs/tables/RQ5_coping_mechanisms.csv")
  
  p_rq5 <- ggplot(coping_df, aes(x = fct_reorder(mechanism, count), y = percent, fill = mechanism)) +
    geom_col(show.legend = FALSE) + coord_flip() + scale_fill_cmu() +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = "Frequency of Coping Mechanism Use", x = "", y = "Percent of Respondents")
  ggsave("thesis_outputs/plots/RQ5_coping_mechanisms.png", p_rq5, width = 7, height = 5, dpi = 300)
}

# ---------------------------
# ---- RQ6: Binary Logistic Regression (BLR) ----
# NOTE: RQ6 was designed to use the *justification* variables (RQ5 old list) as IVs.
# Assuming you want to KEEP the *original justification* variables as IVs for RQ6,
# as they relate to limited_healthseeking_capacity. If you intended to use the *coping*
# variables as IVs for RQ6, you must let me know.
# ---------------------------
just_vars_for_rq6 <- c("cannot_afford_treatment", "hosp_is_too_far", "afraid_to_know_illness", "afraid_of_huge_bill",
                       "distrust_nearby_doctors", "distrust_clinical_medicine", "symptom_not_too_serious",
                       "no_free_time", "long_queues")

if("limited_healthseeking_capacity" %in% names(DSA_df)) {
  missing_iv <- setdiff(just_vars_for_rq6, names(DSA_df))
  if(length(missing_iv) > 0) {
    message("RQ6 BLR: missing IVs: ", paste(missing_iv, collapse = ", "))
  }
  # Build formula using the justification vars that exist
  ivs_present <- intersect(just_vars_for_rq6, names(DSA_df))
  if (length(ivs_present) > 0) {
    blr_formula <- as.formula(paste("limited_healthseeking_capacity ~", paste(ivs_present, collapse = " + ")))
    # Ensure DV is a factor if it's 0/1 numeric for correct binomial modeling
    DSA_df$limited_healthseeking_capacity <- as.factor(DSA_df$limited_healthseeking_capacity)
    blr_model <- glm(blr_formula, data = DSA_df, family = binomial)
    saveRDS(blr_model, "thesis_outputs/RQ6_BLR_model.rds")
    tidy_blr <- broom::tidy(coeftest(blr_model, vcov = sandwich))
    safe_write_csv(tidy_blr, "thesis_outputs/tables/RQ6_BLR_coefs.csv")
    # Odds ratios
    ORs <- tryCatch(exp(cbind(Estimate = coef(blr_model), confint(blr_model))), error = function(e) exp(coef(blr_model)))
    safe_write_csv(as.data.frame(ORs), "thesis_outputs/tables/RQ6_BLR_ORs.csv")
    # Assumptions: VIF
    asum_blr <- blr_assumptions(blr_model, DSA_df)
    safe_write_csv(as.data.frame(t(asum_blr)), "thesis_outputs/tables/RQ6_BLR_assumptions.csv")
    # Linearity of logit: Box-Tidwell on continuous predictors (if any)
    cont_predictors <- ivs_present[sapply(DSA_df[ivs_present], is.numeric)]
    if(length(cont_predictors) > 0) {
      bt_res <- list()
      for(pred in cont_predictors) {
        safe <- tryCatch({
          car::boxTidwell(as.formula(paste("limited_healthseeking_capacity ~", pred)), data = DSA_df)
        }, error = function(e) NULL)
        bt_res[[pred]] <- safe
      }
      saveRDS(bt_res, "thesis_outputs/tables/RQ6_BoxTidwell_results.rds")
    }
    # Separation check
    sep <- check_separation(blr_model)
    saveRDS(sep, "thesis_outputs/tables/RQ6_separation_check.rds")
    # Coefficient plot
    p_blr_coef <- ggplot(tidy_blr %>% filter(term != "(Intercept)"), aes(x = reorder(term, estimate), y = estimate)) +
      geom_point() + geom_errorbar(aes(ymin = estimate - 1.96*std.error, ymax = estimate + 1.96*std.error), width = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      coord_flip() + labs(title = "BLR coefficients (log-odds) predicting limited_healthseeking_capacity", x = "", y = "Log-odds")
    ggsave("thesis_outputs/plots/RQ6_blr_coefs.png", p_blr_coef, width = 7, height = 5, dpi = 300)
  } else {
    message("RQ6 BLR skipped: No valid IVs present.")
  }
} else {
  message("RQ6 skipped: 'limited_healthseeking_capacity' variable missing in DSA_df.")
}

# ---------------------------
# ---- RQ7: Frequency counts (symptom severity & observation thresholds) ----
# ---------------------------
cat("\n---- RQ7: Symptom severity and observation thresholds ----\n")
rq7_vars <- c("symptom_severity_to_consult", "symptom_observed_to_consult")
for(v in rq7_vars) {
  if(v %in% names(DSA_df)) {
    ftab <- DSA_df %>% count(.data[[v]]) %>% mutate(percent = n / sum(n, na.rm=TRUE) * 100)
    safe_write_csv(ftab, paste0("thesis_outputs/tables/RQ7_freq_", v, ".csv"))
    
    # Check if we can label the severity variable
    v_title <- v
    if (v == "symptom_severity_to_consult") {
      DSA_df$severity_fct <- factor(DSA_df$symptom_severity_to_consult,
                                    levels = 1:5,
                                    labels = c("Not severe", "Slightly not severe", "Neutral", "Slightly severe", "Severe"))
      p7 <- ggplot(DSA_df, aes(x = severity_fct, y = after_stat(count) / sum(after_stat(count)) * 100, fill = severity_fct)) +
        geom_bar(show.legend = FALSE) + scale_fill_cmu() +
        labs(title = paste("Distribution:", v), x = v, y = "Percent")
    } else {
      p7 <- ggplot(ftab, aes(x = as.factor(.data[[v]]), y = percent, fill = as.factor(.data[[v]]))) +
        geom_col(show.legend = FALSE) + scale_fill_cmu() + labs(title = paste("Distribution:", v), x = v, y = "Percent")
    }
    
    ggsave(paste0("thesis_outputs/plots/RQ7_freq_", v, ".png"), p7, width = 7, height = 4, dpi = 300)
  } else {
    message("RQ7 skipped for: '", v, "' (variable missing in DSA_df).")
  }
}

# ---------------------------
# ---- RQ8: Othering Index creation, reliability, MLR, assumptions, viz ----
# ---------------------------
other_items <- c("limited_healthseeking_capacity", "treatment_better_in_city",
                 "consultation_better_in_city", "loc_gov_assistance")
missing_other <- setdiff(other_items, names(DSA_df))
if(length(missing_other) > 0) message("RQ8: missing items: ", paste(missing_other, collapse = ", "))

# create z-scored items and index
for(it in other_items) {
  if(it %in% names(DSA_df) && is.numeric(DSA_df[[it]])) {
    DSA_df[[paste0("z_", it)]] <- as.numeric(scale(DSA_df[[it]]))
  } else if (it %in% names(DSA_df)) {
    message("RQ8: Item '", it, "' is not numeric and cannot be z-scored for the index.")
  }
}
z_items <- intersect(paste0("z_", other_items), names(DSA_df))
if(length(z_items) > 0) DSA_df$othering_index <- rowMeans(DSA_df[, z_items, drop=FALSE], na.rm = TRUE) else stop("RQ8: no numeric items present to compute index.")

# reliability
if(all(other_items %in% names(DSA_df)) && all(sapply(DSA_df[, other_items], is.numeric))) {
  alpha_res <- psych::alpha(DSA_df[, other_items], check.keys = TRUE)
  write.csv(as.data.frame(alpha_res$total), "thesis_outputs/tables/RQ8_alpha_total.csv")
  message("RQ8 Cronbach's alpha (raw): ", round(alpha_res$total$raw_alpha, 3))
}

# distribution plot
p_othering <- ggplot(DSA_df, aes(x = othering_index, fill = community_type)) +
  geom_density(alpha = 0.5) + scale_fill_cmu() +
  labs(title = "Distribution of Othering Index by Community Type", x = "Othering index", y = "Density")
ggsave("thesis_outputs/plots/RQ8_othering_dist.png", p_othering, width = 7, height = 4, dpi = 300)

# MLR: predictors are base_ivars facility-agnostic? use community-level predictors (distance, waiting_time, affordability, awareness, available_services, neighbordhood_perception) if present (non-facility-prefixed)
mlr_predictors <- c("distance", "waiting_time", "affordability", "awareness", "available_services", "neighbordhood_perception")
present_preds <- intersect(mlr_predictors, names(DSA_df))
if(length(present_preds) == 0) message("RQ8: no generic predictors found; check variable naming.")
rq8_formula <- as.formula(paste("othering_index ~", paste(present_preds, collapse = " + "), "+ community_type"))
rq8_model <- lm(rq8_formula, data = DSA_df)
saveRDS(rq8_model, "thesis_outputs/RQ8_MLR_model.rds")
tidy_rq8 <- broom::tidy(coeftest(rq8_model, vcov = sandwich))
safe_write_csv(tidy_rq8, "thesis_outputs/tables/RQ8_MLR_coefs.csv")
# assumptions
rq8_assum <- mlr_assumptions(rq8_model)
safe_write_csv(as.data.frame(t(rq8_assum)), "thesis_outputs/tables/RQ8_assumptions.csv")
# dot-whisker plot
p_rq8_dw <- dotwhisker::dwplot(rq8_model) + labs(title = "RQ8: Predictors of Othering Index") + scale_color_cmu()
ggsave("thesis_outputs/plots/RQ8_dw.png", p_rq8_dw, width = 7, height = 4, dpi = 300)

# ---------------------------
# Save master results objects
# ---------------------------
saveRDS(list(rq2a = results_rq2a, rq2b = results_rq2b), "thesis_outputs/rq2_results.rds")
# Save other models & tables list by listing directory contents (they were saved as files)
cat("\nFinished running RQ1-RQ8. All tables and plots are in thesis_outputs/tables and thesis_outputs/plots\n")
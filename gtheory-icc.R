


# theming for apa
theme_apa2 <- function(.data, pad = 0, spacing = 2) {
  apa.border <- list("width" = 1, color = "black", style = "solid")
  font(.data, part = "all", fontname = "Times New Roman") |>
    line_spacing(space = spacing, part = "all") |>
    padding(padding=pad) |>
    hline_top(part = "head", border = apa.border) |>
    hline_bottom(part = "head", border = apa.border) |>
    hline_top(part = "body", border = apa.border) |>
    hline_bottom(part = "body", border = apa.border) |>
    align(align = "center", part = "all", j=-1) |>
    valign(valign = "center", part = "all") |>
    colformat_double(digits = 2) |>
    # Add footer formatting
    italic(part = "footer") |>
    hline_bottom(part = "footer", border = fp_border(width = 0, color = "transparent")) |>
    fix_border_issues()
}

# Set global defaults for all flextables
apa_caption <- function(text) {
  as_paragraph(as_chunk(text, props = fp_text_default(font.family = "Times New Roman")))
}

# Calculates dependability and generalizability coefficients for a given model speification, 
# returns bootstrapped 95% CI for the coefficient estimates
# returns variance components for the model
# cleans up data at the end
# includes model results for overallmodel
# uses insight::get_variance to get the variance components
# note: for this function (get_variance), the results are validated against the solutions provided by Nakagawa et al. (2017), in particular examples shown in the Supplement 2 of the paper.
# Binomial: For other binomial models, the distribution-specific variance for Bernoulli models is used, divided by a weighting factor based on the number of trials and successes.
calc_icc_formula_binomial <- function(formula, icc_facet, data, R = 1000, decimal = 2) {
  
  # Extract random effect facets from formula
  random_terms <- findbars(formula)
  facets <- map_chr(random_terms, ~deparse(.x[[3]]))
  
  # Helper function to extract variance components consistently
  extract_variance_components <- function(model, verbose = TRUE) {
    var_decomp <- insight::get_variance(model, verbose = verbose, tolerance = 1e-12, approximation = "observation_level")
    
    # Extract components
    var_random <- var_decomp$var.intercept
    var_residual <- var_decomp$var.residual
    var_total <- sum(var_random) + var_residual
    
    list(
      var_random = var_random,
      var_residual = var_residual,
      var_total = var_total,
      full_decomp = var_decomp
    )
  }
  
  # Fit initial model
  warnings_initial <- character(0)
  model <- suppressWarnings(withCallingHandlers({
    glmer(formula, data = data, family = binomial("logit"),
          control = glmerControl(check.conv.singular = "ignore"))
  }, warning = function(w) {
    warnings_initial <<- c(warnings_initial, w$message)
  }))
  
  if(isSingular(model)) {warnings_initial <- c(warnings_initial, "Singular fit detected")}
  
  # Extract variance components
  vc <- extract_variance_components(model)
  
  # Calculate ICC and percentages
  icc <- vc$var_random[icc_facet] / vc$var_total
  var_components <- c(vc$var_random, residual = vc$var_residual)
  var_percentages <- var_components / vc$var_total
  
  # Bootstrap
  message("Start Bootstrap")
  unique_levels <- unique(data[[icc_facet]])
  n_levels <- length(unique_levels)
  
  # Store bootstrap results for ICC and all variance components
  boot_results <- matrix(NA, nrow = R, ncol = length(var_components) + 1)
  colnames(boot_results) <- c(names(var_components), "icc")
  warnings_bootstrap <- character(0)
  
  for(i in 1:R) {
    if(i %% 100 == 0) message(paste0("  Iteration ", i, "/", R))
    
    # Resample participants with replacement
    boot_sample <- sample(1:n_levels, size = n_levels, replace = TRUE)
    sampled_participants <- unique_levels[boot_sample]
    
    # Build bootstrap dataset
    d <- map_dfr(seq_along(sampled_participants), function(idx) {
      participant_id <- sampled_participants[idx]
      participant_data <- data |> filter(!!sym(icc_facet) == participant_id)
      
      # Create unique bootstrap ID
      participant_data[[icc_facet]] <- paste0(participant_id, "_boot", idx)
      
      # Update interaction terms if present
      if("stimuli:participant" %in% facets) {
        interaction_col <- grep(":", names(participant_data), value = TRUE)
        if(length(interaction_col) > 0) {
          participant_data[[interaction_col]] <- paste0(
            participant_data$stimuli, ":", participant_data[[icc_facet]]
          )
        }
      }
      
      participant_data
    })
    
    # Fit bootstrap model and extract components
    tryCatch({
      m <- suppressWarnings(
        glmer(formula, data = d, family = binomial("logit"),
              control = glmerControl(check.conv.singular = "ignore",
                                     optimizer = "bobyqa",
                                     optCtrl = list(maxfun = 2e6)))
      )
      
      vc_boot <- extract_variance_components(m, verbose = FALSE)
      
      # Store variance components
      boot_results[i, names(vc_boot$var_random)] <- vc_boot$var_random
      boot_results[i, "residual"] <- vc_boot$var_residual
      boot_results[i, "icc"] <- vc_boot$var_random[icc_facet] / vc_boot$var_total
      
    }, error = function(e) {
      warnings_bootstrap <<- c(warnings_bootstrap, e$message)
    })
  }
  
  # Calculate CIs for ICC
  valid_boots <- sum(!is.na(boot_results[, "icc"]))
  boot_success_rate <- valid_boots / R
  boot_ci <- quantile(boot_results[, "icc"], probs = c(0.025, 0.975), na.rm = TRUE)
  
  # Create warnings table
  all_warnings <- c(warnings_initial, warnings_bootstrap)
  unique_warnings <- unique(all_warnings)
  
  warnings_table <- if(length(unique_warnings) > 0) {
    warning_counts <- table(all_warnings)
    tibble(
      warning = names(warning_counts),
      count = as.integer(warning_counts)
    )
  } else {
    tibble(warning = character(0), count = integer(0))
  }
  
  # Create formula string
  formula_str <- paste0("σ²", icc_facet, " / (", 
                        paste0("σ²", c(names(vc$var_random), "residual"), collapse = " + "), ")")
  
  # Format results
  results <- list(
    icc_results = tibble(
      icc = round(icc, decimal),
      ci_lower = round(boot_ci[1], decimal),
      ci_upper = round(boot_ci[2], decimal),
      ci_95 = paste0("[", round(boot_ci[1], decimal), ", ", 
                     round(boot_ci[2], decimal), "]"),
      formula = formula_str,
      bootstrap_success_rate = round(boot_success_rate, 2)
    ),
    
    variance_components = tibble(
      component = names(var_components),
      variance = round(var_components, decimal),
      percentage = scales::label_percent(accuracy = (10^-decimal))(var_percentages)
    ),
    
    variance_decomposition = vc$full_decomp,
    warnings = warnings_table,
    parameters = parameters::model_parameters(model, verbose = FALSE ),
    model = model,
    bootstrap_raw = boot_results
  )
  
  return(results)
}

# changes the number of averaged stimuli given teh models that have already been run above
# so it just uses the icc_results object retuned by the above function
# caluclates the coefficient and makes the right adjustments
# not ethis excludes participant:stimuli interaction, per the paper
calc_decision_study_from_bootstrap <- function(icc_results, n_stimuli = 1:5, decimal = 2) {
  
  # Extract stored bootstrap results
  boot_matrix <- icc_results$bootstrap_raw
  
  # Identify relevant columns (single random occasion, n stimuli from fixed set)
  participant_col <- "participant"
  time_col <- "time"
  time_p_col <- "time:participant"
  stimuli_col <- "stimuli"
  time_stim_col <- "time:stimuli"
  residual_col <- "residual"
  # stimuli:participant excluded (fixed stimulus set)
  
  # Calculate Φ for each bootstrap iteration and each n
  phi_matrix <- matrix(NA, nrow = nrow(boot_matrix), ncol = length(n_stimuli))
  
  for(i in 1:nrow(boot_matrix)) {
    if(!is.na(boot_matrix[i, "icc"])) {
      sigma2_p <- boot_matrix[i, participant_col]
      sigma2_time <- boot_matrix[i, time_col]
      sigma2_time_p <- boot_matrix[i, time_p_col]
      sigma2_s <- boot_matrix[i, stimuli_col]
      sigma2_time_s <- boot_matrix[i, time_stim_col]
      sigma2_e <- boot_matrix[i, residual_col]
      
      for(j in seq_along(n_stimuli)) {
        n <- n_stimuli[j]
        phi_matrix[i, j] <- sigma2_p / 
          (sigma2_p + sigma2_time + sigma2_time_p + sigma2_s/n + sigma2_time_s/n + sigma2_e/n)
      }
    }
  }
  
  # Get original variance components
  vc <- icc_results$variance_components
  sigma2_p <- vc$variance[vc$component == participant_col]
  sigma2_time <- vc$variance[vc$component == time_col]
  sigma2_time_p <- vc$variance[vc$component == time_p_col]
  sigma2_s <- vc$variance[vc$component == stimuli_col]
  sigma2_time_s <- vc$variance[vc$component == time_stim_col]
  sigma2_e <- vc$variance[vc$component == residual_col]
  
  original_phi <- map_dbl(n_stimuli, function(n) {
    sigma2_p / (sigma2_p + sigma2_time + sigma2_time_p + sigma2_s/n + sigma2_time_s/n + sigma2_e/n)
  })
  
  # Calculate CIs from bootstrap distribution
  results <- map_dfr(seq_along(n_stimuli), function(i) {
    boot_phi <- phi_matrix[, i]
    boot_phi <- boot_phi[!is.na(boot_phi)]
    
      ci <- quantile(boot_phi, probs = c(0.025, 0.975), na.rm = TRUE)
      tibble(
        n_stimuli = n_stimuli[i],
        phi_coefficient = round(original_phi[i], decimal),
        ci_lower = round(ci[1], decimal),
        ci_upper = round(ci[2], decimal),
        ci_95 = paste0("[", round(ci[1], decimal), ", ", round(ci[2], decimal), "]")
      )

  })
  
  return(results)
}


# calculates the mdc from the binomial model
# mdc is on the logit scale
# option to eclude from error certain terms depending on the g-university
# not used though
calc_mdc_formula_binomial <- function(formula, icc_facet, data, 
                                      exclude_from_error = NULL,
                                      decimal = 2) {
  
  # Fit binomial model
  model <- suppressWarnings(
    glmer(formula, data = data, family = binomial("logit"),
          control = glmerControl(check.conv.singular = "ignore"))
  )
  
  # Extract variance components
  vc <- insight::get_variance(model, tolerance = 1e-12, verbose = FALSE)
  
  # Calculate error variance (exclude specified components)
  error_facets <- setdiff(names(vc$var.intercept), c(icc_facet, exclude_from_error))
  error_variance <- vc$var.residual + sum(vc$var.intercept[error_facets])
  
  # Calculate MDC on logit scale
  sem_logit <- sqrt(error_variance)
  mdc_logit <- 1.96 * sqrt(2) * sem_logit
  
  return(mdc_logit)
}

# calculates the mdc from the binomial model
# mdc is on the logit scale
# option to eclude from error certain terms depending on the g-university
# not used though
calc_sem <- function(formula, icc_facet, data, 
                                      exclude_from_error = NULL,
                                      decimal = 2) {
  
  # Fit binomial model
  model <- suppressWarnings(
    glmer(formula, data = data, family = binomial("logit"),
          control = glmerControl(check.conv.singular = "ignore"))
  )
  
  # Extract variance components
  vc <- insight::get_variance(model, tolerance = 1e-12, verbose = FALSE)
  
  # Calculate error variance (exclude specified components)
  error_facets <- setdiff(names(vc$var.intercept), c(icc_facet, exclude_from_error))
  error_variance <- vc$var.residual + sum(vc$var.intercept[error_facets])
  
  # Calculate MDC on logit scale
  sem_logit <- sqrt(error_variance)

  return(sem_logit)
}

# uses logit MDC to derive the MDC on the probabiltiy scale
# calcualtes this for strata of size ~ 5 but makes the ones on the lower end larger
# if the max score is not divisible by 5, because there are few participants at this
# end of the scale. 
# for each strata, uses the score furthest from the overall scale midpoint to 
# calcualte the MDC for that strata, ensuring that estimatse are conservative, espcially
# at the tails of the distribution
# we were reluctant to calcualte the MDC for each individual 1-increment score because
# these estimates were too small at the tails. 
calc_mdc_by_baseline <- function(mdc_df) {
  
  mdc_df |>
    mutate(
      baseline_data = pmap(list(mdc_logit, max_score), function(mdc_log, max_sc) {
        
        # Calculate 5-point bins with remainder at bottom
        total_values <- max_sc + 1  # 0 to max_sc inclusive
        n_complete_bins <- floor(total_values / 5)
        remainder <- total_values %% 5
        
        # Create breaks
        if (remainder == 0) {
          # Divisible - all bins size 5
          breaks <- seq(0, max_sc, by = 5)
        } else {
          # Put remainder in first bin
          first_bin_size <- 5 + remainder
          breaks <- c(0, first_bin_size, seq(first_bin_size + 5, max_sc, by = 5))
        }
        
        n_strata <- length(breaks) - 1
        scale_midpoint <- max_sc / 2
        
        tibble(stratum = 1:n_strata) |>
          mutate(
            # Define stratum ranges
            min_score = breaks[stratum],
            max_score_range = breaks[stratum + 1],
            
            # Get all integer scores in this range
            scores_in_range = map2(min_score, max_score_range, function(min_s, max_s) {
              ceiling(min_s):floor(max_s - 1)  # -1 because max_s is exclusive
            }),
            
            # Pick score closest to scale midpoint
            baseline_score = map_dbl(scores_in_range, function(scores) {
              scores[which.min(abs(scores - scale_midpoint))]
            })
          ) |>
          mutate(
            # Calculate MDC at this baseline
            baseline_pct = baseline_score / max_sc,
            baseline_logit = qlogis(pmax(pmin(baseline_pct, 0.9999), 0.0001)),
            upper_logit = baseline_logit + mdc_log,
            upper_pct = plogis(upper_logit),
            upper_score = upper_pct * max_sc,
            
            mdc_items = ceiling(upper_score - baseline_score),
            mdc_percentage = upper_pct - baseline_pct,
            
            score_range = paste0(floor(min_score), "-", floor(max_score_range - 1))
          ) |>
          select(stratum, score_range, baseline_score, mdc_items, mdc_percentage)
      })
    ) |>
    select(stimuli, dx, max_score, baseline_data) |>
    unnest(baseline_data) |> 
    select(stimuli, dx, max_score, stratum, score_range, 
           `mdc (items)` = mdc_items, mdc_percentage)
}

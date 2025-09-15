calc_icc <- function(y, fixed_effects = NULL, facets, icc_facet, data, R = 1000, decimal = 2) {
  
  # Build formula
  fixed_part <- if(is.null(fixed_effects)) "1" else paste(fixed_effects, collapse = " + ")
  random_part <- paste0("(1|", facets, ")", collapse = " + ")
  model_formula <- as.formula(paste(y, "~", fixed_part, "+", random_part))
  
  # Fit model
  model <- lmer(model_formula, data = data)
  
  # Extract variance components
  vc <- VarCorr(model)
  sigma2_residual <- sigma(model)^2
  
  # Get all random effect variances
  sigma2_facets <- map_dbl(facets, ~as.numeric(vc[[.x]][1]))
  names(sigma2_facets) <- facets
  
  # Calculate total variance and percentages
  total_var <- sum(sigma2_facets) + sigma2_residual
  var_components <- c(sigma2_facets, residual = sigma2_residual)
  var_percentages <- (var_components / total_var) * 100
  
  # Calculate ICC
  sigma2_icc <- sigma2_facets[icc_facet]
  icc <- sigma2_icc / total_var
  
  # Bootstrap participants for CI
  boot_icc <- function(data, indices) {
    participants <- unique(data[[icc_facet]])[indices]
    d <- data |> filter(!!sym(icc_facet) %in% participants)
    tryCatch({
      m <- lmer(model_formula, data = d)
      vc <- VarCorr(m)
      s_icc <- as.numeric(vc[[icc_facet]][1])
      s_res <- sigma(m)^2
      s_others <- map_dbl(setdiff(facets, icc_facet), ~as.numeric(vc[[.x]][1]))
      s_icc / (s_icc + sum(s_others) + s_res)
    }, error = function(e) NA)
  }
  
  unique_levels <- unique(data[[icc_facet]])
  boot_results <- boot(data, boot_icc, R = R)
  ci <- boot.ci(boot_results, type = "perc")$percent[4:5]
  
  # Create formula string
  formula_str <- paste0("σ²", icc_facet, " / (", 
                        paste0("σ²", c(facets, "residual"), collapse = " + "), ")")
  
  list(
    icc_results = tibble(
      icc = icc,
      ci_lower = ci[1], 
      ci_upper = ci[2],
      formula = formula_str
    ),
    variance_components = tibble(
      component = names(var_percentages),
      variance = var_components,
      percentage = var_percentages
    )
  )
}

calc_icc_formula <- function(formula, icc_facet, data, R = 1000, decimal = 2) {
  
  # Extract facets from formula random effects
  random_terms <- findbars(formula)
  if(is.null(random_terms)) stop("No random effects found in formula")
  
  facets <- map_chr(random_terms, ~deparse(.x[[3]]))
  
  # Fit model and capture warnings
  warnings_initial <- character(0)
  model <- suppressWarnings(withCallingHandlers({
    lmer(formula, data = data, control = lmerControl(check.conv.singular = "ignore"))
  }, warning = function(w) {
    warnings_initial <<- c(warnings_initial, w$message)
  }))
  
  # Check for singular fit
  if(isSingular(model)) {
    warnings_initial <- c(warnings_initial, "Singular fit detected")
  }
  
  # Extract variance components
  vc <- VarCorr(model)
  sigma2_residual <- sigma(model)^2
  
  # Get all random effect variances
  sigma2_facets <- map_dbl(facets, ~as.numeric(vc[[.x]][1]))
  names(sigma2_facets) <- facets
  
  # Calculate total variance and percentages
  total_var <- sum(sigma2_facets) + sigma2_residual
  var_components <- c(sigma2_facets, residual = sigma2_residual)
  var_percentages <- (var_components / total_var)
  
  # Calculate ICC
  sigma2_icc <- sigma2_facets[icc_facet]
  icc <- sigma2_icc / total_var
  
  # Bootstrap participants for CI
  warnings_bootstrap <- character(0)
  boot_icc <- function(data, indices) {
    participants <- unique(data[[icc_facet]])[indices]
    d <- data |> filter(!!sym(icc_facet) %in% participants)
    tryCatch({
      m <- suppressWarnings(withCallingHandlers({
        lmer(formula, data = d, control = lmerControl(check.conv.singular = "ignore"))
      }, warning = function(w) {
        warnings_bootstrap <<- c(warnings_bootstrap, w$message)
      }))
      
      # Check for singular fit in bootstrap
      if(isSingular(m)) {
        warnings_bootstrap <<- c(warnings_bootstrap, "Singular fit detected")
      }
      
      vc <- VarCorr(m)
      s_icc <- as.numeric(vc[[icc_facet]][1])
      s_res <- sigma(m)^2
      s_others <- map_dbl(setdiff(facets, icc_facet), ~as.numeric(vc[[.x]][1]))
      s_icc / (s_icc + sum(s_others) + s_res)
    }, error = function(e) NA)
  }
  
  unique_levels <- unique(data[[icc_facet]])
  boot_results <- boot(data, boot_icc, R = R)
  ci <- boot.ci(boot_results, type = "perc")$percent[4:5]
  
  # Create warnings table
  all_warnings <- c(warnings_initial, warnings_bootstrap)
  unique_warnings <- unique(all_warnings)
  
  warnings_table <- if(length(unique_warnings) > 0) {
    tibble(
      warning = unique_warnings,
      count = map_int(unique_warnings, ~sum(all_warnings == .x))
    )
  } else {
    tibble(warning = character(0), count = integer(0))
  }
  
  # Create formula string
  formula_str <- paste0("σ²", icc_facet, " / (", 
                        paste0("σ²", c(facets, "residual"), collapse = " + "), ")")
  
  # Create results list
  results <- list(
    icc_results = tibble(
      icc = round(icc, decimal),
      ci_lower = ci[1],
      ci_upper = ci[2],
      ci_95 = paste0("[", round(ci[1], decimal), ", ", round(ci[2], decimal), "]"),
      formula = formula_str
    ),
    variance_components = tibble(
      component = names(var_percentages),
      variance = round(var_components, decimal),
      percentage = scales::label_percent(accuracy = (10^-decimal))(var_percentages)
    ),
    warnings = warnings_table,
    paramters = parameters::model_parameters(model)
  )
  
  # Print message if warnings detected
  if(nrow(warnings_table) > 0) {
    message("Model fit warnings detected, inspect $warnings for more information")
  }
  
  return(results)
}


calc_mdc_formula <- function(formula, icc_facet, data, error_type = "residual", decimal = 2) {
  
  # Extract facets from formula random effects
  random_terms <- findbars(formula)
  if(is.null(random_terms)) stop("No random effects found in formula")
  
  facets <- map_chr(random_terms, ~deparse(.x[[3]]))
  
  # Fit model and capture warnings
  warnings_initial <- character(0)
  model <- suppressWarnings(withCallingHandlers({
    lmer(formula, data = data, control = lmerControl(check.conv.singular = "ignore"))
  }, warning = function(w) {
    warnings_initial <<- c(warnings_initial, w$message)
  }))
  
  # Check for singular fit
  if(isSingular(model)) {
    warnings_initial <- c(warnings_initial, "Singular fit detected")
  }
  
  # Extract variance components
  vc <- VarCorr(model)
  sigma2_residual <- sigma(model)^2
  
  # Get all random effect variances
  sigma2_facets <- map_dbl(facets, ~as.numeric(vc[[.x]][1]))
  names(sigma2_facets) <- facets
  
  # Calculate total variance and percentages
  total_var <- sum(sigma2_facets) + sigma2_residual
  var_components <- c(sigma2_facets, residual = sigma2_residual)
  var_percentages <- (var_components / total_var)
  
  # Calculate error variance based on error_type
  if(error_type == "residual") {
    error_variance <- sigma2_residual
    error_components <- "residual"
  } else if(error_type == "non_icc") {
    non_icc_facets <- setdiff(facets, icc_facet)
    error_variance <- sigma2_residual + sum(sigma2_facets[non_icc_facets])
    error_components <- c(non_icc_facets, "residual")
  } else {
    stop("error_type must be 'residual' or 'non_icc'")
  }
  
  # Calculate MDC
  sem <- sqrt(error_variance)
  mdc <- 1.96 * sqrt(2) * sem
  
  # Create warnings table
  unique_warnings <- unique(warnings_initial)
  warnings_table <- if(length(unique_warnings) > 0) {
    tibble(
      warning = unique_warnings,
      count = map_int(unique_warnings, ~sum(warnings_initial == .x))
    )
  } else {
    tibble(warning = character(0), count = integer(0))
  }
  
  # Create results
  results <- list(
    mdc_results = tibble(
      mdc = round(mdc, decimal),
      sem = round(sem, decimal),
      error_variance = round(error_variance, decimal),
      error_type = error_type
    ),
    warnings = warnings_table
  )
  
  if(nrow(warnings_table) > 0) {
    message("Model fit warnings detected, inspect $warnings for more information")
  }
  
  return(results)
}

calc_decision_study <- function(model, reliability_facet = "participant", 
                                n_stimuli = 1:5, decimal = 2, R = 1000) {
  
  # Get original data and formula from model
  data <- model@frame
  formula <- formula(model)
  
  # Extract original variance components
  vc <- VarCorr(model)
  sigma2_residual <- sigma(model)^2
  sigma2_participant <- as.numeric(vc[[reliability_facet]][1])
  sigma2_stimuli <- as.numeric(vc[["stimuli"]][1])
  sigma2_interaction <- as.numeric(vc[["stimuli:participant"]][1])
  
  # Calculate original G coefficients
  original_g <- map_dbl(n_stimuli, function(n) {
    sigma2_participant / 
      (sigma2_participant + sigma2_stimuli/n +  sigma2_residual/n) #sigma2_interaction/n +
  })
  
  # Bootstrap function
  boot_g_coeff <- function(data, indices, formula, n_stimuli, reliability_facet) {
    
    # Resample participants
    participant_ids <- unique(data[[reliability_facet]])
    boot_participants <- sample(participant_ids, replace = TRUE)
    
    # Create bootstrap dataset
    boot_data <- map_dfr(boot_participants, function(pid) {
      data[data[[reliability_facet]] == pid, ]
    })
    
    # Refit model, handling potential failures
    tryCatch({
      boot_model <- lmer(formula, data = boot_data, 
                         control = lmerControl(check.conv.singular = "ignore"))
      
      # Extract variance components
      boot_vc <- VarCorr(boot_model)
      boot_sigma2_residual <- sigma(boot_model)^2
      boot_sigma2_participant <- as.numeric(boot_vc[[reliability_facet]][1])
      boot_sigma2_stimuli <- as.numeric(boot_vc[["stimuli"]][1])
      boot_sigma2_interaction <- as.numeric(boot_vc[["stimuli:participant"]][1])
      
      # Calculate G coefficients for all n
      map_dbl(n_stimuli, function(n) {
        boot_sigma2_participant / 
          (boot_sigma2_participant + boot_sigma2_stimuli/n + 
             + boot_sigma2_residual/n) #boot_sigma2_interaction/n 
      })
      
    }, error = function(e) {
      rep(NA, length(n_stimuli))
    })
  }
  
  # Run bootstrap
  boot_results <- replicate(R, boot_g_coeff(data, NULL, formula, n_stimuli, reliability_facet))
  
  # Calculate confidence intervals
  results <- map_dfr(seq_along(n_stimuli), function(i) {
    
    boot_values <- boot_results[i, ]
    boot_values <- boot_values[!is.na(boot_values)]  # Remove failed fits
    
    if(length(boot_values) > 0) {
      ci_lower <- quantile(boot_values, 0.025, na.rm = TRUE)
      ci_upper <- quantile(boot_values, 0.975, na.rm = TRUE)
      convergence_rate <- length(boot_values) / R
    } else {
      ci_lower <- ci_upper <- NA
      convergence_rate <- 0
    }
    
    tibble(
      n_stimuli = n_stimuli[i],
      g_coefficient = round(original_g[i], decimal),
      ci_lower = round(ci_lower, decimal),
      ci_upper = round(ci_upper, decimal),
      ci_95 = paste0("[", round(ci_lower, decimal), ", ", round(ci_upper, decimal), "]"),
      convergence_rate = round(convergence_rate, 2)
    )
  })
  
  return(results)
}


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
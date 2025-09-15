## Functions to compute ICCs for IRR 
## Author: Debby ten Hove
## UPDATED 20 June 2025 by Terrence D. Jorgensen: bug fixes on Lines 84:103
##          (change data$rater to data[[raters]], likewise for subjects,
##           calculate a single table, fix how (in)complete is checked)

## IF the packages below are not installed yet, install them first:
list.of.packages <- c("brms", "rstan", "lme4", "merDeriv", "car", "semTools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages) 

############################################################################
## FUNCTIONS TO COMPUTE KHAT AND Q
############################################################################
## Function to compute khat and Q from obs. design
computeQkhat <- function(data, subjects = "subjects", raters = "raters") {
  RR <- data[[raters]]
  SS <- data[[subjects]]
  tabRxS <- table(RR, SS)
  uSub <- ncol(tabRxS) # Number of unique subjects
  khat <- uSub / sum(1/colSums(tabRxS)) # harmonic mean number of raters per subject
  
  share <- 0 # initialize the proportion shared raters
  for (i in 1:uSub) {
    k_s <- colSums(tabRxS)[i]
    
    for (j in (1:uSub)[-i]) {
      k_sprime <- colSums(tabRxS)[j]
      k_s.sprime <- sum(RR[SS == i] %in% RR[SS == j])
      share <- share + (k_s.sprime / (k_s*k_sprime))/(uSub * (uSub-1))
    }
  }
  Q <- round(1/khat - share, 3)
  names(Q) <- "Q"
  
  return(list(Q = Q, khat = khat))
}

## Function to compute only khat (saves time when Q is not needed)
computeKhat <- function(data, subjects = "subjects", raters = "raters") {
  RR <- data[[raters]]
  SS <- data[[subjects]]
  tabRxS <- table(RR, SS)
  uSub <- ncol(tabRxS) # Number of unique subjects
  khat <- uSub / sum(1/colSums(tabRxS)) # harmonic mean number of raters per subject
  
  return(khat)
}
############################################################################
## FUNCTIONS TO ESTIMATE ICCs from CONTINUOUS DATA 
############################################################################

## MCMC with BRMS 
estICC_MCMC <- function(data, Y, subjects, raters, level, k, khat, Q){
  ## ICC and sigma names, for indexing and renaming
  ICCnames <- c("ICCa1", "ICCak", "ICCakhat",
                "ICCc1", "ICCck", "ICCqkhat")
  sigmanames <- c("S_s", "S_r", "S_sr") 
 
  ## Estimate model
  
  modForm <- paste(Y, "~ 1 + (1|", subjects, ") + (1 |", raters, ")")
  brmOUT <- brms::brm(as.formula(modForm), 
                      data   = data, 
                      warmup = 500, 
                      iter   = 1000, 
                      chains = 3, 
                      init  = "random",
                      cores  = 3)
  # IF non-converged: Return NA for everything 
  if(any( abs(brms::rhat(brmOUT)[2:4] - 1) > .1)){  
    return("Sorry, the model did not converge: Rhat > 1.1")
    } else {
    ## If converged: Give me results 
    ## Extract posterior distribution of SDs
    SDs <- rstan::extract(brmOUT$fit, c(paste0("sd_", subjects, "__Intercept"),
                                        paste0("sd_",   raters, "__Intercept"),
                                        "sigma"))
    # List SDs to later get MAPs
    names(SDs) <- c("SD_s", "SD_r", "SD_sr")
    
    ## Convert to variances 
    S_s <- SDs[["SD_s"]]^2
    S_r <- SDs[["SD_r"]]^2  
    S_sr <- SDs[["SD_sr"]]^2
    
    # List variances to later get MAPs
    sigmas <- list(S_s = S_s, S_r = S_r, S_sr = S_sr)
    
    ## Obtain ICCs
    ICCa1 <- as.numeric(S_s / (S_s + S_r + S_sr))
    ICCak <- as.numeric(S_s / (S_s + (S_r + S_sr)/k))
    ICCakhat <- as.numeric(S_s / (S_s + (S_r + S_sr)/khat))
    ICCc1 <- as.numeric(S_s / (S_s + S_sr))
    ICCck <- as.numeric(S_s / (S_s + S_sr/k))
    ICCqkhat <- as.numeric(S_s / (S_s + Q*S_r + S_sr/khat))
    
    
    ICCs <- list(ICCa1 = ICCa1, ICCak = ICCak, ICCakhat = ICCakhat,
                 ICCc1 = ICCc1, ICCck = ICCck, ICCqkhat = ICCqkhat)
    
    ## Confidence levels 
    ## Note. Percentiles-based BCIs based on hyperprior paper Ten Hove et al. (2021)
    ICC_cis <- unlist(lapply(ICCs, quantile, probs = c((1-level)/2, level + (1-level)/2)))
    ICC_cis <- matrix(ICC_cis, ncol = 2, byrow = T, dimnames = list(ICCnames, c("lower", "upper")))
    sigma_cis <- unlist(lapply(sigmas, quantile, probs = c((1-level)/2, level + (1-level)/2)))
    sigma_cis <- matrix(sigma_cis, ncol = 2, byrow = T, dimnames = list(sigmanames, c("lower", "upper")))

    ## SEs (posterior SDs)
    sigma_se <- unlist(lapply(sigmas, sd))
    names(sigma_se) <- paste0(names(sigmas), "_se")
    ICC_se <- unlist(lapply(ICCs, sd))
    names(ICC_se) <- paste0(names(ICCs), "_se")
    
    ## Point estimates (last, to not overwrite sigmas and ICCs sooner)
    # function to estimate posterior modes
    Mode <- function(x) {
      d <- density(x)
      d$x[which.max(d$y)]
    }
    # Point estimates (MAPs: Ten Hove et al., 2021) 
    est <- mapply(Mode, ICCs) # ICCs 
    ICCs = cbind(est, ICC_cis, ICC_se)
    
    est <- mapply(Mode, sigmas) # Variances
    variances = cbind(est, sigma_cis, sigma_se)
    
    # Return results
    OUT <- list(ICCs = ICCs,
                variances = variances,
                raterDesign = c("k" = k, "khat" = khat, "Q" = Q))
    return(OUT)
    }
  }

## MLE with LME4
estICC_LME4 <- function(data, Y, subjects, raters, level, k, khat, Q){
  ## ICCand sigma names, for indexing and renaming
  ICCnames <- c("ICCa1", "ICCak", "ICCakhat",
                "ICCc1", "ICCck", "ICCqkhat")
  sigmanames <- c("S_s", "S_r", "S_sr") 

  ## Define model
  modForm <- paste(Y, "~ 1 + (1|", subjects, ") + (1 |", raters, ")")
  ## Estimate model
  mod   <- lme4::lmer(as.formula(modForm), data = data)
  ## Check convergence
  checkConv <- function(mod) { 
    warn <- mod@optinfo$conv$lme4$messages
    !is.null(warn) && grepl('failed to converge', warn) 
  }
  if(checkConv(mod)){
    return("Sorry, the model did not converge.")
  } else {
    # If converged: Give me results
    ## Extract variances
    S_s  <- lme4::VarCorr(mod)[[subjects]][1, 1]  
    S_r  <- lme4::VarCorr(mod)[[raters]][1, 1]
    S_sr  <- sigma(mod)^2 
    ## Compute ICC point estimates
    ICCa1 <- S_s / (S_s + S_r + S_sr)
    ICCak <- S_s / (S_s + (S_r + S_sr)/k)
    ICCakhat <- S_s / (S_s + (S_r + S_sr)/khat)
    ICCc1 <- S_s / (S_s + S_sr) 
    ICCck <- S_s / (S_s + (S_sr)/k) 
    ICCqkhat <- S_s / (S_s + Q*S_r + S_sr/khat)
    
    ## List all
    sigmas <- c(S_s = S_s, S_r = S_r, S_sr = S_sr)
    ICCs <- c(ICCa1 = ICCa1, ICCak = ICCak, ICCakhat = ICCakhat,
              ICCc1 = ICCc1, ICCck = ICCck, ICCqkhat = ICCqkhat)
    print(ICCs)
    names(ICCs) <- ICCnames
    
    ## Asymptotic vcov matrix of sigmas
    suppressWarnings(ACOV <- merDeriv::vcov.lmerMod(mod, full = TRUE))
    Sidx <- grep(pattern = subjects, colnames(ACOV), fixed = TRUE) 
    Ridx <- grep(pattern = raters, colnames(ACOV), fixed = TRUE)
    SRidx <- which(colnames(ACOV) == "residual")
    idx      <- c(   Sidx  ,  Ridx  ,  SRidx  )
    newNames <- c("subject", "rater", "interaction")
    VCOV <- ACOV[idx, idx]
    dimnames(VCOV) <- list(newNames, newNames)
    vars <- c(subject = S_s, rater = S_r, interaction = S_sr)
    
    ## CIs and SEs of ICCs using asymptotic vcov matrix
    ## All info of all ICCs in one list
    ICCdefs <- c("subject / (subject + rater + interaction)", 
                 "subject / (subject + (rater + interaction)/k)",
                 "subject / (subject + (rater + interaction)/khat)", 
                 "subject / (subject + interaction)",
                 "subject / (subject + interaction/k)",
                 "subject / (subject + Q*rater + interaction/khat)"
    )
    names(ICCdefs) <- ICCnames
    
    
    ## Delta-method SEs of variances and ICCs
    ICCs_dm <- do.call("rbind", lapply(ICCdefs, FUN = function(x){
      car::deltaMethod(vars, vcov. = VCOV, level = level,g. = x)
    }))
    
    ICC_se <- ICCs_dm[,"SE"]
    names(ICC_se) <- paste0(ICCnames, "_se")
    sigma_se <- do.call("rbind", lapply(newNames, FUN = function(x){
      car::deltaMethod(vars, vcov. = VCOV, level = level,g. = x)
    }))$SE 
    names(sigma_se) <- paste0(names(sigmas), "_se")
    
    ## Monte-Carlo cis of variances and ICCs
    dimnames(VCOV) <- list(names(sigmas), names(sigmas))
    ICCs <- semTools::monteCarloCI(expr = c(ICCa1 = "S_s / (S_s + S_r + S_sr)", 
                                                  ICCak = paste0("S_s / (S_s + (S_r + S_sr)/", k, ")"),
                                                  ICCakhat = paste0("S_s / (S_s + (S_r + S_sr)/", khat, ")"), 
                                                  ICCc1 = "S_s / (S_s + S_sr)",
                                                  ICCck = paste0("S_s / (S_s + S_sr/", k, ")"),
                                                  ICCqkhat = paste0("S_s / (S_s + ", Q, "*S_r + S_sr/", khat, ")")),
                                         coefs = sigmas, ACM = VCOV)
    sigmas <- semTools::monteCarloCI(expr = c(S_s = 'S_s', S_r = "S_r", S_sr = "S_sr"),
                                     coefs = sigmas, ACM = VCOV)
    
    # Return results
    OUT <- list(ICCs = cbind(ICCs, ICC_se),
                variances = cbind(sigmas, sigma_se),
                raterDesign = c("k" = k, "khat" = khat, "Q" =  Q))
    return(OUT)
  }
}

############################################################################
## FUNCTIONS TO ESTIMATE ICCs from BINARY DATA 
############################################################################

## MCMC with brms
estICC_MCMC_bin <- function(data, Ybin, subjects, raters, level, k, khat, Q){
  
  ## ICC and sigma names, for indexing and renaming
  ICCnames <- c("ICCa1", "ICCak", "ICCakhat",
                "ICCc1", "ICCck", "ICCqkhat")
  sigmanames <- c("S_s", "S_r", "S_sr") 

  
  ## Estimate model  
  fam <- brms::bernoulli(link = "logit")
  modForm <- paste(Ybin, "~ 1 + (1|", subjects, ") + (1 |", raters, ")")
  brmOUT <- brms::brm(as.formula(modForm), 
                      data   = data, 
                      family = fam,
                      warmup = 500, 
                      iter   = 1000, 
                      chains = 3, 
                      init  = "random",
                      cores  = 3)
  
  # IF non-converged: Return NA for everything 
  if(any( abs(brms::rhat(brmOUT)[2:4] - 1) > .1)){ 
    return("Sorry, the model did not converge.")
    } else {
    ## If converged: Give me results 
    ## Extract posterior distribution of SDs
    SDs <- rstan::extract(brmOUT$fit, c(paste0("sd_", subjects, "__Intercept"),
                                        paste0("sd_",   raters, "__Intercept")))
    names(SDs) <- c("SD_s", "SD_r")
    ## Transform to variance components
    S_s <- SDs[["SD_s"]]^2
    S_r <- SDs[["SD_r"]]^2  
    sigmas <- list(S_s = S_s, S_r = S_r)
    S_sr <- pi^2 / 3
      
    ## Obtain ICCs
    ICCa1 <- as.numeric(S_s / (S_s + S_r + S_sr))
    ICCak <- as.numeric(S_s / (S_s + (S_r + S_sr)/k))
    ICCakhat <- as.numeric(S_s / (S_s + (S_r + S_sr)/khat))
    ICCc1 <- as.numeric(S_s / (S_s + S_sr))
    ICCck <- as.numeric(S_s / (S_s + S_sr/k))
    ICCqkhat <- as.numeric(S_s / (S_s + Q*S_r + S_sr/khat))
    
    ## List posteriors of ICCs 
    ICCs <- list(ICCa1 = ICCa1, ICCak = ICCak, ICCakhat = ICCakhat,
                 ICCc1 = ICCc1, ICCck = ICCck, ICCqkhat = ICCqkhat)
    
    ## Confidence levels 
    ## Note. I Changed HDI to percentiles-based BCIs based on hyperprior paper 
    ICC_cis <- unlist(lapply(ICCs, quantile, probs = c((1-level)/2, level + (1-level)/2)))
    ICC_cis <- matrix(ICC_cis, ncol = 2, byrow = T, dimnames = list(ICCnames, c("lower", "upper")))
    sigma_cis <- c(unlist(lapply(sigmas, quantile, probs = c((1-level)/2, level + (1-level)/2))), NA, NA)
    sigma_cis <- matrix(sigma_cis, ncol = 2, byrow = T, dimnames = list(sigmanames, c("lower", "upper")))

    ## SEs (posterior SDs)
    sigma_se <- unlist(lapply(sigmas, sd))
    names(sigma_se) <- paste0(names(sigmas), "_se")
    sigma_se  <- c(sigma_se, S_sr_se = NA) # Add NA for SE for S_sr
    ICC_se <- unlist(lapply(ICCs, sd))
    names(ICC_se) <- paste0(names(ICCs), "_se")
    
    ## Point estimates (last, to not overwrite sigmas and ICCs sooner)
    # function to estimate posterior modes
    Mode <- function(x) {
      d <- density(x)
      d$x[which.max(d$y)]
    }
    # MAPs 
    sigmas <- list(S_s = S_s, S_r = S_r)
    sigmas <- c(mapply(Mode, sigmas), S_sr = S_sr) # Variances
    ## YIELD MAPS FOR sigmas < Exclude S_sr because not estimated from the data      
    ICCs <- mapply(Mode, ICCs) # ICCs 
    
    # Return results
    OUT <- list(ICCs = cbind(ICCs, ICC_cis, ICC_se),
                variances = cbind(sigmas, sigma_cis, sigma_se),
                raterDesign = c("k" = k, "khat" = khat, "Q" = Q))
    return(OUT)
  }
}

## MLE with LME4
estICC_LME4_bin <- function(data, Ybin, subjects, raters, level, k, khat, Q){
  ## ICC and sigma names, for indexing and renaming
  ICCnames <- c("ICCa1", "ICCak", "ICCakhat",
                "ICCc1", "ICCck", "ICCqkhat")
  sigmanames <- c("S_s", "S_r", "S_sr") 
  
  ## Define model
  modForm <- paste(Ybin, "~ 1 + (1|", subjects, ") + (1 |", raters, ")")
  ## Estimate model
  mod   <- lme4::glmer(as.formula(modForm), data = data, family = binomial("logit"))
  ## Check convergence
  checkConv <- function(mod) { 
    warn <- mod@optinfo$conv$lme4$messages
    !is.null(warn) && grepl('failed to converge', warn) 
  }
  if(checkConv(mod)){
    # If nonconverged: Return NAs for everything
    return("Sorry, the model did not converge.")
  } else {
    # If converged: Give me results
    ## Extract variances
    S_s  <- lme4::VarCorr(mod)[[subjects]][1, 1]  
    S_r  <- lme4::VarCorr(mod)[[raters]][1, 1]
    S_sr  <- 1 
    
    ## Compute ICC point estimates
    ICCa1 <- S_s / (S_s + S_r + S_sr)
    ICCak <- S_s / (S_s + (S_r + S_sr)/k)
    ICCakhat <- S_s / (S_s + (S_r + S_sr)/khat)
    ICCc1 <- S_s / (S_s + S_sr) 
    ICCck <- S_s / (S_s + (S_sr)/k) 
    ICCqkhat <- S_s / (S_s + Q*S_r + S_sr/khat)
    
    ## List all 
    sigmas <- c(S_s = S_s, S_r = S_r, S_sr = S_sr)
    ICCs <- c(ICCa1 = ICCa1, ICCak = ICCak, ICCakhat = ICCakhat,
              ICCc1 = ICCc1, ICCck = ICCck, ICCqkhat = ICCqkhat)
    names(ICCs) <- ICCnames
    
    ## list results
    OUT <- list(ICCs = ICCs, variances = sigmas,  
                raterDesign = c("k" = k, "khat" = khat, "Q" = Q)) 
    return(OUT)
  }
}



############################################################################
## COMBINED FUNCTIONS FOR SIMULATION
############################################################################
estICCs <- function(data, 
                    Rep = Rep,
                    Y = "Y", 
                    subjects = "subjects", 
                    raters = "raters", 
                    level = .95,
                    k = NULL, 
                    khat = NULL, 
                    Q = NULL, 
                    estimator = "MLE",
                    response = "continuous") {
  
  ## Number of raters
  if(is.null(k)){
    k <- length(unique(data[[raters]]))
  }
  
  ## Check type of design
  tabSxR <- table(data[[subjects]], by = data[[raters]])
  
  # Balanced or unbalanced (== number raters per subject)
  if (length(unique(rowSums(tabSxR))) == 1) {
    balanced <- TRUE
  } else {
    balanced <- FALSE
  }
  # Complete or incomplete (number of subjects per rater == number of subjects)
  if (all(colSums(tabSxR) == length(unique(data[[subjects]])))) {
    complete <- TRUE 
  } else {
    complete <- FALSE
  }
  # Two-Way (crossed design)?  One-Way (nested) implies each rater rates 1 subject
  if (all(colSums(tabSxR) == 1)) {
    twoWay <- FALSE
  } else {
    twoWay <- TRUE
  }
  
  if(is.null(khat) | is.null(Q)){
    ## Decide on values for khat and q 
    if(balanced == T & complete == T){ 
      khat <- k
      Q <- 0 
    } else {
      if(balanced == T & complete == F){
        khat <- unique(rowSums(table(data[[subjects]], by = data[[raters]])))
        Q <- computeQkhat(data, subjects = subjects, raters = raters)$Q[["Q"]]
      } else {
        if(balanced == F & complete == F){
          Qkhat <- computeQkhat(data, subjects = subjects, raters = raters)
          khat <- Qkhat$khat
          Q <- Qkhat$Q[["Q"]]
        } else {
          if(twoWay == F){
            khat <- computeKhat(data, subjects = subjects, raters = raters)
            Q <- 1/k # But not needed, since sigmaR cannot be distinguished
          }
        }
      }
    }
  }
  

  ### ESTIMATE ICCs: Continuous data
  if(response == "continuous"){
    if(estimator == "MCMC"){
      ### Estimate two-way model using MCMC with BRMS, silent = TRUE)
      OUT <- estICC_MCMC(data = data, Y = Y, subjects = subjects, raters = raters, 
                        level = level, k = k, khat = khat, Q = Q)
      }
    if(estimator == "MLE"){
      ### MLE with lme4 (random-effects model)
      OUT <- estICC_LME4(data = data, Y = Y, subjects = subjects, raters = raters, 
                             level = level, k = k, khat = khat, Q = Q)
      }
    }
  
  ### ESTIMATE ICCs: Binary data
  if(response == "binary"){
    if(estimator == "MCMC"){
      ### MCMC with BRMS
      OUT <- estICC_MCMC_bin(data = data, Ybin = Ybin, subjects = subjects, raters = raters, 
                            level = level, k = k, khat = khat, Q = Q)
      }
    if(estimator == "MLE"){
      ### MLE with lme4 (random-effects model)
      OUT <- estICC_LME4_bin(data = data, Ybin = Ybin, subjects = subjects, raters = raters, 
                            level = level, k = k, khat = khat, Q = Q)
      }
  }
  return(OUT)
  
}


# Details ----
#' 01d_outliers.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: 
#' ++ Original diagnostics code from (Sarakinis, 2024)
#' Date: 
#' - Origin: 6 Jul 2025
#' - Updated: 20 Feb 2026
#' 
#' Content: Looking for outliers
#' ++ part II, diagnostic plots
#' 
#' Note:
#' + xlsx prepped for analysis in 00_prep_spot
#' -----------

# load:
source('helpers/microHelp.R')
micro <- readxl::read_xlsx("data/output/icpms_clean.xlsx") 

# function:
get_outliers <- function(data, response, predictors) {
  outlier_summary <- data.frame()
  
  for (col in response) {
    cat("\n=== Model Outliers for", col, "===\n")
    
    # Create model formula
    formula_str <- paste(col, "~", paste(predictors, collapse = " + "))
    
    # Fit the model
    model <- lm(as.formula(formula_str), data = data)
    
    # Create diagnostic plots
    par(mfrow = c(2, 2))
    plot(model, main = paste("Model diagnostics:", col))
    par(mfrow = c(1, 1))
    
    # Identify influential outliers
    cooks_d <- cooks.distance(model)
    n <- nobs(model)
    leverage <- hatvalues(model)
    
    # Standard thresholds
    cooks_threshold <- 4/n
    leverage_threshold <- 2 * length(coef(model))/n
    
    influential_points <- which(cooks_d > cooks_threshold | leverage > leverage_threshold)
    
    cat("Model summary:\n")
    cat("Observations:", n, "\n")
    cat("R-squared:", round(summary(model)$r.squared, 3), "\n")
    cat("Influential points (Cook's D or leverage):", length(influential_points), "\n")
    
    if (length(influential_points) > 0) {
      cat("Row indices of influential points:", influential_points, "\n")
    }
    
    # Store summary
    outlier_summary <- rbind(outlier_summary, data.frame(
      Variable = col,
      N_obs = n,
      R_squared = summary(model)$r.squared,
      N_influential = length(influential_points),
      Max_cooks_D = max(cooks_d, na.rm = TRUE)
    ))
    
    if (interactive()) {
      readline(prompt = "Press [Enter] to continue to next variable...")
    }
  }
  
  return(outlier_summary)
}

# application:
get_outliers(data = micro,
            response = el,
            c("Species")) # not normal - lets fix and run again!
#smart!

get_outliers(data = micro, # go again
             response = ln_rats,
             c("Species", "Region")) # lol Ni and Pb!

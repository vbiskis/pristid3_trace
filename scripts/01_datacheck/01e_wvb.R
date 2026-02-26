# Details ----
#' 01e_wvb.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 
#' - Origin: 18 Aug 2024
#' - Updated: 20 Feb 2026
#' 
#' Content: 
#' + checking how signature has changed bleach vs water preparation
#' 
#' Note:
#' + xlsx prepped for analysis in 00_prep_spot
#' + style supplied in microHelp.R
#' -----------

source('helpers/microHelp.R')
micro <- readxl::read_xlsx("data/output/icpms_clean.xlsx") 

library(lme4)
library(emmeans)

BWset <- micro %>%
  group_by(ID) %>%  
  filter(all(c("B", "W") %in% BW),
         In_Ca <= 0.5,
         Type == "Tooth",
         Cut != "Trans") %>% #i bet its the verts
  ungroup() 

results_list <- list() #of course initialise a list to store results in

#check 1 - t test B to W----

for (element in ln_rats) {
  cat("Processing:", element, "\n")
  
  # Check if column exists
  if (!element %in% colnames(BWset)) {
    cat("Column", element, "not found in data\n")
    next
  }
  
  # Create separate datasets for B and W
  B_data <- BWset %>% 
    dplyr::filter(BW == "B") %>%
    dplyr::select(ID, Type, Cut, Zone_2, Value = !!sym(element))
  
  W_data <- BWset %>% 
    dplyr::filter(BW == "W") %>%
    dplyr::select(ID, Type, Cut, Zone_2, Value = !!sym(element))
  
  # Join them to create paired data
  paired_data <- B_data %>%
    inner_join(W_data, by = c("ID", "Type", "Cut", "Zone_2"), 
               suffix = c("_B", "_W")) %>%
    filter(!is.na(Value_B) & !is.na(Value_W))
  
  cat("Paired observations (same location, B vs W):", nrow(paired_data), "\n")
  
  if (nrow(paired_data) >= 2) {
    tryCatch({
      # Paired t-test: B condition vs W condition for same anatomical location
      t_result <- wilcox.test(paired_data$Value_B, paired_data$Value_W, paired = TRUE)
      
      results_list[[element]] <- tibble(
        Element = element,
        n_pairs = nrow(paired_data),
        p_value = t_result$p.value,
        statistic = as.numeric(t_result$statistic),
        conf_int_low = t_result$conf.int[1],
        conf_int_high = t_result$conf.int[2],
        mean_diff = t_result$estimate,  # Mean difference (B - W)
        mean_B = mean(paired_data$Value_B),
        mean_W = mean(paired_data$Value_W)
      )
      
      cat("  Success! n =", nrow(paired_data), "pairs, p =", round(t_result$p.value, 4), "\n")
      
    }, error = function(e) {
      cat("Error with", element, ":", e$message, "\n")
    })
  } else {
    cat("  Insufficient paired data for", element, "\n")
  }
}

# Results
if (length(results_list) > 0) {
  results_summary <- bind_rows(results_list)
  print(results_summary, n = Inf)
} else {
  cat("No results generated.\n")
}

#check 2 - t test for var between zones----
results_list <- list()

for (element in ln_rats) {
  cat("Processing:", element, "\n")
  
  if (!element %in% colnames(BWset)) {
    cat("Column", element, "not found in data\n")
    next
  }
  
  # Test variability within each BW condition separately
  for (bw_condition in c("B", "W")) {
    cat("  Processing", bw_condition, "condition\n")
    
    # Get all data for this BW condition
    bw_data <- BWset %>%
      dplyr::filter(BW == bw_condition) %>%
      dplyr::select(ID, Type, Cut, Zone_2, Value = !!sym(element)) %>%
      dplyr::filter(!is.na(Value))
    
    # Find locations with exactly 2 replicates for paired comparison
    paired_replicates <- bw_data %>%
      group_by(ID, Type, Cut, Zone_2) %>%
      dplyr::filter(n() == 2) %>%  # Only groups with exactly 2 measurements
      arrange(ID, Type, Cut, Zone_2) %>%
      summarise(
        rep1 = Value[1],
        rep2 = Value[2],
        .groups = 'drop'
      ) %>%
      filter(!is.na(rep1) & !is.na(rep2))
    
    cat("    Paired replicates found:", nrow(paired_replicates), "\n")
    
    if (nrow(paired_replicates) >= 2) {
      tryCatch({
        # Paired t-test: replicate 1 vs replicate 2 within same condition
        t_result <- t.test(paired_replicates$rep1, paired_replicates$rep2, paired = TRUE)
        
        results_list[[paste(element, bw_condition, sep = "_")]] <- tibble(
          Element = element,
          BW_condition = bw_condition,
          Comparison = paste(bw_condition, "vs", bw_condition, "(replicates)"),
          n_pairs = nrow(paired_replicates),
          p_value = t_result$p.value,
          statistic = as.numeric(t_result$statistic),
          conf_int_low = t_result$conf.int[1],
          conf_int_high = t_result$conf.int[2],
          mean_diff = t_result$estimate,  # Mean difference between replicates
          mean_rep1 = mean(paired_replicates$rep1),
          mean_rep2 = mean(paired_replicates$rep2)
        )
        
        cat("    Success! n =", nrow(paired_replicates), "pairs, p =", round(t_result$p.value, 4), "\n")
        
      }, error = function(e) {
        cat("    Error with", element, bw_condition, ":", e$message, "\n")
      })
    } else {
      cat("    Insufficient replicate pairs for", element, bw_condition, "\n")
    }
  }
}

# Results
if (length(results_list) > 0) {
  results_summary <- bind_rows(results_list)
  print(results_summary, n = Inf)
} else {
  cat("No results generated.\n")
}

#lmer!
model_li <- lmer(lnLi ~ BW + Zone_2 + Species + (1|ID), data = BWset)
summary(model_li)
anova(model_li)

model_na <- lmer(lnNa ~ BW + Zone_2 + Species + (1|ID), data = BWset)
summary(model_na) #yes it is!
anova(model_na)

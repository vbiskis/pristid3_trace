# Details ----
#' 01f_data_filter.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 
#' - Origin: 24 Jan 2024
#' - Updated: 20 Feb 2026
#' 
#' Content: 
#' need to filter E for only one of each before going further.
#' 
#' Note:
#' + xlsx prepped for analysis in 00_prep_spot
#' + style supplied in microHelp.R
#' -----------

source('helpers/microHelp.R')
micro <- readxl::read_xlsx("data/output/icpms_clean.xlsx") 
library(writexl)

microE <- micro[micro$Zone_2 %in% c("Edge"),] # grouped with both edges together as base

edge_sum <- microE %>%
  group_by(ID) %>%
  summarise(
    has_E1 = any(Zone == "EdgeE1"),
    has_E2 = any(Zone == "EdgeE2"),
    .groups = 'drop'
  ) %>%
  mutate(
    category = case_when(
      has_E1 & has_E2 ~ "Both E1 and E2",
      has_E1 & !has_E2 ~ "E1 only",
      !has_E1 & has_E2 ~ "E2 only", 
      !has_E1 & !has_E2 ~ "Neither E1 nor E2"
    )
  )

edge_counts <- edge_sum %>%
  count(category) %>%
  arrange(desc(n))

print(edge_counts)
edge_sum %>%
  filter(has_E1 & !has_E2) %>%
  pull(ID)

# ah! it's just a bunch of the CEQ region + one pristis and one v important NQ anoxy
# does it even matter - lets see:

library(broom)

results_list <- list()

for (element in ln_rats){
  
    E1_data <- microE %>% 
    filter(Zone == "EdgeE1") %>%
    select(ID, Value = !!sym(element)) %>%
    filter(!is.na(Value))
  
    E2_data <- microE %>% 
    filter(Zone == "EdgeE2") %>%
    select(ID, Value = !!sym(element)) %>%
    filter(!is.na(Value))
  
  # PAIRED
  paired_data <- E1_data %>%
    inner_join(E2_data, by = "ID", suffix = c("_E1", "_E2")) %>%
    filter(!is.na(Value_E1) & !is.na(Value_E2))
  
  if (nrow(paired_data) >= 2) {
    paired_test <- t.test(paired_data$Value_E1, 
                          paired_data$Value_E2, paired = TRUE)
    
    results_list[[paste(element, "paired", sep = "_")]] <- tibble(
      Element = element,
      Test_type = "Paired t-test",
      n_pairs = nrow(paired_data),
      mean_E1 = mean(paired_data$Value_E1),
      mean_E2 = mean(paired_data$Value_E2),
      mean_diff = paired_test$estimate,
      p_value = paired_test$p.value,
      statistic = paired_test$statistic,
      conf_int_low = paired_test$conf.int[1],
      conf_int_high = paired_test$conf.int[2]
    )
  }

  # UNPAIRED
  if (nrow(E1_data) >= 2 && nrow(E2_data) >= 2) {
    unpaired_test <- t.test(E1_data$Value, E2_data$Value, paired = FALSE)
    
    results_list[[paste(element, "unpaired", sep = "_")]] <- tibble(
      Element = element,
      Test_type = "Unpaired t-test",
      n_E1 = nrow(E1_data),
      n_E2 = nrow(E2_data),
      mean_E1 = mean(E1_data$Value),
      mean_E2 = mean(E2_data$Value),
      mean_diff = unpaired_test$estimate[1] - unpaired_test$estimate[2],
      p_value = unpaired_test$p.value,
      statistic = unpaired_test$statistic,
      conf_int_low = unpaired_test$conf.int[1],
      conf_int_high = unpaired_test$conf.int[2]
    )
  }
}

# Combine
if (length(results_list) > 0) {
  all_results <- bind_rows(results_list)
    all_results <- all_results %>%
    group_by(Test_type) %>%
    mutate(p_adjusted = p.adjust(p_value, method = "fdr")) %>%
    ungroup()
  
  print(all_results, n = Inf)
  
  significant_results <- all_results %>%
    filter(p_value < 0.05) %>%
    arrange(Test_type, p_value)
  
    print(significant_results, n = Inf)
  
} else {
  cat("No results generated - check data structure\n")
}

#retain those 2xE1 for low sample size but otherwise chuck the extras - no need.
e1pals <- c("CB001", "QBFP004468")

microET <- microE %>% 
  group_by(ID) %>% 
  filter(Zone == "EdgeE2" | (Zone == "EdgeE1" & !any(Zone == "EdgeE2"))) %>%
  ungroup()

write_xlsx(microET, path = "data/output/microE.xlsx" ) # kept

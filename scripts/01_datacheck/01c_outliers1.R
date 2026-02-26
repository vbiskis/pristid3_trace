# Details ----
#' 01c_outliers.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 
#' - Origin: 23 Aug 2024
#' - Updated: 20 Feb 2026
#' 
#' Content: Looking for outliers 
#' + first visualisation (a)
#' 
#' Note:
#' + xlsx prepped for analysis in 00_prep_spot
#' + style supplied in microHelp.R
#' -----------

source('helpers/microHelp.R')
set_theme(mytheme)
micro <- readxl::read_xlsx("data/output/icpms_clean.xlsx") 

#time to write a code for outliers

for(element in el){
  outliers <- micro %>%
    group_by(Species, Region, Zone) %>%
    mutate(mean = mean(.data[[element]], na.rm = TRUE),
           sd = sd(.data[[element]], na.rm = TRUE),
           is_out = abs(.data[[element]] - mean) > 3 * sd) 
  # Outliers defined as 3 SD from mean
  
  z <- ggplot(micro, 
              aes(Species,
                  .data[[element]],
                  color = Zone)) +
    facet_wrap(~Region, scales = "free") +
    geom_boxplot(position = position_dodge(width = 0.8),
                 outlier.colour = "red", 
                 outlier.shape = 16, 
                 outlier.size = 2) +
    geom_jitter(position = position_dodge(width = 0.8)) +
    geom_text(data = outliers %>% filter(is_out), 
              aes(label = File_ID),
              position = position_dodge(width = 0.8), 
              vjust = -0.5, 
              size = 3, 
              color = "black") + 
    labs(title = "El var by Species/Region", 
         x= bquote("Species"),
         y= paste("Mean", gsub("_", "/", element), "Ratio")) + 
    theme(axis.text.x = element_text(face = "italic"))
  print(z)
}

# a few standing out here. can remove 3SD values and look.
micro_cln <- micro # new one so dont mess up old guy

# Loop through each element in el
for(element in ln_rats) {
  micro_cln <- micro_cln %>% #same thing but now remove
    group_by(Species, Region, Zone) %>%
    mutate(
      mean = mean(.data[[element]], na.rm = TRUE),
      sd = sd(.data[[element]], na.rm = TRUE),
      is_out = abs(.data[[element]] - mean) > 3 * sd,  # Outliers defined as > 3 SD from mean
      !!element := ifelse(is_out, NA, .data[[element]])  # Replace outliers with NA
    ) %>%
    ungroup()
  
  micro_cln <- micro_cln %>%  # those columns were never here
    dplyr::select(-mean, -sd, -is_out)
}

#n is same, didnt even throw any out
for(element in ln_rats) {
  outliers <- micro_cln %>%
    group_by(Species, Region, Zone) %>%
    mutate(
      mean = mean(.data[[element]], na.rm = TRUE),
      sd = sd(.data[[element]], na.rm = TRUE),
      is_out = abs(.data[[element]] - mean) > 3 * sd  # Outliers defined as 3 SD from mean
    )
  
  z <- ggplot(micro_cln, 
              aes(Species, .data[[element]], color = Zone)) +
    facet_wrap(~Region, scales = "free") +
    geom_boxplot(position = position_dodge(width = 0.8),
                 outlier.colour = "red", 
                 outlier.shape = 16, 
                 outlier.size = 2) +
    geom_jitter(position = position_dodge(width = 0.8)) +
    geom_text(data = outliers %>% filter(is_out), 
              aes(label = File_ID),
              position = position_dodge(width = 0.8), 
              vjust = -0.5, 
              size = 3, 
              color = "black") + 
    labs(title = "El var by Species/Region", 
         x= "Species",
         y= paste("Mean", gsub("_", "/", element), "Ratio")) + 
    theme_minimal() +
    theme(axis.text.x = element_text(face = "italic"))
  print(z)
}

#really just created new outliers. Will see if we want to use this...
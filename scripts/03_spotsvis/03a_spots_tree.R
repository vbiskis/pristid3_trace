# Details ----
#' 03_spots_tree.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 18 Aug 2024
#' 
#' Content: 
#' looking at separation of signatures across key terms (vis)
#' + Update 20 Feb 2026:
#' ++ made simplified ctree + plots for repetition 
#' ++ looped it
#' ++ output autosaved
#' -----------

source('helpers/microHelp.R')
set_theme(mytheme)
micro <- read_xlsx('data/output/icpms_clean.xlsx')
str(micro)
micro <- micro %>% 
  mutate(across(where(is.character) & !all_of(c("ID", "File_ID")), as.factor))
#das better

#TREE----
library(partykit)

for(element in ln_rats){
    png(paste0("figures/ext/trees/", 
               element, 
               "plot.png"), 
        width = 10, height = 6, units = "in", res = 300)
  
    tree(element, dat = micro)
    
    dev.off()
}

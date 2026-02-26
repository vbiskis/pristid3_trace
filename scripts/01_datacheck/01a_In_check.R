# Details ----
#' 01a_In_check.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 
#' - Origin: 24 Jan 2024
#' - Updated: 20 Feb 2026
#' 
#' Content: Quality Control - First Look
#' + Indium vs other trace metals / linearity
#' + trace metals could just be higher due to being components of the epoxy
#' ++ also check whether IDs mistyped or didn't run (contaminated all)
#' 
#' Note:
#' + xlsx prepped for analysis in 00_prep_spot
#' + style supplied in microHelp.R
#' -----------

source('MicroHelp.R')
set_theme(mytheme)
micro <- readxl::read_xlsx("data/output/icpms_clean.xlsx") 
sawinf <- readxl::read_xlsx("data/Saw_Summ.xlsx", sheet = "Summary")

#check first if IDs mistyped or didn't run (contaminated all)
idmicro <- unique(micro$ID) # looking at all micro, unfiltered by In (inc V)
idsawinf <- unique(sawinf$ID)
setdiff(idsawinf, idmicro) # IDs in sawinf, but not micro, filtered out bc year

#Scatter----
#allspots
for(element in ln_cols){
  p <- ggplot(micro, 
         aes(x=lnIn, 
             y = .data[[element]],
             color = Zone,
             shape = Type)) +
    geom_point() +
    labs(title = paste("lnIn vs.", gsub("_", "/", element)), 
         x= bquote("lnIn Ratio"),
         y= paste(gsub("_", "/", element), "Ratio"))
  print(p)
}

#looks good
#now will investigate the ceiling and floor values of each trace
#filtered spots

#BoxPlots----
for(element in el){
  
  # filter - remove SC saws
  MicroKnown <- subset(MicroFil, r <= 50)
  
  q <- ggplot(MicroKnown, 
              aes(Species,
                  .data[[element]],
                  color = Zone)) +
    facet_wrap(~Region, scales = "free") +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    geom_jitter(position = position_dodge(width = 0.8)) +
    labs(title = "El var by Species/Region", 
         x= bquote("Species"),
         y= paste("Mean", gsub("_", "/", element), "Ratio")) + 
    theme(axis.text.x = element_text(face = "italic"))
  print(q)
}

#can see there are some ceiling values, will use plotly to identify

library(plotly)
# plotly won't do these all at once as is. 
# Needs encouragement (storage into a list), then we can call them

plot_list <- list()

for(element in el){
  r <- ggplot(MicroFil, 
              aes(Species,
                  .data[[element]],
                  color = Zone)) +
    facet_wrap(~Region) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    geom_jitter(aes(shape = ID), 
                position = position_dodge(width = 0.8)) +
    labs(title = "El var by Species/Region", 
         x= bquote("Species"),
         y= paste("Mean", gsub("_", "/", element), "Ratio")) + 
    theme(axis.text.x = element_text(face = "italic"))
  plot_list[[element]] <- ggplotly(r)
}

plot_list[1] #Li
plot_list[2] #B - all samples of CNSAQ1 look like that...?
plot_list[3] #Mg
plot_list[4] #Mn
plot_list[5] #Ni - maybe remove the 480s and see how it compares.
plot_list[6] #Cu - outlier 3150
plot_list[7] #Zn - cool trend, will investigate
plot_list[8] #Sr - great
plot_list[9] #Ba - yup
plot_list[10] #Pb - huh such a small amount is that real? 3151


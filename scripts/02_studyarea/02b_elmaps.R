# Details ----
#' 02b_elmaps.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 18 Aug 2024
#' - Updated: 20 Feb 2026
#' 
#' Content: 
#' + expl. look at element dist across study area
#' update: pulled out mapping libraries / formatting
#' ++ made source code, style supplied in microHelp.R
#' -----------

source('helpers/mapHelp.R')
source('helpers/microHelp.R')
set_theme(mytheme)
microE <- read_xlsx('data/output/microE.xlsx')
library(writexl)

# LODs input in place of neg vals
# high In samples excluded, and lowest value chosen per site.
# each point reps one tooth

# Element Variation----

elmaps <- list() # call them independently

for (element in el) {
  El_coords <- cbind(microE$Species, microE$LAT, microE$LONG, microE[[element]])
  crd_to_xy(El_coords, name = "ptsel")
  ptsel$V4 <- as.numeric(as.character(ptsel$V4))
  
  el_lab <- gsub("_", "/", element)
  ptsel <- ptsel[order(ptsel$V4), ]
  
  elmap <- ggplot() + 
    annotation_map_tile('cartolight') +
    geom_spatial_point(crs = 4326) +
    geom_sf(data = ptsel, aes(colour = V4), show.legend = TRUE) +  
    scale_color_viridis_c(name = el_lab,
                          trans = "log") +  # Continuous color scale 
    facet_grid(~V1) +  # Facet by species
    labs(x = "Longitude", y = "Latitude")
  
  elmaps[[element]] <- elmap
}

print(elmaps[[1]]) #no mines known
print(elmaps[[2]]) #no mines known
print(elmaps[[3]]) #no mines known
print(elmaps[[4]]) #big ass mine in Stanage - interesting the high dot there
print(elmaps[[5]]) #no mines known
print(elmaps[[6]]) #Ni mine in area of large point
print(elmaps[[7]]) #Cu mines in this area
print(elmaps[[8]]) #Zn outlier??
print(elmaps[[9]])
print(elmaps[[11]])
print(elmaps[[12]])
# Details ----
#' 07_env_map.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 17 Jun 2026
#' 
#' Content: 
#' + mapping study areas from DAF and DETSI
#' + prelim identification of trends
#' -----------

source('helpers/mapHelp.R')
source('helpers/microHelp.R')
set_theme(mytheme)

# Daf first----
daf <- read_xlsx('data/output/daf_clean.xlsx')
des <- read_xlsx('data/output/des_clean.xlsx')

# Let's look at maps of variation quick ----

elmaps <- list() 
for (element in el) {
  daf_loop <- daf %>% 
    filter(!is.na(Latitude) & !is.na(Longitude) & !is.na(.data[[element]]))
  
  ptsel <- st_as_sf(daf_loop, coords = c("Longitude", "Latitude"), crs = 4326)
  el_lab <- gsub("_", "/", element)
  
  elmap <- ggplot() + 
    annotation_map_tile('cartolight') +
    geom_spatial_point(crs = 4326) +
    geom_sf(data = ptsel, aes(colour = .data[[element]]), show.legend = TRUE) +  
    scale_color_viridis_c(name = el_lab, trans = "log") +
    labs(x = "Longitude", y = "Latitude")
  
  elmaps[[element]] <- elmap
}

print(elmaps[[1]]) 
print(elmaps[[2]]) 
print(elmaps[[3]]) 
print(elmaps[[4]]) 
print(elmaps[[5]]) 
print(elmaps[[6]]) 
print(elmaps[[7]]) 
print(elmaps[[8]]) 
print(elmaps[[9]])
print(elmaps[[10]])
print(elmaps[[11]]) 

# and mitchell-----

elmaps <- list() 
for (element in el) {
  des_loop <- des %>% 
    filter(!is.na(Lat) & !is.na(Long) & !is.na(.data[[element]]))
  
  ptsel <- st_as_sf(des_loop, coords = c("Long", "Lat"), crs = 4326)
  el_lab <- gsub("_", "/", element)
  
  elmap <- ggplot() + 
    # annotation_map_tile('cartolight', zoom = 6) +
    geom_spatial_point(crs = 4326) +
    geom_sf(data = ptsel, aes(colour = .data[[element]]), show.legend = TRUE) +  
    scale_color_viridis_c(name = el_lab, trans = "log") +
    # coord_sf(lims_method = "geometry_bbox") +
    labs(x = "Longitude", y = "Latitude")
  
  elmaps[[element]] <- elmap
}

print(elmaps[[1]]) 
print(elmaps[[2]]) 
print(elmaps[[3]]) 
print(elmaps[[4]]) 
print(elmaps[[5]]) 
print(elmaps[[6]]) 
print(elmaps[[7]]) 
print(elmaps[[8]]) 
print(elmaps[[9]])
print(elmaps[[10]])
print(elmaps[[11]]) 
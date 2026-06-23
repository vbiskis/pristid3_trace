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

# combine data sources for comparability-----

daf <- daf %>% 
  mutate(source = 'daf', ID = CbottleID, r = 0, Month = match(Month, month.name),
         Species = 'water', riv = substr(ShortName, 1, 2), 
         Region = 'NEQ', Rsplit = 'NEQ', RivGrp = ifelse(Latitude <= -19, 'RG05', 'RG06'),
         EW = 'E', NS = 'N',
         ) %>% select(-AbottleID, -CbottleID, -ShortName)

des <- des %>% 
  mutate(source = 'des', ID = make.unique(as.character(Site_Number), sep = "."),
         Year = year(Date), Month = month(Date), r = 0, salty = "fresh",
         Latitude = Lat, Longitude = Long,
         Species = 'water', riv = 'Mitch', Region = 'GoC', Rsplit = 'NGoC', 
         RivGrp = 'RG10', EW = 'W', NS = 'N',
  ) %>% select(-Site_Number, -Site_Name, -Date, -Lat, -Long)

enviro <- bind_rows(daf, des)

## now get in the same units----

enviro1 <- enviro %>% 
  # create new data columns with PPM to Element:Ca conversion 
  mutate(Li_Ca = Li_Ca/6.94*40.08*1000000,
         B_Ca = B_Ca/10.81*40.08*1000000,
         Na_Ca = Na_Ca/22.99*40.08*1000000,
         Mg_Ca = Mg_Ca/24.31*40.08*1000000,
         Mn_Ca = Mn_Ca/54.94*40.08*1000000,
         Ni_Ca = Ni_Ca/58.69*40.08*1000000,
         Cu_Ca = Cu_Ca/63.55*40.08*1000000,
         Zn_Ca = Zn_Ca/65.38*40.08*1000000,
         Sr_Ca = Sr_Ca/87.62*40.08*1000000,
         Ba_Ca = Ba_Ca/137.33*40.08*1000000,
         Pb_Ca = Pb_Ca/207.2*40.08*1000000) # now in umol/mol 

enviro1 <- enviro1 %>% 
  mutate(across(all_of(el), ~log(.x + 1), .names = "ln{str_remove(.col, '_Ca')}"))

library(writexl)
write_xlsx(enviro1, path = "data/output/water_dat.xlsx")

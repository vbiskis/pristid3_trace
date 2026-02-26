# Details ----
#' mapHelp.R
#' 
#' Adapted from prev paper: [Putting sawfishes back on the map]
#' Author: Nikki Biskis
#' Date: 
#' - 18 Aug 2024
#' Content: General settings for mapping
#'  avoid repetition
#' -----------

pacman::p_load('dplyr', 'tidyverse', 
               'terra', 'sf', 'sp', 
               'rnaturalearth', 'maps', 'mapsf',
               'adehabitatHR',
               'tmap', 'tmaptools', 'devtools',
               'ggspatial', 'ggthemes', 
               'geodata', 'geosphere')

QLD_crs <- crs("EPSG:4326")
aus <- ne_countries(scale = "medium", country = "Australia", returnclass = "sf")%>%
  st_transform(QLD_crs)
studyarea <- ext(c(137, 154, -29, -10))
qldbox <- st_bbox(c(xmin = 138, xmax = 154, ymin = -28, ymax = -10), crs = QLD_crs)
grid0 <- st_make_grid(qldbox, n = c(1, 1)) #no grid
qld0 <- st_intersection(aus, grid0)

#qldboard <- vect("FetchEnviro/QLDSpatial/Coastline_and_State_border.shp")

# this thing i always do!
crd_to_xy <- function(coords, lat = "V2", lon = "V3", name = "pts.sp") { # make it so it can be def.
  sf_obj <- as.data.frame(coords) %>%
    st_as_sf(coords = c(lon, lat), 
             crs = st_crs("+proj=longlat +datum=WGS84")) # could swap?
  assign(name, sf_obj, envir = .GlobalEnv)
}

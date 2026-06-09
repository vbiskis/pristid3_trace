# Details ----
#' 02_study_area.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 8 Aug 2024
#' - Updated: 20 Feb 2026
#' 
#' Content: 
#' + mapping study regions
#' + pulled out mapping libraries / formatting
#' ++ made source code, style supplied in microHelp.R
#' -----------

source('helpers/mapHelp.R')
source('helpers/microHelp.R')
set_theme(mytheme)
microE <- read_xlsx('data/output/microE.xlsx')
library(writexl)
library(RColorBrewer)

# LODs input in place of neg vals
# high In samples excluded, and lowest value chosen per site.
# each point reps one tooth

# Species Points----
saw_coords <- cbind(microE$Species, microE$LAT, microE$LONG)
crd_to_xy(saw_coords)
pts.sp$V1 <- as.factor(pts.sp$V1)

specmap <- ggplot() + 
  annotation_map_tile('cartolight') +
  geom_spatial_point(crs = 4326) +
  geom_sf(data = pts.sp, aes(colour = V1), show.legend = TRUE) +  
  scale_color_manual(name = "Species", 
                     values = spec_col, labels = spec_labs) +  
  labs(x = "Longitude", y = "Latitude")
specmap

## Region----
microEF <- microE %>% 
  filter(RivGrp != "RG04", RivGrp != "RG11", RivGrp != "RG13")

coords_ef <- cbind(microEF$Species, 
                    microEF$LAT, microEF$LONG,
                    microEF$RivGrp, microEF$Rsplit, microEF$Region)

crd_to_xy(coords_ef, name = "ptsef")
ptsef <- ptsef %>% 
  rename(Species = V1, RivGrp = V4, Rsplit = V5, Region = V6)

#as points
sitemap <- ggplot() + 
  annotation_map_tile('cartolight') +
  geom_spatial_point(crs = 4326) +
  geom_sf(data = ptsef, aes(colour = RivGrp), show.legend = TRUE) +  
  labs(x = "Longitude", y = "Latitude") 
sitemap

# nice!
# now lets do this with boxes around the centroid instead

# fig 1----
pts_centroids <- ptsef %>%
  group_by(RivGrp) %>%
  summarise(geometry = st_centroid(st_union(geometry)), .groups = "drop")

n_groups <- length(unique(microEF$RivGrp[!microEF$RivGrp %in% c("RG04", "RG11", "RG13")]))
rdylgn_colors <- RColorBrewer::brewer.pal(9, "RdYlGn") # mmmm k I like most of those
# want to remove that yellow one tho
map_cols <- c("#A50026", #WCY
              "#D73027", #NG
              "#F46D43", #SG
              "#FDAE61", "#FEE08B", #FNQ
              "#D9EF8B", "#A6D96A", #NQ
              "#66BD63", "#006837", #CEQ
              "#6495ED") #SEQ

sitemap <- ggplot() +
  geom_sf(data = qld0, fill = "gray90", color = "black") +
  geom_sf(data = pts_centroids, aes(color = RivGrp), size = 3) +
  labs(x = "Longitude", y = "Latitude", color = "River Groups") +
  scale_color_brewer(palette = "RdYlGn", name = "River Groups") +
  scale_color_manual(values = rev(map_cols)) +
  annotation_custom(
    ggplotGrob(
      ggplot() +
        geom_sf(data = aus, fill = "gray85", color = "black") + 
        theme_map() +
        theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
    ),
    xmin = 148, xmax = 154, ymin = -16, ymax = -10  # Position of inset
  ) +
  geom_rect(aes(xmin = 148.1, xmax = 153.6, ymin = -14.6, ymax = -10.2), 
            color = "black", fill = "transparent", size = .5) +
  theme_bw() +
  theme(
    legend.position = "none"
  )

ggsave(
  "fig1a_study_area.png",
  plot = sitemap,
  device = NULL,
  path = 'figures/fig1/',
  scale = 1,
  width = 5,
  height = 6,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
)

# need bar chart!
map_cols_w <- c("#A50026","white", "#D73027", "#F46D43", "#FDAE61", "#FEE08B", "#D9EF8B", "#A6D96A", "#66BD63", "#006837", "#6495ED") 
map_cols_o <- c("#A50026","#D73027", "#D73027", "#F46D43", "#FDAE61", "#FEE08B", "#D9EF8B", "#A6D96A", "#66BD63", "#006837", "#6495ED") 

ggplot() +
  geom_bar(data = microE %>% 
             filter(RivGrp != "RG04", RivGrp != "RG13"), 
           aes(x = Species, fill = RivGrp, color = RivGrp),
           alpha = 0.9,
           width = 0.75) +
  scale_x_discrete(labels = c("A. cuspidata", "P. pristis", "P. zijsron")) +
  scale_fill_manual(name = "River Groups",
                    values = rev(map_cols_w)) +
  scale_color_manual(name = "River Groups",
                     values = rev(map_cols_o)) +
  labs(y = "Number of Rostra") +
  theme_bw() +
  theme(
    axis.text.x = element_text(face = "italic", size = 11),
    axis.text.y = element_text(size = 11), axis.title = element_text(size = 12),
    legend.title = element_text(size = 12), legend.text = element_text(size = 11)
  ) 

ggsave(
  "fig1b_groups.png",
  plot = last_plot(),
  device = NULL,
  path = 'figures/fig1/',
  scale = 1,
  width = 4.5,
  height = 6,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
)

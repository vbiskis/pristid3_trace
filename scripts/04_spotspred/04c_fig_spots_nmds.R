# Details ----
#' 04c_spots_nmds.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 7 Jul 2025
#' 
#' Content: nmds visualisation
#' Update 24 Feb 2026: theme update w/ microHelp.R
#' -----------

library(vegan)
source('helpers/microHelp.R')
set_theme(mytheme)

microE <- read_xlsx('data/output/microE.xlsx')
microE <- microE %>% 
  mutate(across(where(is.character) & !all_of(c("ID", "File_ID")), as.factor))

microE$Rsplit <- factor(microE$Rsplit, 
                         levels = c("AQ", "WCY", "NGoC", "SGoC", "UGoC",
                                    "FNEQ", "NEQ", "CEQ", "SEQ")) # right order

#AllAges----
dist_matrix <- vegdist(microE[, c('lnLi', 'lnZn', 'lnPb', 'lnB',
                                  #'lnNi', 'lnCu',  #bad
                                  'lnMn', 'lnMg', 'lnSr', 'lnBa')], method = "bray", 
                       na.rm = TRUE)

#Side Plot!----
disp <- betadisper(dist_matrix, microE$Species) # wait i love him!
plot(disp)

##nMDS----
MD <- metaMDS(dist_matrix, k = 3)
print(paste("Stress:", round(MD$stress, 3)))
plot(MD, type = "n")
points(MD, col = as.factor(microE$Species), 
       pch = as.numeric(as.factor(microE$Region)))

##Ordiplot----
md_scores <- scores(MD)
md_df <- data.frame(
  NMDS1 = md_scores[,1],
  NMDS2 = md_scores[,2],
  Species = microE$Species,
  Region = microE$Rsplit
)

nma <- ggplot(md_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Species, shape = Region), size = 3) +
  stat_ellipse(aes(color = Species), level = 0.95) +  # Add confidence ellipses
  scale_shape_manual(values = c(9, 18, 2, 25, 4, 22, 15, 1, 19)) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  labs(title = paste("NMDS Ordination (Stress =", round(MD$stress, 3), ")")) +
  scale_color_manual(values = spec_col, labels = spec_labs)

ggsave(
  "nmds_a.png",
  plot = nma,
  device = NULL,
  path = 'figures/fig4/',
  scale = 1,
  width = 6,
  height = 5,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  bg = "transparent",
  limitsize = TRUE
)

#Juvenile Only----
microjuv <- microE %>% filter(Class != "ADULT", Class != "BTWN") 

dist_matrix <- vegdist(microjuv[, c('lnLi', 'lnZn', 'lnPb', 'lnB',
                                    #'lnNi', 'lnCu', #bad
                                    'lnMn', 'lnMg', 'lnSr', 'lnBa')], 
                       method = "bray",  na.rm = TRUE)

## nMDS----
MD <- metaMDS(dist_matrix, k = 3)
print(paste("Stress:", round(MD$stress, 3)))
plot(MD, type = "n")
points(MD, col = as.factor(microjuv$Species), 
       pch = as.numeric(as.factor(microjuv$Region)))

##Ordiplot----
md_scores <- scores(MD)
md_df2 <- data.frame(
  NMDS1 = md_scores[,1],
  NMDS2 = md_scores[,2],
  Species = microjuv$Species,
  Region = microjuv$Rsplit
)

nmb <- ggplot(md_df2, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Species, shape = Region), size = 3) +
  stat_ellipse(aes(color = Species), level = 0.95) +  # Add confidence ellipses
  scale_shape_manual(values = c(9, 2, 25, 4, 15, 1, 19)) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  labs(title = paste("NMDS Ordination (Stress =", round(MD$stress, 3), ")")) +
  scale_color_manual(values = spec_col, labels = spec_labs)

ggsave(
  "nmds_b.png",
  plot = nmb,
  bg = "transparent",
  path = 'figures/fig4/',
  scale = 1,
  width = 7,
  height = 5,
  dpi = 300)

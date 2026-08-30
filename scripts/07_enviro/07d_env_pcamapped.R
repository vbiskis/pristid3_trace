# Details ----
#' 07d_env_pcamapped.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 22 Jun 2026
#' 
#' Content: 
#' Fig p3 & p4: water data mapped to tooth data in same pcs axes
#' -----------

source('helpers/microHelp.R')
set_theme(mytheme)
microE <- read_xlsx('data/output/microE.xlsx')
waterdat <- read_xlsx('data/output/water_dat.xlsx')
pca_all <- readRDS(file = 'data/output/pca_spots_base.rds')

# project water samples into existing teeth PCA space----
water_scores <- predict(pca_all, newdata = waterdat %>% 
                          dplyr::select(all_of(lnrats)) %>% 
                          as.matrix())

water_scores_df <- data.frame(water_scores, type = "water",
                              RivGrp = waterdat$RivGrp, 
                              salty = factor(waterdat$salty,
                                             levels = c("fresh", "brackish", "sea")))
glimpse(water_scores_df)

## species graph----
PCA_all <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(data = data.frame(pca_all$x, microE),
             aes(x = PC1, y = PC2, 
                 color = Species, shape = RivGrp)) +
  scale_shape_manual(values = 1:13) + 
  scale_color_manual(values = spec_col, labels = spec_labs) +
  stat_ellipse(data = data.frame(pca_all$x, microE),
               aes(x = PC1, y = PC2, color = Species, group = Species),
               level = 0.95) + 
  labs(x = "PC1", y = "PC2") 

## Overlay onto existing plot----
library(ggrepel)
library(ggnewscale)

loadings <- pca_all$rotation  
loadings_df <- as.data.frame(loadings)
loadings_df$Element <- rownames(loadings_df)

scale_factor <- max(abs(pca_all$x[,1:2]))/ max(abs(loadings[,1:2]))
loadings_scaled <- as.data.frame(loadings*scale_factor*0.75)
loadings_scaled$Element <- rownames(loadings_scaled)

p3 <- PCA_all +
  geom_segment(data = loadings_scaled,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.3, "cm")), color = "grey10") +
  geom_text_repel(data = loadings_scaled, aes(x = PC1, y = PC2, label = Element),
                  size = 4, color = "grey10", family = "optima") +
  new_scale_color() +
  geom_point(data = water_scores_df,
             aes(x = PC1, y = PC2, shape = RivGrp, color = salty)) +
  scale_color_manual(values = c("fresh" = "steelblue", "brackish" = "purple4", "sea" = "orange3"),
                     name = "Water") +
  geom_segment(data = loadings_scaled,
               aes(x = -10, y = 10, xend = 2*PC1-10, yend = 2*PC2+10),
               arrow = arrow(length = unit(0.3, "cm")), color = "grey10") +
  geom_text_repel(data = loadings_scaled, aes(x = 2*PC1-10, y = 2*PC2+10, label = Element),
                  size = 4, color = "grey10", family = "optima") +
  guides(shape = guide_legend(order = 3, ncol = 2),
         color = guide_legend(order = 2, ncol = 1),
         color_new = guide_legend(order = 1, ncol = 1)) +
  theme(legend.box = "vertical")

ggsave(
  "pca_tog.png",
  plot = p3,
  bg = "transparent",
  path = 'figures/grantrep/',
  scale = 1,
  width = 8,
  height = 6,
  dpi = 300)

### or facet----
teeth_scores_df <- data.frame(pca_all$x, microE) %>%
  mutate(type = "teeth")

# bind cols for plotting
combined_df <- bind_rows(
  teeth_scores_df %>% dplyr::select(PC1, PC2, RivGrp, Species, type),
  water_scores_df %>% mutate(Species = NA_character_) %>% 
    dplyr::select(PC1, PC2, RivGrp, type, salty, Species)
) %>% mutate(type = factor(type, levels = c("water", "teeth")))

loadings_scaled_teeth <- loadings_scaled %>% mutate(type = "teeth",
                                                    xs = 0, ys = 0)
loadings_scaled_water <- loadings_scaled %>% 
  mutate(type = "water", PC1 = PC1*2.1, PC2 = PC2*2.1,
         xs = -2.5, ys = 8.5) # for the arrows

loadings_both <- bind_rows(loadings_scaled_teeth, loadings_scaled_water)

p4 <- ggplot(combined_df, aes(x = PC1, y = PC2, color = Species)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(data = . %>% filter(type == "teeth"),
             aes(shape = RivGrp, color = Species)) +
  stat_ellipse(data = . %>% filter(type == "teeth"), level = 0.95) +
  scale_color_manual(values = spec_col, labels = spec_labs, na.translate = FALSE) +
  new_scale_color() +
  geom_point(data = . %>% filter(type == "water"),
             aes(x = PC1, y = PC2, shape = RivGrp, color = salty), inherit.aes = FALSE) +
  geom_segment(data = loadings_both,
               aes(x = xs, y = ys, xend = PC1 + xs, yend = PC2 + ys),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "grey10", inherit.aes = FALSE) +
  geom_text_repel(data = loadings_both, 
                  aes(x = PC1 + xs, y = PC2 + ys, label = Element),
                  size = 4, color = "grey10", 
                  family = "optima", inherit.aes = FALSE) +
  guides(shape = guide_legend(order = 1, ncol = 2),
         color = guide_legend(order = 2, ncol = 1)) +
  scale_color_manual(values = c("fresh" = "steelblue",  "brackish" = "purple4", 
                                "sea" = "orange3"),
                     labels = c("Fresh", "Brackish", "Salty"),
                     name = "Water Salinity") +
  facet_wrap(~ factor(type, levels = c("water", "teeth")), scales = "free") +
  scale_shape_manual(values = 1:13) + 
  labs(x = "PC1", y = "PC2", shape = "River Group") +
  theme(legend.box = "vertical", strip.text = element_blank())

ggsave(
  "pca_facet.png",
  plot = p4,
  bg = "transparent",
  path = 'figures/grantrep/',
  scale = 1,
  width = 9,
  height = 6,
  dpi = 300)

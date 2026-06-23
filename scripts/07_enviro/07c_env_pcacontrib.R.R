# Details ----
#' 07c_env_pcacontrib.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 23 Jun 2026
#' 
#' Content: 
#' Fig p2: water data pca w/ region
#' Contribution table
#' Fig p3/p4: mapped to tooth data in 7d
#' -----------

source('helpers/microHelp.R')
set_theme(mytheme)
microE <- read_xlsx('data/output/microE.xlsx')
waterdat <- read_xlsx('data/output/water_dat.xlsx')
pca_all <- readRDS(file = 'data/output/pca_spots_base.rds')

# just looking at the water samples----

PCA2 <- waterdat %>%
  dplyr::select(all_of(lnrats)) %>% 
  na.omit() %>%
  as.matrix()

pca_wat <- prcomp(PCA2, scale. = TRUE)
summary(pca_wat) # 65% explained in 2 axes

loadingsw <- pca_wat$rotation  
loadingsw_df <- as.data.frame(loadingsw)
loadingsw_df$Element <- rownames(loadingsw_df)

scale_factorw <- max(abs(pca_wat$x[,1:2]))/ max(abs(loadingsw[,1:2]))
loadingsw_scaled <- as.data.frame(loadingsw*scale_factorw*0.75)
loadingsw_scaled$Element <- rownames(loadingsw_scaled)

PCA_wat <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(data = data.frame(pca_wat$x, waterdat),
             aes(x = PC1, y = PC2, shape = RivGrp, color = salty), size = 3) +
  scale_color_manual(values = c("purple4", "steelblue", "orange3"),
                     labels = c("Brackish", "Fresh", "Salty")) + 
  scale_shape_manual(values = c(5,6,10)) + 
  stat_ellipse(data = data.frame(pca_wat$x, waterdat),
               aes(x = PC1, y = PC2, group = RivGrp),
               level = 0.95) + 
  labs(x = "PC1", y = "PC2", shape = "River Group",
       color = "Conductivity") 

p2 <- PCA_wat +
  geom_segment(data = loadingsw_scaled,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.3, "cm")), color = "grey10") +
  geom_text(data = loadingsw_scaled, aes(x = PC1 - 1.25, y = PC2 + 0.25, label = Element),
            hjust = -0.2, vjust = 0.2, size = 4, 
            color = "grey10", family = "optima")

ggsave(
  "pca_wat.png",
  plot = p2,
  bg = "transparent",
  path = 'figures/grantrep/',
  scale = 1,
  width = 8,
  height = 6,
  dpi = 300)

# contribution----

library(factoextra)
fviz_eig(pca_wat)
fviz_pca_var(pca_wat, col.var = "cos2",
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) # ah remove those then

get_pca_var(pca_wat)$cos2

wat_var <- get_pca_var(pca_wat)
# str(wat_var)
water_contrib <- as.data.frame(wat_var$contrib) %>%
  rownames_to_column("Element") %>%
  dplyr::select(Element, Dim.1, Dim.2) %>%
  rename(W_contrib_PC1 = Dim.1, W_contrib_PC2 = Dim.2)

water_cor <- as.data.frame(wat_var$cor) %>%
  rownames_to_column("Element") %>%
  dplyr::select(Element, Dim.1, Dim.2) %>%
  rename(W_cor_PC1 = Dim.1, W_cor_PC2 = Dim.2)

## compare to teeth----
toof_var <- get_pca_var(pca_all)
teeth_contrib <- as.data.frame(toof_var$contrib) %>%
  rownames_to_column("Element") %>%
  dplyr::select(Element, Dim.1, Dim.2) %>%
  rename(T_contrib_PC1 = Dim.1, T_contrib_PC2 = Dim.2)

teeth_cor <- as.data.frame(toof_var$cor) %>%
  rownames_to_column("Element") %>%
  dplyr::select(Element, Dim.1, Dim.2) %>%
  rename(T_cor_PC1 = Dim.1, T_cor_PC2 = Dim.2)

## join----
contrib_table <- left_join(teeth_contrib, teeth_cor, by = "Element") %>%
  left_join(water_contrib, by = "Element") %>%
  left_join(water_cor, by = "Element") %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))

contrib_table #meh

# compare----
# sum cos2 across Dim1 and Dim2 for each dataset
cos2_compare <- data.frame(
  Element = rownames(get_pca_var(pca_all)$cos2),
  Teeth = get_pca_var(pca_all)$cos2[,1] + get_pca_var(pca_all)$cos2[,2],
  Water = get_pca_var(pca_wat)$cos2[,1] + get_pca_var(pca_wat)$cos2[,2]
)

cosplot <- ggplot(cos2_compare, aes(x = Teeth, y = Water, label = Element)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  geom_point(size = 3, color = "steelblue") +
  geom_text_repel(family = "optima", size = 3.5) +
  labs(x = "Teeth cos² (PC1+PC2)", 
       y = "Water cos² (PC1+PC2)") +
  xlim(0, 1) + ylim(0, 1)

ggsave(
  "cosplot.png",
  plot = cosplot,
  bg = "transparent",
  path = 'figures/grantrep/',
  scale = 1,
  width = 5,
  height = 5,
  dpi = 300)

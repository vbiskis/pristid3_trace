# Details ----
#' 03f_pca_cap.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 18 Aug 2024
#' 
#' Content: original pca w/ region
#' includes which elements are irrelevant to predictive success
#' + Update 24 Feb 2026
#' ++ output cleaned up with theme in microHelp.R
#' -----------

source('helpers/microHelp.R')
set_theme(mytheme)
microE <- read_xlsx('data/output/microE.xlsx')

# capture loc----

#all base
microE <- microE %>%
  na.omit()
  
PCA1 <- microE %>%
  dplyr::select(all_of(lnrats)) %>% 
  na.omit() %>%
  as.matrix()

pca_all <- prcomp(PCA1, scale. = TRUE)
summary(pca_all) # 53% explained in 2 axes

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

loadings <- pca_all$rotation  
loadings_df <- as.data.frame(loadings)
loadings_df$Element <- rownames(loadings_df)

scale_factor <- max(abs(pca_all$x[,1:2]))/ max(abs(loadings[,1:2]))
loadings_scaled <- as.data.frame(loadings*scale_factor*0.75)
loadings_scaled$Element <- rownames(loadings_scaled)

PCA_all +
  geom_segment(data = loadings_scaled,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.3, "cm")), color = "grey10") +
  geom_text(data = loadings_scaled, aes(x = PC1, y = PC2, label = Element),
            hjust = -0.2, vjust = 0.2, size = 4, 
            color = "grey10", family = "optima")

library(factoextra)
fviz_eig(pca_all)
fviz_pca_var(pca_all, 
             col.var = "cos2",
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) # ah remove those then

contrib <- get_pca_var(pca_all)

# Details ----
#' 03f_pca_spec.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' Presented at: [IPFC 2024] [OCS 2025]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date:  18 Aug 2024
#' 
#' Content: pcas/umaps sep by species (for presentations)
#' + Update 24 Feb 2026:
#' ++ output/calls cleaned up with microHelp.R
#' -----------

source('helpers/microHelp.R')
library(factoextra)
set_theme(mytheme)
microE <- read_xlsx('data/output/microE.xlsx')

# Filter Species ----
AC <- microE %>%
  filter(Species == 'Anoxypristis cuspidata') %>% 
  mutate(Maturity = case_when(Age_SP >= 2 ~ "Mature",
                              .default = "Immature"))
AC <- microE %>%
  filter(Species == 'Pristis pristis')

AC <- microE %>%
  filter(Species == 'Pristis zijsron') 

#PCA----
## anoxy----
AC <- AC %>% 
  filter(Maturity == "Immature")

goodrats <- setdiff(lnrats, "lnPb") # from prev run

keep_all <- AC %>% filter(r <= 101)
sample_20 <- AC %>% filter(r > 101) %>% slice_sample(n = 20)

AC <- bind_rows(keep_all, sample_20)

PCA_A <- AC %>%
  dplyr::select(all_of(goodrats)) %>% 
  as.matrix()

pca_a <- prcomp(PCA_A, scale. = TRUE)
summary(pca_a)

PCA_anoxy <- ggplot() + # can adjust to be maturity, etc
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(data = data.frame(pca_a$x, AC),
             aes(x = PC1, y = PC2, 
                 shape = Region, size = Age_SP)) +
  scale_shape_manual(values = 1:10) + 
  stat_ellipse(data = data.frame(pca_a$x, AC),
               aes(x = PC1, y = PC2, 
                   group = Region, color = Region),
               level = 0.9) + 
  labs(x = "PC1", y = "PC2", size = 'Age') #separating on age!

fviz_eig(pca_a, 
         choice = "eigenvalue",  # or "variance" 
         addlabels = TRUE, #main = "Scree Plot",
         xlab = "Principal Components", ylab = "Eigenvalue")
fviz_contrib(pca_a, choice = "var", axes = 1:2)
fviz_pca_var(pca_a, 
             col.var = "cos2",
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) # Zn for anoxy maturity!

ggsave(
  "ac_pca.png",
  plot = PCA_anoxy,
  device = NULL,
  path = 'figures/ext/',
  scale = 1,
  width = 6,
  height = 6,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
)

### regions? ----
elemA <- AC %>%
  dplyr::select(all_of(els))

# Standardize the elemental signatures
scaled_data <- scale(elemA)
umap_result <- umap(scaled_data)
umap_A <- data.frame(umap_result$layout)

set.seed(123)  
kmeans_result <- kmeans(scaled_data, centers = 4)
AC$Cluster <- as.factor(kmeans_result$cluster)
umap_A <- cbind(umap_A, AC)

UMAP_anoxy <- ggplot(umap_A) +
  geom_point(aes(x = X1, y = X2, color = Year, # swap out maturity?
                 shape = RivGrp)) +
  scale_shape_manual(values = 1:10) + 
  stat_ellipse(aes(x = X1, y = X2, 
                   group = Cluster),
               level = 0.95) + 
  labs(x = "UMAP1", y = "UMAP2") 
# should just do with juveniles - i.e. cap*

ggsave(
  "ac_umap.png",
  plot = UMAP_anoxy,
  device = NULL,
  path = 'figures/ext/',
  scale = 1,
  width = 6,
  height = 6,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
)

## freshie----
PCA_P <- PP %>%
  dplyr::select(all_of(els)) %>% 
  as.matrix()

pca_p <- prcomp(PCA_P, scale. = TRUE)
summary(pca_p)

ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(data = data.frame(pca_p$x, ACP),
             aes(x = PC1, y = PC2, 
                 color = Year, shape = RivGrp)) +
  stat_ellipse(data = data.frame(pca_p$x, ACP),
               aes(x = PC1, y = PC2, 
                   group = Rsplit),
               level = 0.95) + 
  labs(x = "PC1", y = "PC2") # checks out

elemP <- PP %>%
  dplyr::select(all_of(els))

scaled_data <- scale(elemP)
umap_result <- umap(scaled_data)
umap_P <- data.frame(umap_result$layout)

set.seed(123)  
kmeans_result <- kmeans(scaled_data, centers = 2)
PP$Cluster <- as.factor(kmeans_result$cluster)
umap_P <- cbind(umap_P, PP)

UMAP_pristis <- ggplot(umap_P) +
  geom_point(aes(x = X1, y = X2, color = Year, shape = RivGrp)) +
  scale_shape_manual(values = 1:10) + 
  stat_ellipse(aes(x = X1, y = X2, 
                   group = Cluster),
               level = 0.95) + 
  labs(x = "UMAP1", y = "UMAP2") 

ggsave(
  "pp_umap.png",
  plot = UMAP_pristis,
  device = NULL,
  path = 'figures/ext/',
  scale = 1,
  width = 6,
  height = 4,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
)

## green----
PCA_Z <- PZ %>%
  dplyr::select(all_of(els)) %>% 
  as.matrix()

pca_z <- prcomp(PCA_Z, scale. = TRUE)

ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(data = data.frame(pca_z$x, PZ),
             aes(x = PC1, y = PC2, 
                 color = Year, shape = Region)) +
  stat_ellipse(data = data.frame(pca_z$x, PZ),
               aes(x = PC1, y = PC2, 
                   group = Region),
               level = 0.95) + 
  labs(x = "PC1", y = "PC2") #eek!

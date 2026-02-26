# Details ----
#' 01b_batch_check.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 
#' - Origin: 12 Sep 2024
#' - Updated: 20 Feb 2026
#' 
#' Content: Quality Control 2
#' + make sure results aren't influenced by batch (random)
#' 
#' Note:
#' + xlsx prepped for analysis in 00_prep_spots
#' + style supplied in microHelp.R
#' -----------

source('MicroHelp.R')
set_theme(mytheme)
micro <- readxl::read_xlsx("data/output/icpms_clean.xlsx") 

#PCA----
Batch <- micro %>%
  dplyr::select(all_of(ln_rats)) %>% 
  as.matrix() %>% 
  na.omit()
  
pca_b <- prcomp(Batch, scale. = TRUE)
summary(pca_b)

##Biplot----
biplot(pca_b)

micro_all <- micro %>%
  filter(complete.cases(dplyr::select(., all_of(el))))

#lets color it by rand effects
ggplot() +
  facet_wrap(~ Species, nrow = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(data = data.frame(pca_b$x, micro_all),
             aes(x = PC1, y = PC2, color = Batch, shape = Rsplit)) +
  stat_ellipse(data = data.frame(pca_b$x, micro_all),
               aes(x = PC1, y = PC2, group = Batch, color = Batch),
               level = 0.95, type = "norm") +
  labs(x = "PC1", y = "PC2") 

# can see that I put related P. pristis on same slides, due to size
# otherwise its fine
# Details ----
#' 04d_spots_anosim.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 12 Sep 2024
#' 
#' Content: 
#' permanova etc. for spec diffs (mainly)
#' + Update 25 Feb 2026
#' ++ cleaned up for functional independence, otherwise same
#' -----------

source('helpers/microHelp.R')
library(vegan)
microE <- read_xlsx('data/output/microE.xlsx')
microE <- microE %>% 
  mutate(across(where(is.character) & !all_of(c("ID", "File_ID")), as.factor),
         Rsplit = factor(Rsplit, levels = c("AQ", "WCY", "NGoC", "SGoC", "UGoC",
                                            "FNEQ", "NEQ", "CEQ", "SEQ")),
         Region = factor(Region, levels = c('AQ', "WCY", "GoC", "NEQ", "CEQ", "SEQ")),
         Class = factor(Class, levels = c('YOY', 'IMMATURE', 'BTWN', 'ADULT')))

sigrats <- c('lnLi', 'lnZn', 'lnPb', 'lnB', 'lnMn',
             #'lnNi', #'lnCu', # adjust if necessary, not sig in the pca
             'lnMg', 'lnSr', 'lnBa')

dist_matrix <- vegdist(microE[, sigrats], method = "bray", na.rm = TRUE)

# CAP----
cap_resAll <- capscale(dist_matrix ~ Region + Species, data = microE)
summary(cap_resAll)
plot(cap_resAll, type = "n")  # Empty plot
points(cap_resAll, display = "sites",
       col = spec_col[as.numeric(as.factor(microE$Species))],
       pch = as.numeric(as.factor(microE$Maturity)))
text(cap_resAll, display = "bp", col = "grey50")  # Biplot arrows
legend("topright", legend = levels(as.factor(microE$Species)), 
       col = spec_col, pch = 19)

anova(cap_resAll)  
anova(cap_resAll, by = "terms")
anova(cap_resAll, by = "axis")  

cap_resAll$CCA$tot.chi / cap_resAll$tot.chi  # Proportion constrained variance
scores(cap_resAll, display = "bp") 

#PERMANOVA
permanova <- adonis2(dist_matrix ~ Region + Species, data = microE)
permanova

#SIMPER
simper_spec <- simper(microE[, lnrats], group = microE$Species)
summary(simper_spec) # between spec groups

simper_reg <- simper(microE[, lnrats], group = microE$Region)
summary(simper_reg) # between regions

# ANOSIM test
# species----
dist_mat2 <- dist(microE[, lnrats], method = "euclidean")
ano_res <- anosim(dist_mat2, grouping = microE$Species, permutations = 999)
print(ano_res)

# Pairwise ANOSIM between all species
species_list <- unique(microE$Species)
n_species <- length(species_list) # 3 obvs

# Create matrices for R-values and p-values
r_matrix <- matrix(NA, n_species, n_species)
p_matrix <- matrix(NA, n_species, n_species)
rownames(r_matrix) <- colnames(r_matrix) <- species_list
rownames(p_matrix) <- colnames(p_matrix) <- species_list

# Fill diagonal with dashes
diag(r_matrix) <- "–"
diag(p_matrix) <- "–"

# Calculate pairwise comparisons
for(i in 1:(n_species-1)) {
  for(j in (i+1):n_species) {
    sp1 <- species_list[i]
    sp2 <- species_list[j]
    
    # Subset data for these two species
    subset_indices <- which(microE$Species %in% c(sp1, sp2))
    subset_dist <- as.dist(as.matrix(dist_matrix)[subset_indices, subset_indices])
    subset_groups <- microE$Species[subset_indices]
    
    # Run ANOSIM
    anosim_result <- anosim(subset_dist, grouping = subset_groups, permutations = 999)
    
    # Store R-value (upper triangle) and p-value (lower triangle)
    r_matrix[i, j] <- round(anosim_result$statistic, 3)
    p_matrix[j, i] <- round(anosim_result$signif, 4)
  }
}

# Combine into final table
final_matrix <- r_matrix
final_matrix[lower.tri(final_matrix)] <- p_matrix[lower.tri(p_matrix)]

print(final_matrix) # that makes sense!

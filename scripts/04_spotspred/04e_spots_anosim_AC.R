# Details ----
#' 04e_spots_anosim_AC.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 7 Jul 2025
#' 
#' Content: 
#' permanova etc: focus on region diffs, just A. cuspidata
#' + Update 25 Feb 2026
#' ++ removed redundant tests
#' -----------

source('helpers/microHelp.R')
library(vegan)
microE <- read_xlsx('data/output/microE.xlsx')
microE <- microE %>% 
  mutate(across(where(is.character) & !all_of(c("ID", "File_ID")), as.factor),
         Rsplit = factor(Rsplit, levels = c("AQ", "WCY", "NGoC", "SGoC", "UGoC",
                                            "FNEQ", "NEQ", "CEQ", "SEQ")),
         Region = factor(Region, levels = c('AQ', "WCY", "GoC", "NEQ", "CEQ", "SEQ")))

AC <- microE %>%
  mutate(Region2 = case_when(Region == "WCY" ~ "GoC",
                             .default = Region),
         Rsplit2 = case_when(Region == "WCY" | Region == "GoC" ~ "GoC",
                             .default = Rsplit),
         Maturity = case_when(Age_SP >= 2 ~ "Mature", .default = "Immature")) %>% 
  filter(Species == 'Anoxypristis cuspidata') 

ACJ <- AC %>%
  filter(Maturity == 'Immature') 
ACJ$RivGrp <- droplevels(ACJ$RivGrp)

# Go 1: All ages----
elsA <- AC[, lnrats]
distA <- vegdist(elsA, method = "bray")

#PERMANOVA
permA <- adonis2(distA ~ Maturity + Rsplit2, data = AC)

#SIMPER
simperAC <- simper(elsA, group = AC$Rsplit2)
summary(simperAC) # between groups

simperAC2 <- simper(elsA, group = AC$Maturity)
summary(simperAC2) # between groups

# ANOSIM test
distA2 <- dist(elsA, method = "euclidean")
anosimA <- anosim(distA2, grouping = AC$Rsplit2, permutations = 999)
print(anosimA) # k its sig too

#Output
species_list <- unique(AC$Region2)
n_species <- length(species_list)

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
    subset_indices <- which(AC$Region2 %in% c(sp1, sp2))
    subset_dist <- as.dist(as.matrix(distA)[subset_indices, subset_indices])
    subset_groups <- AC$Region2[subset_indices]
    
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

print(final_matrix) #ooh looks like it!

# do have those adults in ther tho
# so lets just pull juv

## 2: Juveniles----
elsJ <- ACJ[, lnrats]
distJ <- vegdist(elsJ, method = "bray")

#PERMANOVA
permJ <- adonis2(distJ ~ Region2, data = ACJ)

#SIMPER
simperACJ <- simper(elsJ, group = ACJ$Region2)
summary(simperACJ) # between groups

# ANOSIM test
distAJ2 <- dist(elsJ, method = "euclidean")
anosimJ <- anosim(distAJ2, grouping = ACJ$Region2, permutations = 999)
print(anosimJ) # still sig

# dont have enough replicates within river groups for only juveniles


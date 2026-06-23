# Details ----
#' 07e_env_tests.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 23 Jun 2026
#' 
#' Content: 
#' permanova, beta dispersion (for a check), and partial cap comparison
#' -----------

source('helpers/microHelp.R')
set_theme(mytheme)
microE <- read_xlsx('data/output/microE.xlsx')
waterdat <- read_xlsx('data/output/water_dat.xlsx')

waterdat <- waterdat %>% 
  mutate(Season = if_else(Month <= 5, 'wet', 'dry'))

library(vegan)

# permanova-----
# water
water_mat <- waterdat %>% 
  dplyr::select(all_of(lnrats)) %>% 
  na.omit()

adonis2(water_mat ~ Region, data = waterdat, 
        method = "euclidean", permutations = 999) # looks like ti but

adonis2(water_mat ~ Region + as.factor(Year) + as.factor(Season) + salty, 
        data = waterdat,
        method = "euclidean",
        by = "margin",
        permutations = 999)

# check dispersion----

dist_water <- dist(water_mat, method = "euclidean")
disp_water <- betadisper(dist_water, waterdat$Region) 
# dispersed if using RivGrp bc of the low n iin RG06
permutest(disp_water)
plot(disp_water)
boxplot(disp_water)

# cap----

# constrain on region
cap_water_partial <- capscale(water_mat ~ Region + Condition(as.factor(Year) + Season),
                              data = waterdat,
                              distance = "euclidean")
unique(waterdat$Year)

cap_teeth_partial <- capscale(teeth_mat_sub ~ Region + Condition(as.factor(Year), Species),
                              data = teeth_sub,
                              distance = "euclidean")
unique(teeth_sub$Year) # awesome so similar

# extract CAP1 scores from partial CAPs
water_cap1_partial <- scores(cap_water_partial, display = "species") %>%
  as.data.frame() %>%
  rownames_to_column("Element") %>%
  dplyr::select(Element, CAP1) %>%
  rename(Water_CAP1 = CAP1)

teeth_cap1_partial <- scores(cap_teeth_partial, display = "species") %>%
  as.data.frame() %>%
  rownames_to_column("Element") %>%
  dplyr::select(Element, CAP1) %>%
  rename(Teeth_CAP1 = CAP1)

# merge
cap_compare_partial <- left_join(water_cap1_partial, teeth_cap1_partial, by = "Element")

# to long
cap_compare_long <- cap_compare_partial %>%
  pivot_longer(-Element, names_to = "Dataset", values_to = "CAP1") %>%
  mutate(Dataset = factor(Dataset, levels = c("Teeth_CAP1", "Water_CAP1"),
                          labels = c("Teeth", "Water")),
         Element = reorder(Element, ifelse(Dataset == "Teeth", CAP1, NA), 
                           na.rm = TRUE))

# plot----
pcap <-  ggplot(cap_compare_long, aes(x = CAP1, y = Element, fill = Dataset)) +
  geom_col(position = "dodge") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_fill_manual(values = c("Teeth" = "steelblue", "Water" = "orange3")) +
  labs(x = "CAP1 Score", y = NULL, fill = NULL) +
  theme(legend.position = c(0.2, 0.9),
        legend.direction = "horizontal")

ggsave(
  "partial_cap.png",
  plot = pcap,
  bg = "transparent",
  path = 'figures/grantrep/',
  scale = 1,
  width = 6,
  height = 5,
  dpi = 300)

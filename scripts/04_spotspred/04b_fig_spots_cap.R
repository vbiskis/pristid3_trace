# Details ----
#' 04b_spots_cap.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 7 Jul 2025
#' 
#' Content: cap visualisation
#' Update 24 Feb 2026 with helper scripts (microHelp.R)
#' -----------

source('helpers/microHelp.R')
set_theme(mytheme)
cap_data <- readxl::read_xlsx("data/primer/prim_results/CAP Summary.xlsx", sheet = "CAP1") 

microE <- read_xlsx('data/output/microE.xlsx')
microE <- microE %>% 
  mutate(across(where(is.character) & !all_of(c("ID", "File_ID")), as.factor))

# merge----
cap_data <- cap_data %>% 
  left_join(microE %>% 
              dplyr::select(ID, Region, Rsplit, RivGrp, Class, all_of(lnrats)), 
            by = c("Sample" = "ID")) %>% 
  mutate(Region = case_when(is.na(Region) ~ "NEQ", .default = Region),
         Rsplit = case_when(is.na(Rsplit) ~ "NEQ", .default = Rsplit),
         RivGrp = case_when(is.na(RivGrp) ~ "RG05", .default = RivGrp),
         Class = case_when(is.na(Class) ~ "YOY", .default = Class))

distAC <- vegdist(cap_data[, lnrats], method = "bray", 
                       na.rm = TRUE)

cap_resA <- capscale(distAC ~ Region, data = cap_data)
summary(cap_resA)

# Create the CAP plot
cyoy <- ggplot(cap_data, aes(x = CAP1, y = CAP2, 
                     color = factor(Region, 
                                    levels = c("GoC", "NEQ", "CEQ", "SEQ")), 
                     shape = factor(Region, 
                                    levels = c("GoC", "NEQ", "CEQ", "SEQ")))) +
  geom_point(size = 3, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    legend.position = c(0.1, 0.2),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  labs(x = "CAP1",  y = "CAP2") +
  guides(color = guide_legend(title = "Region"),
         shape = guide_legend(title = "Region")) +
  scale_shape_manual(
    values = c("CEQ" = 1, "NEQ" = 15, "GoC" = 17, "SEQ" = 16)
  ) +
  scale_color_manual(values = c("#D73027", "#FDAE61", "#66BD63", "#6495ED")) +
  coord_cartesian(ylim = c(-0.8, 0.8))

ggsave(
  "CAP-AC-YOY.png",
  plot = cyoy,
  bg = "transparent",
  device = NULL,
  path = 'figures/fig5/',
  scale = 1,
  width = 5,
  height = 4,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE
)

# numero 2----
cap_data2 <- readxl::read_xlsx("data/primer/prim_results/CAP Summary.xlsx", sheet = "CAP2") 

cap_data2 <- cap_data2 %>% 
  left_join(microE %>% 
              dplyr::select(ID, Region, Rsplit, RivGrp, Class, all_of(lnrats)), 
              by = c("Sample" = "ID")) %>% 
  mutate(Region = case_when(Region == "WCY" ~ "GoC", .default = Region))

cap_data2$Rsplit <- factor(cap_data2$Rsplit, 
                           levels = c("WCY", "NGoC", "SGoC", "FNEQ", "NEQ", "CEQ"))
cap_data2$Region <- factor(cap_data2$Region, 
                           levels = c("GoC", "FNEQ", "NEQ", "CEQ"))

corr_sq <- c(0.639, 0.4585, 0.0939) # from sheet

# Calculate % of fitted variation
total_fitted <- sum(corr_sq)
percent_fitted <- (corr_sq / total_fitted) * 100 # there's our % of fitted for cap 1, 2, 3

# If we want on plot labels? (unlikely)
cap1_percent <- round(percent_fitted[1], 1)
cap2_percent <- round(percent_fitted[2], 1)

# Create the CAP plot
cad <- ggplot(cap_data2, aes(x = CAP1, y = CAP2, color = Rsplit, shape = Region)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    legend.position = c(0.1, 0.2),
    legend.box = "horizontal", 
    plot.margin = margin(5, 5, 5, 5)
  ) +
  labs(
    # subtitle = "CAP1 correlation: 0.7994 | CAP2 correlation: 0.6771",
    x = "CAP1 (53.6% of fitted variation)", y = "CAP2 (38.5% of fitted variation)"
  ) +
  guides (color = "none", shape = guide_legend(title = "Region")) +
  scale_shape_manual(values = c("CEQ" = 1, "NEQ" = 15, "GoC" = 17)) +
  scale_color_manual(values = c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#A6D96A", "#66BD63")) +
  coord_cartesian(ylim = c(-0.8, 0.8))

ggsave(
  "CAP-AC-Adult.png",
  plot = cad,
  bg = "transparent",
  device = NULL,
  path = 'figures/fig5/',
  scale = 1,
  width = 5,
  height = 4,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE
)

# Details ----
#' 03b_spots_box.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 18 Aug 2024
#' 
#' Content: 
#' separation of signatures across key terms (but boxplot)
#' 
#' Update 24 Feb 2026:
#' pulled out the models section and moved to own script
#' + output cleaned up with theme in microHelp.R
#' ++ also specified which labels print to fit more graphs/page
#' -----------

source('helpers/microHelp.R')
set_theme(mytheme)
micro <- read_xlsx('data/output/icpms_clean.xlsx')
micro <- micro %>% 
  mutate(across(where(is.character) & !all_of(c("ID", "File_ID")), as.factor))

#Grouped Box----
melt <- micro %>%
  mutate(Region2 = case_when(Region == "WCY" ~ "GoC",
                             Rsplit == "FNEQ" ~ "FNQ",
                             .default = Region),
         Zone = factor(Zone, levels = c('EdgeE1', 'EdgeE2', 'Core', 'Tip'),
                       labels = c('E1', 'E2', 'C', 'T'))) %>% 
  dplyr::select(Species, Region2, Zone, 
                all_of(lnrats)) %>%
  gather(key = "Element", value = "Concentration", 
         -Species, -Region2, -Zone)

zone <- subset(melt, complete.cases(Zone))

zone_clean <- zone %>%
  group_by(Element, Species, 
           Region2, Zone) %>%
  mutate(
    Q1 = quantile(Concentration, 0.25, na.rm = TRUE),
    Q3 = quantile(Concentration, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower = Q1 - 3 * IQR,  # More conservative
    upper = Q3 + 3 * IQR
  ) %>%
  filter(Concentration >= lower & Concentration <= upper) %>%
  ungroup()

str_show <- c('lnBa')
xt_show <- c('lnLi')
lg_el <- c('lnMg', 'lnSr') # for axes

for (el in unique(zone_clean$Element)) {
  
  element_data <- subset(zone_clean, Element == el)
  clean_element <- gsub("^ln", "", el)
  element_data$Region2 <- factor(element_data$Region2, 
                                levels = c("AQ", "GoC", "FNQ", "NEQ", "CEQ", "SEQ"))
  element_data$Species <- factor(element_data$Species, 
                                 levels = c("Anoxypristis cuspidata", "Pristis pristis", "Pristis zijsron"),
                                 labels = c("A. cuspidata", "P. pristis", "P. zijsron"))
  
  p <- ggplot(element_data %>% 
                filter(Region2 != "SEQ"),
              aes(x = Zone, y = Concentration, fill = Region2)) +
    facet_wrap( ~ Species) +
    geom_boxplot(position = "dodge", na.rm = TRUE, outliers = FALSE) +
    scale_fill_manual(name = "Region", values = c("gray30", "#A50026", 
                                                  "#FEE08B", "#A6D96A",
                                                  "#006837", "#6495ED")) +
    labs(x = "Tooth Zone", 
         y = paste0(clean_element, ":Ca (", if(el %in% lg_el) "mmol" else "\u03bcmol", "/mol)")) +
    scale_y_continuous(labels = function(x) {
      exp_vals <- exp(x) - 1
      vals <- if(el %in% lg_el) exp_vals / 1000 else exp_vals
      signif(vals, digits = 2)
    }) +
    theme(
      strip.text = if(el %in% str_show) element_text(size = 12, face = "italic") else element_blank(),
      axis.text.x = if(el %in% xt_show) element_text(size = 11) else element_blank(),
      axis.title.x = if(el %in% xt_show) element_text(size = 12) else element_blank(),
      legend.position = if(el %in% xt_show) "bottom" else "none", 
      plot.margin = margin(5, 5, 5, 5)
    )
  
  ggsave(
    filename = paste0(clean_element, "_bp.png"),
    path = "figures/fig3/",
    plot = p,
    width = 8, 
    height = if(el %in% str_show) 2.25 else if(el %in% xt_show) 3.5 else 2,
    dpi = 300,
    bg = "transparent"
  )
  
  print(p)
}

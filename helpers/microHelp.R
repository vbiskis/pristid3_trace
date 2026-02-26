# Details ----
#' microHelp.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 20 Feb 2026
#' Content: General settings for plotting and stats
#' -----------

# libraries----
pacman::p_load('dplyr', 'tidyverse', 'readxl', 'plotly',
               'ggplot2', 'cowplot', 'reshape2',
               'patchwork', 'ggpubr', 
                # 'grid', 'gridExtra',
                # 'ggfortify', 'ggmosaic',
               'zoo', 'purrr', 'tsfeatures',
               'factoextra', 'lattice', 'cluster',
               'tm', 'psych', 'stringr', 'stats',
               'pilot', 'MetBrewer')

# repeated lists----

el <- c('Li_Ca', 'B_Ca', 'Na_Ca', 'Mg_Ca', 'Mn_Ca', 'Ni_Ca', 'Cu_Ca', 'Zn_Ca', 'Sr_Ca', 'Ba_Ca', 'Pb_Ca')
el_all <- c('In_Ca', el) # with In
els <- setdiff(el, "Na_Ca") # no Na

ln_rats <- paste0("ln", str_remove(el, "_Ca$"))
ln_cols <- c('lnIn', ln_rats) # with In
lnrats <- setdiff(ln_rats, "lnNa") # no Na

# colours----

spec_col <- c(
  "Anoxypristis cuspidata" = "blue4",
  #"Pristis clavata" = "gold2",
  "Pristis pristis" = "red3",
  "Pristis zijsron" = "darkgreen"
)

# labs----

spec_labs <- c(
  "Anoxypristis cuspidata" = expression(italic("A. cuspidata")), 
  #"Pristis clavata" = expression(italic("P. clavata")),
  "Pristis pristis" = expression(italic("P. pristis")),
  "Pristis zijsron" = expression(italic("P. zijsron"))
)

# theme----
suppressWarnings({mytheme <- pilot::theme_pilot(axis_title_size = 12, legend_title_size = 12,
                                                legend_title_family = "Optima", legend_text_family = "Optima",
                                                facet_title_family = "Optima", facet_title_size = 12,
                                                subtitle_family = "Optima", title_family = "Optima",
                                                axis_title_family = "Optima", axis_text_family = "Optima") 
}) # my god shes loud

tree <- function(elem, vars = c("Species", "Zone", "Class", "Region"), dat){
  f <- as.formula(paste(elem, "~", paste(vars, collapse = " + "))) # most common in there, but to allow editing
  t1 <- ctree(f, data = dat)
  plot(t1, gp = gpar(fontsize = 8), type = "simple", main = elem)
}

# Details ----
#' 05_profiles_aq.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 16 Jul 2024
#' Content: 
#' vert and teeth pics for aquarium samples
#' this one's manual, doesnt use the loop, can adjust
#' -----------

source('helpers/microHelp.R')
source('helpers/functions/fxn_mavs.R')
source('helpers/functions/fxn_profile.R')
profclean <- readRDS('data/output/prof_processed.RDS')
transsplit <- readRDS('data/output/transsplit.RDS')
names(transsplit)[54:60] 
# want 58-60 for SMA002
#58 tooth bleach
#59 tooth water
#60 vert bleach 

#issues with vert water (in general not as high performing)

#SMA002----
SMA002TB <- transsplit[[58]] 

# bleach ----
# lets remove In spikes
b2 <- plot_holes(SMA002TB, trigger =  3, release = 1.5,
                 elements = c("Sr_Ca", "Ba_Ca", "Li_Ca", "Mg_Ca")) # lol welp

SMA002TB <- SMA002TB %>% 
  filter(# true_dist >= min(true_dist[which(In_Ca < 5)]), # this tooth is rough
         true_dist <= 49907)

SMA002TB_out <- profs(SMA002TB, el_all, "tooth", span = 0.01, window_size = 11, # stats
                      show_In = TRUE, In_style = "points", 
                      rect_threshold = 50) 
SMA002TB_out$plots$Ca

sma02a2 <- SMA002TB_out$plots$Sr + theme(axis.title.x = element_blank())
sma02b2 <- SMA002TB_out$plots$Ba + theme(axis.title.x = element_blank())
sma02c2 <- SMA002TB_out$plots$Li + theme(axis.title.x = element_blank())
sma02d2 <- SMA002TB_out$plots$Mg + theme(axis.title.x = element_blank())
sma02e2 <- SMA002TB_out$plots$Mn 

# water----
SMA002TH <- transsplit[[59]]

h2 <- plot_holes(SMA002TH, trigger =  3, release = 1.5,
                 elements = c("Sr_Ca", "Ba_Ca", "Li_Ca", "Mg_Ca"))

SMA002TH <- h2$data[!h2$data$.hole, ] # lol! k

## holes and smoothing----
# red = one loess across everything (how it was); blue = each piece on its own (now)
find_gaps(SMA002TH$true_dist) # holes > 50 um (normal spacing ~8.6 um)
plot_smoothing(SMA002TH, c("Sr_Ca", "Ba_Ca", "Li_Ca", "Mg_Ca"), span = 0.05)
plot_smoothing(SMA002TH, c("Sr_Ca", "Ba_Ca", "Li_Ca", "Mg_Ca"), span = 0.05, xlim = c(13000, 25000)) # the big hole
# data right at the hole edges still spike: trim_gaps() drops points within buffer um of a hole
plot_smoothing(trim_gaps(SMA002TH, buffer = 400), c("Sr_Ca", "Ba_Ca", "Li_Ca", "Mg_Ca"), span = 0.05, xlim = c(13000, 25000))

SMA002TH_out <- profs(SMA002TH, el_all, "tooth", span = 0.01, window_size = 11,
                      # img_path = "data/transects/imgs/SMA002/SMA002-Trsct.jpg", # still big
                      show_In = TRUE, In_style = "points", ) 
SMA002TH_out$plots$Ca

sma02a <- SMA002TH_out$plots$Sr + theme(axis.title.x = element_blank())
sma02b <- SMA002TH_out$plots$Ba + theme(axis.title.x = element_blank())
sma02c <- SMA002TH_out$plots$Li + theme(axis.title.x = element_blank())
sma02d <- SMA002TH_out$plots$Mg + theme(axis.title.x = element_blank())
sma02e <- SMA002TH_out$plots$Mn 

# vert----
SMA002V <- transsplit[[60]]
SMA002V <- SMA002V %>% 
  filter(true_dist >= true_dist[which(In_Ca < 10)][1], # adjust to see window move (2, 5, 10)
         true_dist <= 7950) # bc it comes back down, manual here

SMA002V_out <- profs(SMA002V, el_all, "vert", span = 0.025, show_In = TRUE) # show Patrick!
SMA002V_out$plots$Ca

sma02va <- SMA002V_out$plots$Sr + theme(axis.title.x = element_blank())
sma02vb <- SMA002V_out$plots$Ba + theme(axis.title.x = element_blank())
sma02vc <- SMA002V_out$plots$Li + theme(axis.title.x = element_blank())
sma02vd <- SMA002V_out$plots$Mg + theme(axis.title.x = element_blank())
sma02ve <- SMA002V_out$plots$Mn 

sma002 <- (sma02a + sma02va) / (sma02b + sma02vb) / (sma02c + sma02vc) / 
  (sma02d + sma02vd) & theme(plot.margin = margin(5, 5, 5, 5))

ggsave(
  "SMA002VT.png",
  plot = sma002,
  device = NULL,
  path = 'figures/supp/s4',
  scale = 1,
  width = 8,
  height = 8,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
)

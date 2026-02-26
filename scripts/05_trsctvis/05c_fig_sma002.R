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
sapply(transsplit, function(df) unique(df$ID))

# want 58-60 for SMA002
#58 tooth bleach
#59 tooth water
#60 vert bleach 

#issues with vert water (in general not as high performing)

#SMA002----
SMA002TB <- transsplit[[58]] 
SMA002TB <- SMA002TB %>% 
  filter(true_dist >= min(true_dist[which(In_Ca < 5)]), # this tooth is rough
         true_dist <= max(true_dist[which(In_Ca < 5)])) # so bad we get very little of it, no matter the val

SMA002TB_filt <- profs(SMA002TB, el_all, "tooth",
                      show_In = TRUE) 

SMA002TH <- transsplit[[59]] 
SMA002TH <- SMA002TH %>% 
  filter(true_dist >= 130, #looks good! just to take the first 5 sec off
         true_dist <= max(true_dist[which(In_Ca < 5)])) 

SMA002TH_filt <- profs(SMA002TH, el_all, "tooth", span = 0.2, thin = 5,
                      # img_path = "data/transects/imgs/SMA002/SMA002-Trsct.jpg", # still big
                      show_In = TRUE) 

sma02a <- a + theme(axis.title.x = element_blank())
sma02b <- b + theme(axis.title.x = element_blank())
sma02c <- cc + theme(axis.title.x = element_blank())
sma02d <- d + theme(axis.title.x = element_blank())
sma02e <- e
sma02f <- f 

#vert specs
SMA002V <- transsplit[[60]]
SMA002V <- SMA002V %>% 
  filter(true_dist >= true_dist[which(In_Ca < 10)][1], # adjust to see window move (2, 5, 10)
         true_dist <= 7900) # bc it comes back down, manual here

SMA002V_filt <- profs(SMA002V, el_all, "vert", span = 0.1, show_In = TRUE) # show Patrick!

sma02va <- a + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
sma02vb <- b + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
sma02vc <- cc + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
sma02vd <- d + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
sma02ve <- e 
sma02vf <- f

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

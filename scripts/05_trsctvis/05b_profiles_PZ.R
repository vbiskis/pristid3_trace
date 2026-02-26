# Details ----
#' 05_profiles_PZ.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 18 Aug 2024
#' Content: 
#' comparison of section planes for P. zijsron
#' -----------
#' 
source('helpers/microHelp.R')
source('helpers/functions/fxn_mavs.R')
source('helpers/functions/fxn_profile.R')
profclean <- readRDS('data/output/prof_processed.RDS')

#PZ----
##nb028----
PZPro <- profclean %>% 
  filter(Species == "zijsron")
unique(PZPro$SpecID)

NB028w <- PZPro %>% 
  filter(SpecID == "NB028", Bath == "water", Loc == "full",
         true_dist >= ca_st, true_dist <= ca_end)

#reverse these before upload
nb028_w <- profs(NB028w, el_all, "tooth", 
                        img_path = 'data/transects/imgs/NB028/NB028H-Full-flip.jpg')

nb028_wa <- a + theme(axis.title.x = element_blank())
nb028_wb <- b + theme(axis.title.x = element_blank())
nb028_wc <- c + theme(axis.title.x = element_blank())
nb028_wd <- d 
nb028_we <- e 
nb028_wi <- f

NB028b <- PZPro %>% 
  filter(SpecID == "NB028", Bath == "bleach", Loc == "full",
         true_dist >= ca_st, true_dist <= ca_end,
         Li_Ca <= 150) # big spike ?

nb028_b <- profs(NB028b, el_all, "tooth") # well its not as bad as I thought?

##437----

pz437 <- PZPro %>% 
  filter(SpecID == "437", Loc == "mid",
         true_dist >= ca_st, true_dist <= ca_end) 

pz437_w <- profs(pz437, el_all, "tooth", span = 0.1) # is that even a mirror?

##3150----

pz3150 <- PZPro %>% 
  filter(SpecID == "3150", Loc == "mid",
         true_dist >= ca_st, true_dist <= ca_end)

pz3150_w <- profs(pz3150, el_all, "tooth", span = 0.1) 

##3149----

pz3149 <- PZPro %>% 
  filter(SpecID == "3149",
         true_dist >= ca_st, true_dist <= ca_end)

pz3149_w <- profs(pz3149, el_all, "tooth", span = 0.1)

##3151----

pz3151 <- PZPro %>% 
  filter(SpecID == "3151",
         true_dist >= ca_st, true_dist <= ca_end)

pz3151_w <- profs(pz3151, el_all, "tooth", span = 0.1)

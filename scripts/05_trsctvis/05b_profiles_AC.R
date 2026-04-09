# Details ----
#' 05_profiles_AC.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 18 Aug 2024
#' Content: 
#' comparison of section planes for A. cuspidata
#' + verts/teeth 
#' -----------

source('helpers/microHelp.R')
source('helpers/functions/fxn_mavs.R')
source('helpers/functions/fxn_profile.R')
profclean <- readRDS('data/output/prof_processed.RDS')

# AC----
ACPro <- profclean %>% 
  filter(Species == "anoxy")
unique(ACPro$ID) #all profiles, including horizontal

# run profiles----
## a4224----
a4224 <- ACPro %>% 
  filter(SpecID == "4224")

a4224 <- a4224 %>% 
  filter(true_dist >= ca_st, true_dist <= ca_end) # it works! going off values of Ca better than In

a4224t <- profs(a4224, el_all, "tooth", span = 0.1
                #img_path = 'data/transects/imgs/QBFP004224/004224-Full-fwd.jpg'
                ) 
#can leave out if not available

## vert comparison ----
# def going for a loop (there's 6)

an_ids <- c("480", "486", "482", "483", "484", "485")
In_thresh <- 20

bbanoxy <- list()

for(id in an_ids){
  for(bath in c("water", "bleach")){
    
    dat_v <- ACPro %>%
      filter(SpecID == id, Bath == bath, Type == "vert") %>%
      filter(true_dist <= true_dist[tail(which(In_Ca < In_thresh), 1)])
    if(nrow(dat_v) > 0){
      profs(dat_v, el_all, "vert")
      bbanoxy[[paste0(id, "v", bath)]] <- list(a=a, b=b, cc=cc, d=d, e=e, f=f)
    }
    
    dat_t <- ACPro %>%
      filter(SpecID == id, Bath == bath, Type == "tooth",
             true_dist >= ca_st, true_dist <= ca_end)
    if(nrow(dat_t) > 0){
      profs(dat_t, el_all, "tooth")
      bbanoxy[[paste0(id, "t", bath)]] <- list(a=a, b=b, cc=cc, d=d, e=e, f=f)
    }
  }
}

### compare! ----
(bbanoxy[["480twater"]]$a + bbanoxy[["480vwater"]]$a) /
  (bbanoxy[["480tbleach"]]$a + bbanoxy[["480vbleach"]]$a) &
  theme(axis.title.x = element_blank())

(bbanoxy[["482twater"]]$a + bbanoxy[["482vwater"]]$a) /
  (bbanoxy[["482tbleach"]]$a + bbanoxy[["482vbleach"]]$a) &
  theme(axis.title.x = element_blank())

(bbanoxy[["483twater"]]$a + bbanoxy[["483vwater"]]$a) /
  (bbanoxy[["483tbleach"]]$a + bbanoxy[["483vbleach"]]$a) &
  theme(axis.title.x = element_blank())

(bbanoxy[["484twater"]]$a + bbanoxy[["484vwater"]]$a) /
  (bbanoxy[["484tbleach"]]$a + bbanoxy[["484vbleach"]]$a) &
  theme(axis.title.x = element_blank())

(bbanoxy[["485twater"]]$a + bbanoxy[["485vwater"]]$a) /
  (bbanoxy[["485tbleach"]]$a + bbanoxy[["485vbleach"]]$a) &
  theme(axis.title.x = element_blank())

(bbanoxy[["486twater"]]$a + bbanoxy[["486vwater"]]$a) /
  (bbanoxy[["486tbleach"]]$a + bbanoxy[["486vbleach"]]$a) &
  theme(axis.title.x = element_blank())

### 486----
#bc he has pics

a486 <- ACPro %>% 
  filter(SpecID == "486", Type == "tooth", Bath == "water", 
         true_dist >= ca_st, true_dist <= ca_end)

a486t <- profs(a486, el_all, "tooth", 
                img_path = 'data/transects/imgs/QBFP000486/486HTMTrsct.jpg')

a486 <- ACPro %>% # ya shows a mirror with In
  filter(SpecID == "486", Type == "vert", Bath == "water") %>% # note - need bleach ones
  filter(true_dist <= true_dist[tail(which(In_Ca < 10), 1)])

a486v <- profs(a486, el_all, "vert", 
               img_path = 'data/transects/imgs/QBFP000486/486HV-6.3M.jpg')

## just teeth----

an_ids2 <- c("3158", "055", "448", "4428") # all the partials (mid etc.)
anoxyteeth <- list()

for(id in an_ids2){
    
    dat_m <- ACPro %>%
      filter(SpecID == id, Loc == "mid",
             distance >= ca_st, distance <= ca_end)
    if(nrow(dat_m) > 0){
      profs(dat_m, el_all, "tooth")
      anoxyteeth[[paste0(id, "tmid")]] <- list(a=a, b=b, cc=cc, d=d, e=e, f=f)
    }
    
    dat_b <- ACPro %>%
      filter(SpecID == id, Loc == "base",
             distance >= ca_st, distance <= ca_end)
    if(nrow(dat_b) > 0){
      profs(dat_b, el_all, "tooth")
      anoxyteeth[[paste0(id, "tbase")]] <- list(a=a, b=b, cc=cc, d=d, e=e, f=f)
    }
}

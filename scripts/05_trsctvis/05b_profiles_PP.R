# Details ----
#' 05_profiles_PP.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 18 Aug 2024
#' Content: 
#' comparison of teeth, P. pristis
#' ( excl. aquaria)
#' -----------

source('helpers/microHelp.R')
source('helpers/functions/fxn_mavs.R')
source('helpers/functions/fxn_profile.R')
profclean <- readRDS('data/output/prof_processed.RDS')

PPPro <- profclean %>% 
  filter(Species == "pristis")

unique(PPPro$ID) # a cool comparison oppthere!

# just mitchell teeth----

##4462----
p4462 <- PPPro %>% 
  filter(SpecID == "4462", 
         true_dist >= ca_st, true_dist <= ca_end)

px4462 <- profs(p4462, el_all, "tooth", 
                # img_path = 'data/transects/imgs/QBFP004462/004462-crop.jpg',
                span = 0.1)

##4468----
p4468 <- PPPro %>% 
  filter(SpecID == "4468", 
         true_dist >= ca_st, true_dist <= 41500)

px4468<- profs(p4468, el_all, "tooth", span = 0.05)

##4469----
prist_id <- c("10_4469FullT1 - 01.csv", 
              "30_4469FullT2 - 01.csv", 
              "30_4469HTip - 01.csv", 
              "30_4469HMid - 01.csv",
              "30_4469HBase - 01.csv", # lol toss
              "10_4469CrossT1 - 01.csv",
              "30_4469CrossT2 - 01.csv") 

pristeeth <- list()

for(id in prist_id){
  
  dat_p <- PPPro %>%
    filter(ID == id, 
           true_dist >= min(true_dist[which(In_Ca < 5)]), 
           true_dist <= max(true_dist[which(In_Ca < 5)]))
  
  if(nrow(dat_p) > 0){
    profs(dat_p, el_all, "tooth", span = 0.05)
    pristeeth[[gsub("^[0-9]+_| - 01\\.csv$", "", id)]] <- list(a=a, b=b, cc=cc, d=d, e=e, f=f)
  }
}

# slow v fast
# slower / faster
# slower / faster [cross]

(pristeeth[["4469FullT1"]]$a + pristeeth[["4469FullT2"]]$a) /
  (pristeeth[["4469CrossT1"]]$a + pristeeth[["4469CrossT2"]]$a) &
  theme(axis.title.x = element_blank())

(pristeeth[["4469FullT1"]]$b + pristeeth[["4469FullT2"]]$b) /
  (pristeeth[["4469CrossT1"]]$b + pristeeth[["4469CrossT2"]]$b) &
  theme(axis.title.x = element_blank())

(pristeeth[["4469FullT1"]]$cc + pristeeth[["4469FullT2"]]$cc) /
  (pristeeth[["4469CrossT1"]]$cc + pristeeth[["4469CrossT2"]]$cc) &
  theme(axis.title.x = element_blank())

(pristeeth[["4469FullT1"]]$d + pristeeth[["4469FullT2"]]$d) /
  (pristeeth[["4469CrossT1"]]$d + pristeeth[["4469CrossT2"]]$d) &
  theme(axis.title.x = element_blank())

(pristeeth[["4469FullT1"]]$e + pristeeth[["4469FullT2"]]$e) /
  (pristeeth[["4469CrossT1"]]$e + pristeeth[["4469CrossT2"]]$e) &
  theme(axis.title.x = element_blank())

(pristeeth[["4469FullT1"]]$f + pristeeth[["4469FullT2"]]$f) /
  (pristeeth[["4469CrossT1"]]$f + pristeeth[["4469CrossT2"]]$f) &
  theme(axis.title.x = element_blank())

#and: all fast
# full / base
# tip / mid

(pristeeth[["4469FullT2"]]$a + pristeeth[["4469CrossT2"]]$a) /
  (pristeeth[["4469HTip"]]$a + pristeeth[["4469HMid"]]$a) &
  theme(axis.title.x = element_blank())

(pristeeth[["4469FullT2"]]$b + pristeeth[["4469CrossT2"]]$b) /
  (pristeeth[["4469HTip"]]$b + pristeeth[["4469HMid"]]$b)&
  theme(axis.title.x = element_blank())

(pristeeth[["4469FullT2"]]$cc + pristeeth[["4469CrossT2"]]$cc) /
  (pristeeth[["4469HTip"]]$cc + pristeeth[["4469HMid"]]$cc)&
  theme(axis.title.x = element_blank())

(pristeeth[["4469FullT2"]]$d + pristeeth[["4469CrossT2"]]$d) /
  (pristeeth[["4469HTip"]]$d + pristeeth[["4469HMid"]]$d)&
  theme(axis.title.x = element_blank())

(pristeeth[["4469FullT2"]]$e + pristeeth[["4469CrossT2"]]$e) /
  (pristeeth[["4469HTip"]]$e + pristeeth[["4469HMid"]]$e)&
  theme(axis.title.x = element_blank())

(pristeeth[["4469FullT2"]]$f + pristeeth[["4469CrossT2"]]$f) /
  (pristeeth[["4469HTip"]]$f + pristeeth[["4469HMid"]]$f)&
  theme(axis.title.x = element_blank())

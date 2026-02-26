# Details ----
#' 05_get_profiles.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 6 Sep 2024
#' Content: 
#' Extract individual profiles from prepped df
#' -----------

source('helpers/microHelp.R')
profiles <- read_xlsx('data/output/profiles.xlsx')
library(writexl)

#fix df----
transdat <- profiles %>% 
  mutate(SpecID = case_when(str_detect(ID, "_[0-9]") ~ str_extract(ID, "(?<=_)[0-9]+"),
                           str_detect(ID, "SMA") ~ str_extract(ID, "SMA[0-9]+"),
                           .default = "NB028"),
         Type = case_when(str_detect(ID, "V") ~ "vert",
                                      .default = "tooth"),
         Loc = case_when(str_detect(ID, "Base") ~ "base",
                         str_detect(ID, "Mid") ~ "mid",
                         str_detect(ID, "Tip") ~ "tip",
                         str_detect(ID, "Cross") ~ "cross",
                         .default = "full"),
         Bath = case_when(ID == "30_NB028Transect - 01.csv" | str_detect(ID, "[0-9]H") ~ "water",
                          str_detect(ID, "[0-9]B") ~ "bleach",
                          .default = "bleach"),
         Species = case_when(SpecID %in% c("SMA001", "SMA002", "4462", "4468", "4469") ~ "pristis",
                             SpecID %in% c("437", "NB028", "3149", "3150", "3151") ~ "zijsron",
                             .default = "anoxy"),
         Speed = gsub("_.*", "", ID)) %>% 
  dplyr::select(SpecID, Species, Type, Loc, Bath, everything())

transdat %>%
  group_by(SpecID) %>%
  reframe(Spec = unique(Species),
          Speeds = unique(Speed),
          Types = unique(Type),
          Locs = paste(unique(Loc), collapse = ", "),
          Baths = unique(Bath)) %>%
  print(n = Inf) %>% 
  ungroup() #so 43 + 2xspeed 4469 + 10 from bb anoxy + 4 from SMA = 59

unique(transdat$ID) #62 total - OKAY we match

transum <- transdat %>% 
  group_by(ID) %>% 
  summarise(SpecID = unique(SpecID),
            Spec = unique(Species),
            Type = unique(Type),
            Loc = unique(Loc),
            Bath = unique(Bath)) %>% 
  dplyr::select(ID, SpecID, Spec, Type, Loc, Bath)

# earlier teeth were run backwards! lets fix:
revids <- unique(transum %>% filter(Loc == "full", Type == "tooth", Bath == "bleach", Spec != "anoxy"))$ID
#but this misses two
revids_all <- c(revids, "30_NB028Transect - 01.csv", "30_4224Transect - 01.csv") # older green + anoxy

transdat <- transdat %>%
  group_by(ID) %>%
  mutate(true_dist = case_when(ID %in% revids_all ~ max(distance) - distance,
                               .default = distance)) %>%
  ungroup()

#Using inflection In_Ca----
transdatall <- transdat %>% 
  group_by(ID) %>%
  mutate(start = ifelse(any(In_Ca < 10), min(true_dist[which(In_Ca < 10)]), NA),
         end = ifelse(any(In_Ca < 10), max(true_dist[which(In_Ca < 10)]), NA),
         ca_st = ifelse(any(Ca43_CPS > 3300000), min(true_dist[which(Ca43_CPS > 3300000)]), NA),
         ca_end = ifelse(any(Ca43_CPS > 3300000), max(true_dist[which(Ca43_CPS > 3300000)]), NA)) %>% 
  ungroup() %>%
  dplyr::select(ID, SpecID, Species, Type, Loc, Bath, start, ca_st, true_dist, everything())

transsplit <- transdat %>% 
  group_split(SpecID, Type, Loc, Bath) 

saveRDS(object = transdatall, file = "data/output/prof_processed.RDS") # all one, cleaned
saveRDS(object = transsplit, file = "data/output/transsplit.RDS") #individual dfs
write_xlsx(transum, "data/output/transum.xlsx") #summary of all profs


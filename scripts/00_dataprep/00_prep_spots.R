# Details ----
#' 00_prep_spots.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 
#' - Origin: 18 Aug 2024
#' - Updated: 20 Feb 2026
#' 
#' Content: 
#' + combine icpms data with saw information dfs
#' + filtering for spot locations to be used in analysis
#' ++ diagnostics come later
#' -----------

source('helpers/microHelp.R')
library(writexl)

#Load Data (all teeth with B reading oob have been removed)
icpms <- readxl::read_xlsx("data/AllMicro.xlsx", sheet = "Corrected") 
sawinf <- readxl::read_xlsx("data/Saw_Summ.xlsx", sheet = "Summary")

micro <- sawinf %>%
  dplyr::select("ID", "Species", "LAT", "LONG", "r", "Year", "Tip_Year", 
                "Region", "Rsplit", "RivGrp", "EW", "NS", 
                "TL_MP", "Age_SP", "Class") %>%
  dplyr::right_join(icpms, by = "ID")

write_xlsx(micro, path = "data/output/micro.xlsx" ) # keep full version for later ref

# lets make some new sheets for visibility and automate the spot selection

# rows by sawID, teeth only first, bleached condition
# indium vals for each spot type

spots_find <- icpms %>%
  filter(Type == 'Tooth',
         Cut == 'Cor') %>% # all same cut
  group_by(ID, Zone) %>%
  dplyr::summarise(File_min = File_ID[which.min(In_Ca)],
            min_In = min(In_Ca, na.rm = TRUE),
            File_six = File_ID[Size=='60'][which.min(In_Ca)],
            In_six = In_Ca[match(File_six, File_ID)]) %>%
  ungroup() %>%
  mutate(File_fin = case_when(is.na(In_six) ~ File_min, # no 60s available
                       In_six < 0.2 | In_six < 1.5 * min_In ~ File_six, # diff is not sig.
                       TRUE ~ File_min)) # if In is high in 60, we go with lowest.

spots_test <- spots_find %>%
  left_join(icpms %>%
              dplyr::select(File_ID, In_Ca), 
            by = c("File_fin" = "File_ID")) %>%
  rename(In_fin = In_Ca)

find_extra <- icpms %>%
  group_by(ID) %>%
  dplyr::summarise(
    BW = any(BW == "B") & any(BW == "W"),
    Vert = any(Type == "Vert")
  )

spots_look <- spots_test %>%
  dplyr::select(ID, Zone, File_fin, In_fin) %>%
  pivot_wider(
    names_from = Zone,
    values_from = c(In_fin, File_fin))%>%
  right_join(find_extra, by = "ID") %>%
  mutate(
    Dat = case_when(
      !is.na(File_fin_Tip) & !is.na(File_fin_Core) & 
        (!is.na(File_fin_EdgeE1) & In_fin_EdgeE1 < 1 |
           !is.na(File_fin_EdgeE2) & In_fin_EdgeE2 < 1) ~ 'All',
      (is.na(File_fin_Tip) | In_fin_Tip > 1) & 
        ((!is.na(File_fin_EdgeE1) & In_fin_EdgeE1 < 1) |
           (!is.na(File_fin_EdgeE2) & In_fin_EdgeE2 < 1))  ~ 'Capture Only',
      (is.na(File_fin_EdgeE1) | In_fin_EdgeE1 > 1) &
        (is.na(File_fin_EdgeE2) | In_fin_EdgeE2 > 1) &
        !is.na(File_fin_Tip) & In_fin_Tip < 1 ~ 'Natal Only',
      ((!is.na(File_fin_EdgeE1) & In_fin_EdgeE1 < 1) |
         (!is.na(File_fin_EdgeE2) & In_fin_EdgeE2 < 1)) &
        (!is.na(File_fin_Tip) & In_fin_Tip < 1) ~ 'Limited',
      TRUE ~ 'Rerun'
    )
  ) # for viewing success of runs

saw_sum <- sawinf %>%
  dplyr::select("ID", "Species", "Profile", "LAT", "LONG", "Year", "Region", "RivGrp", "TL_MP", "Age_SP", "Class") %>%
  dplyr::right_join(spots_look%>%
                      dplyr::select(ID,
                             File_fin_EdgeE1,
                             File_fin_EdgeE2,
                             File_fin_Core,
                             File_fin_Tip,
                             BW,
                             Vert,
                             Dat),
                    by = "ID") #the determine if missing samples are unique enough that we NEED them

spots_sum <- micro %>%
  filter(File_ID %in% spots_test$File_fin)%>%
  bind_rows(
    micro %>%
      filter(Rsplit == 'AQ', Zone == 'Tip') %>%
      mutate( #adding extra rows for our aquarium specimens to be included in regional study
        Class = "IMMATURE",
        Zone = 'EdgeE2',
        Zone_2 = 'Edge',
        Region = 'GoC',
        Rsplit = 'SGoC',
        RivGrp = 'RG09'
      )
    )

MicroR <- subset(spots_sum, Year >= 1988)
MicroFil <- subset(MicroR, In_Ca < 1)
MicroAdj <- MicroFil %>% 
  mutate(across(all_of(el_all), ~log(.x + 1), .names = "ln{str_remove(.col, '_Ca$')}"),
         Maturity = if_else(Class != "ADULT", "Immature", "Mature"))

#now have a filtered df of spots chosen based on lowest Indium
write_xlsx(MicroAdj, path = "data/output/icpms_clean.xlsx" ) # export clean version

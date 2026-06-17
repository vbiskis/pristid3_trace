# Details ----
#' 07_env_cleandf.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 17 Jun 2026
#' 
#' Content: 
#' + mapping study areas from DAF and DETSI
#' + prelim identification of trends
#' -----------

# source('helpers/mapHelp.R')
source('helpers/microHelp.R')
library(writexl)
set_theme(mytheme)

# Daf first----
daf <- read_xlsx('data/enviro/daf/forUSC_all wild water samples.xlsx')

## Clean it up----
# extract the numeric value after <
daf_lod <- daf %>% 
  mutate(across(12:ncol(.), ~ if_else(
    str_starts(., "<"), 
    as.double(str_remove(., "<")) * 0.5,  # half-LOD substitution
    as.double(.)
  )))

glimpse(daf_lod) # eew

daf_cl <- daf_lod %>% select(where(~!all(is.na(.)))) # k lost 24
# and now ratios 
colnames(daf_cl) # of course they're diff units

ug_map <- c(
  "Lithium" = "Li", "Boron" = "B", "Manganese" = "Mn",
  "Nickel" = "Ni", "Copper" = "Cu", "Zinc" = "Zn",
  "Strontium" = "Sr", "Barium" = "Ba", "Lead" = "Pb"
)

mg_map <- c(
  "Sodium" = "Na", "Magnesium" = "Mg"
)

daf_cl <- daf_cl %>%
  mutate(
    # µg/L ratios
    across(
      all_of(paste0(names(ug_map), "_µg/L")),
      ~ . / (`Calcium_mg/L` * 1000),
      .names = "{ug_map[str_remove(.col, '_µg/L')]}_Ca"
    ),
    # mg/L ratios
    across(
      all_of(paste0(names(mg_map), "_mg/L")),
      ~ . / `Calcium_mg/L`,
      .names = "{mg_map[str_remove(.col, '_mg/L')]}_Ca"
    )
  )

daf_good <- daf_cl %>% 
  select(2:3, 5:8 , 10:11, all_of(el))

daf_good <- daf_good %>%
  mutate(across(all_of(el), ~ if_else(is.na(.), min(., na.rm = TRUE), .)))

glimpse(daf_good)

write_xlsx(daf_good, path = "data/output/daf_clean.xlsx")

# and detsi----

des <- read_xlsx("data/enviro/detsi/NESP WQ - Complete dataset.xlsx", sheet = 2, col_names = FALSE)
# save headers
method_row <- as.character(des[1,])
element_row <- as.character(des[2,])
unit_row <- as.character(des[3,])

# fix col_names
names(des) <- paste0(element_row, "_", unit_row)
des2 <- des[5:nrow(des),] # drop header rows, keep data
utd_cols <- which(method_row == "W_UTD_IMS") # keep site cols + W_UTD_IMS cols
md_icp_cols <- which(method_row == "W_MD_ICP" & element_row %in% c("Calcium", "Magnesium", "Sodium"))
td_icp_cols <- which(method_row == "W_TD_ICP" & element_row == "Strontium")
des2 <- des[5:nrow(des), c(1, 2, 3, md_icp_cols, td_icp_cols, utd_cols)]

names(des2)[1:3] <- c("Site_Number", "Site_Name", "Date")
des2$Date <- as.Date(as.numeric(des2$Date), origin = "1899-12-30")

des2 <- des2 %>% 
  mutate(across(4:ncol(.), ~ if_else(
    str_starts(., "<"), 
    as.double(str_remove(., "<")) * 0.5,  # half-LOD substitution
    as.double(.)
  )))

# and fix sr... 
des2 <- des2 %>% 
  mutate(`Strontium_µg/L` = if_else(is.na(`Strontium_µg/L`), `Strontium_mg/L`*1000, `Strontium_µg/L`))

des2 <- des2 %>% filter(!is.na(`Calcium_mg/L`))

# omg really why no long lat you butt

des_sites <- read_xlsx("data/enviro/detsi/Strontium water samples_NESP Northern hub_project 1.3.1.xlsx")
des_sites <- des_sites %>% 
  select(`QG Site number`, `Collection date`, Lat, Long) %>% 
  mutate(Date = as.Date(as.numeric(`Collection date`), origin = "1899-12-30"),
         Site_Number = `QG Site number`)

des3 <- des2 %>% 
  left_join(des_sites %>% select('Site_Number', 'Date', 'Lat', 'Long'))

# one guy left out from PR goldfields, mine site (-16.092179, 144.321952)
des3 <- des3 %>% 
  mutate(Lat = if_else(is.na(Lat), -16.092179, Lat), Long = if_else(is.na(Long), 144.321952, Long))

# phew!
# k lets get those ratios...

ug_map <- c("Lithium" = "Li", "Boron" = "B", "Manganese" = "Mn", "Nickel" = "Ni", "Copper" = "Cu", "Zinc" = "Zn", "Strontium" = "Sr", "Barium" = "Ba", "Lead" = "Pb")

mg_map <- c("Sodium" = "Na", "Magnesium" = "Mg")

des_cl <- des3 %>%
  mutate(
    # µg/L ratios
    across(
      all_of(paste0(names(ug_map), "_µg/L")),
      ~ . / (`Calcium_mg/L` * 1000),
      .names = "{ug_map[str_remove(.col, '_µg/L')]}_Ca"
    ),
    # mg/L ratios
    across(
      all_of(paste0(names(mg_map), "_mg/L")),
      ~ . / `Calcium_mg/L`,
      .names = "{mg_map[str_remove(.col, '_mg/L')]}_Ca"
    )
  )

des_good <- des_cl %>% 
  select(1:3, all_of(el), Lat, Long)

des_good <- des_good %>%
  mutate(across(all_of(el), ~ if_else(is.na(.), min(., na.rm = TRUE), .)))

#gah finally!!!
write_xlsx(des_good, path = "data/output/des_clean.xlsx")


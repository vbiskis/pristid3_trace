# Details ----
#' 00_prep_trsct.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Dr Patrick Reis-Santos
#' modified by Nikki Biskis
#' Date: 18 Aug 2024
#' 
#' Content: 
#' data wrangling
#' ++ profile df
#' -----------

source('helpers/microHelp.R')
library(writexl)
#transects are in folders based on their speed in microns/s

# we have already processed many:
csv_files <- list.files(path = 'data/transects/processed', 
                        pattern = "\\.csv$", full.names = TRUE)

# the following functions have been adjusted to skip files that are already processed
# Append speed to file name----
file_dir <- "data/transects/profile_dat"
proc_dir <- "data/transects/processed"
folders <- list.dirs(path = file_dir, full.names = TRUE, recursive = FALSE)

# Convert LA-ICPMS files (.d) to .csv----
d_to_csv <- function(folder_path) {
  folder_name <- basename(folder_path)
  files <- list.files(folder_path, full.names = TRUE)
  files <- files[grepl("\\.d$", files)]  # only .d files
  
  processed <- gsub(paste0("^", folder_name, "_"), "", # the prefix we're appending
                    tools::file_path_sans_ext(basename(csv_files))) # what I named processed files above
  files <- files[!tools::file_path_sans_ext(basename(files)) %in% processed]
  
  if (length(files) == 0) {
    cat("No new files to convert in", folder_name, "\n")
    return(invisible(NULL)) # clever
  }
  
  for (file in files) {
    data <- read.delim(file)
    file_name <- paste0(folder_name, "_", tools::file_path_sans_ext(basename(file))) # add the number
    csv_file <- file.path(proc_dir, paste0(file_name, ".csv")) # make a csv
    write.csv(data, file = csv_file, row.names = FALSE) # write it
    cat("Converted", file, "to", csv_file, "\n")
  }
}

subfolders <- c("5", "10", "30")
for (subfolder in subfolders) {
  folder_path <- file.path(file_dir, subfolder)
  d_to_csv(folder_path)
}

# Importing transect data----
# Start here if files renamed already!

data1 <- data.frame()

# Import list of processed CSV files, do again even if ran before
csv_files <- list.files(path = 'data/transects/processed', 
                        pattern = "\\.csv$", full.names = TRUE)

for (file in csv_files) {
  sheet <- read.csv(file, header = TRUE)
  file_name <- basename(file)
  prefix <- strsplit(file_name, "_")[[1]][1]
  multiplier <- switch(prefix,
                       "5" = 5,
                       "10" = 10,
                       "30" = 30,
                       NA) # Default to NA if no match
  if (is.na(multiplier)) {
    warning("Unrecognized prefix in file: ", file_name)
    next
  }
  sheet <- sheet %>%
   mutate(
      ID = file_name,
      ElapsedTime_s = as.numeric(ElapsedTime_s),
      speed = multiplier,
      distance = multiplier * ElapsedTime_s
    ) %>%
    dplyr::select(ID, ElapsedTime_s, speed, distance, 2:15, 17) #removing other Ba + unwanted units
  
  data1 <- bind_rows(sheet, data1)
}

#also remove rows where data = NA
emptyrow <- apply(data1[, 7:19], 1, function(row) all(is.na(row)))
data2 <- data1[!emptyrow, ]

# TO ADDRESS negative concentrations (below limit of detection)
# Explore lowest negative values to determine if there is a problem with the mass spec
# If not, convert all negative values to 0
data2[, 7:19] <- lapply(data2[, 7:19], function(x) {
  x[x < 0] <- 0
  return(x)
})

Ca_valt <- 388000/40.08 #if using generic for teeth
Ca_valv <- 455000/40.08 #if using generic for verts
#swap w/ (Ca_ppm_m44/40.08) for both kinds

# create new data columns with PPM to Element:Ca conversion 
data_rat <- data2 %>% 
  mutate(Li_Ca = Li_ppm_m7/6.94/(Ca_ppm_m44/40.08)*1000000,
         B_Ca = B_ppm_m11/10.81/(Ca_ppm_m44/40.08)*1000000,
         Na_Ca = Na_ppm_m23/22.99/(Ca_ppm_m44/40.08)*1000000,
         Mg_Ca = Mg_ppm_m24/24.31/(Ca_ppm_m44/40.08)*1000000,
         Mn_Ca = Mn_ppm_m55/54.94/(Ca_ppm_m44/40.08)*1000000,
         Ni_Ca = Ni_ppm_m60/58.69/(Ca_ppm_m44/40.08)*1000000,
         Cu_Ca = Cu_ppm_m63/63.55/(Ca_ppm_m44/40.08)*1000000,
         Zn_Ca = Zn_ppm_m66/65.38/(Ca_ppm_m44/40.08)*1000000,
         Sr_Ca = Sr_ppm_m88/87.62/(Ca_ppm_m44/40.08)*1000000,
         In_Ca = In_ppm_m115/114.82/(Ca_ppm_m44/40.08)*1000000,
         Ba_Ca = Ba_ppm_m137/137.33/(Ca_ppm_m44/40.08)*1000000,
         Pb_Ca = Pb_ppm_m208/207.2/(Ca_ppm_m44/40.08)*1000000) %>% 
  dplyr::select(-c(7:19))

write_xlsx(data_rat, path = "data/output/profiles.xlsx" ) # keep clean version

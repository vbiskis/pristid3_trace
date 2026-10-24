# Details ----
#' dtw_prep.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 8 Jul 2026
#' Content: 
#' dynamic time warping strategy
#' -----------

source('helpers/functions/fxn_dtw.R')
source('helpers/functions/fxn_profile.R')
source('helpers/functions/fxn_mavs.R')
source('helpers/microHelp.R')
transsplit <- readRDS('data/output/transsplit.RDS')
names(transsplit)

plot_els <- c('Sr_Ca', 'Ba_Ca', 'Mg_Ca', 'Li_Ca')
elsg <- els[!els %in% c("Cu_Ca", "Pb_Ca", 'Ni_Ca')]
SMA001TB <- subset(flag_holes(transsplit[[54]], trigger = 5, release = 2.5), !.hole)
SMA001T <- subset(flag_holes(transsplit[[55]], trigger = 3, release = 1.5), !.hole)
SMA001VB <- transsplit[[56]]
SMA001V <- transsplit[[57]] 
SMA002TB <- subset(flag_holes(transsplit[[58]], trigger = 3, release = 1.5), !.hole)
SMA002T <- subset(flag_holes(transsplit[[59]], trigger = 3, release = 1.5), !.hole)
SMA002V <- transsplit[[60]] 

# check other ratios!----
plot_els_ba <- c('Sr_Ba', 'Mg_Ba', 'Li_Ba', 'Mn_Ba', 'Zn_Ba')
plot_els_mg <- c('Sr_Mg', 'Ba_Mg', 'Li_Mg', 'Mn_Mg', 'Zn_Mg')
denoms   <- c('Ba', 'Mg')
plot_els_xc <- c('Sr', 'Ba', 'Mg', 'Li', 'Mn', 'Zn') # feels inefficient w/e

dfs <- list(SMA001TB = SMA001TB, SMA001T = SMA001T, SMA001VB = SMA001VB, SMA001V = SMA001V,
            SMA002TB = SMA002TB, SMA002T = SMA002T, SMA002V = SMA002V)
dfs <- lapply(dfs, add_alt, elements = plot_els_xc, denoms = denoms)
list2env(dfs, envir = .GlobalEnv)

# apply mavs----
allmavs <- Reduce(union, list(elsg, plot_els_ba, plot_els_mg, 'In_Ca'))

SMA001TB <- add_mavs(SMA001TB, ratio_columns = allmavs, window_size = 11) # 43 um for all
SMA001T <- add_mavs(SMA001T, ratio_columns = allmavs, window_size = 11)
SMA001VB <- add_mavs(SMA001VB, ratio_columns = allmavs, window_size = 21)
SMA001V <- add_mavs(SMA001V, ratio_columns = allmavs, window_size = 21)

SMA002TB <- add_mavs(SMA002TB, ratio_columns = allmavs, window_size = 11)
SMA002T <- add_mavs(SMA002T, ratio_columns = allmavs, window_size = 11)
SMA002V <- add_mavs(SMA002V, ratio_columns = allmavs, window_size = 21)

names(SMA001TB) 

# filters----
SMA001T_filt <- SMA001T %>% filter(true_dist >= 650, true_dist <= 55583) %>% na.omit()
make_plot(SMA001T_filt, 'Sr_Ca', color = "brown", # had to clip it 3000 extra due to contam
          'Sr_Ca', span = 0.05, show_In = TRUE)

SMA001TB_filt <- SMA001TB %>% filter(true_dist >= 5000, # remove na (overcautious but to be safe)
                                     true_dist <= 55402) %>% na.omit()
make_plot(SMA001TB_filt, 'Sr_Ca', color = "brown", 
          'Sr_Ca', span = 0.05, show_In = TRUE)

SMA001V_filt <- SMA001V %>% filter(true_dist >= 710, true_dist <= 10659) %>% na.omit()
make_plot(SMA001V_filt, 'Sr_Ca', color = "brown", 
          'Sr_Ca', span = 0.0125, show_In = TRUE, In_style = "points")

SMA001VB_filt <- SMA001VB %>% filter(true_dist >= 745, true_dist <= 9968) %>% na.omit() # from measurements
make_plot(SMA001VB_filt, 'Sr_Ca', color = "brown", 
          'Sr_Ca', span = 0.0125, show_In = TRUE, In_style = "points")

SMA002TB_filt <- SMA002TB %>% filter(true_dist <= 49907) %>% na.omit() # yuck - toss
make_plot(SMA002TB_filt, 'Sr_Ca', color = "brown", 
          'Sr_Ca', span = 0.05, show_In = TRUE, In_style = "points") 

SMA002T_filt <- SMA002T %>% na.omit() # he's great tho!
make_plot(SMA002T_filt, 'Sr_Ca', color = "brown", 
          'Sr_Ca', span = 0.05, show_In = TRUE, In_style = "points")

SMA002V_filt <- SMA002V %>% filter(true_dist >= 334, true_dist <= 7952) %>% na.omit() # from measurements
make_plot(SMA002V_filt, 'Sr_Ca', color = "brown", 
          'Sr_Ca', span = 0.0125, show_In = TRUE, In_style = "points")

# and lets save these for later use in 6b-e----


# get the name for mavs cols!
names(SMA001TB_filt) # they arent nearly as good - lose too much detail with smoothing
mavs1 <- c('Sr_Ca_mav', 'Ba_Ca_mav', 'Li_Ca_mav', 'Mg_Ca_mav')
mavs2 <- c('Sr_Ca_mav', 'Ba_Ca_mav', 'Li_Ca_mav', 'Mg_Ca_mav', 'Mn_Ca_mav', 'Zn_Ca_mav')
mavsb <- names(SMA001VB_filt)[43:47]
mavsm <- names(SMA001VB_filt)[48:51]

# structure comparison----
# but be systematic

spans <- c(0.1, 0.075, 0.05, 0.025, 0.01, 0.0075, 0.005) #0.1 is high but
 # and 0.0025 beyond plotting capability
ns <- c(500, 1000, 2000)

## tooth 1----
tooth_comp <- align_grid(SMA001TB_filt, SMA001T_filt, mavs2, spans, ns)
# quick function for summarising these
result <- select_spans(tooth_comp, n_fixed = 2000, r2_threshold = 0.9) # yes!
span_vec <- setNames(result$chosen$span, result$chosen$element) # for ca ratios

# median looks to be 0.01

### test it now----
diff(range(SMA001TB_filt$true_dist, na.rm = TRUE)) # for comparison
result$consensus_end - result$consensus_start # better!

align_els(SMA001TB_filt, SMA001T_filt, elements = names(span_vec), 
          spanT = 0.01, spanV = 0.01, n = 2000) # all start around 20! nice work

finalt <- align_multi(SMA001T_filt, SMA001TB_filt, 
                      elements = mavs1, 
                      open_begin = TRUE, open_end = TRUE,
                      spanT = 0.05, spanV = 0.05, 
                      n = 2000) # and this matches that number well too!

range(SMA001T_filt$true_dist) 
range(SMA001TB_filt$true_dist) 

## vert----
# but most importantly - our verts
vert_comp <- align_grid(SMA001V_filt, SMA001VB_filt, mavs2, spans, ns)
result <- select_spans(vert_comp, n_fixed = 2000, r2_threshold = 0.9) # awesome!
span_vecv <- setNames(result$chosen$span, result$chosen$element) # oh wow

# how much wigglier are we talking 
span_vecv/span_vec2
# that actually makes a lot of sense. 

align_els(SMA001VB_filt, SMA001V_filt, elements = names(span_vec), 
          spanT = 0.001, spanV = 0.001, n = 2000) 
# plots look EXCELLENT!

finalv <- align_multi(SMA001VB_filt, SMA001V_filt, 
                      elements = mavs1, 
                      open_begin = TRUE, open_end = TRUE,
                      spanT = 0.01, spanV = 0.01, 
                      n = 2000)

# straight up... wow

range(SMA001VB_filt$true_dist)
range(SMA002T_filt$true_dist) # 4x; should match spans

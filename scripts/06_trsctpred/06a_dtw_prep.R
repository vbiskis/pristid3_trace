# Details ----
#' dtw.R
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
SMA001TB <- transsplit[[54]] %>% filter(In_Ca < 20) # get em again clean in case
SMA001T <- transsplit[[55]]  %>% filter(In_Ca < 20)
SMA001VB <- transsplit[[56]]
SMA001V <- transsplit[[57]] 
SMA002TB <- transsplit[[58]] %>% filter(In_Ca < 20)
SMA002T <- transsplit[[59]] %>% filter(In_Ca < 20)
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
SMA001T_filt <- SMA001T %>% 
  filter(true_dist >= 220, true_dist <= 54520) %>% na.omit()
make_plot(SMA001T_filt, 'Ba_Ca', color = "brown", # just a check
          'Ba_Ca', span = 0.025, show_In = TRUE, In_style = "points")

SMA001TB_filt <- SMA001TB %>% 
  filter(true_dist >= 7720, true_dist <= 51819) %>% na.omit()
make_plot(SMA001TB_filt, 'Ba_Ca', color = "brown", 
          'Ba_Ca', span = 0.025, show_In = TRUE, In_style = "points")

SMA001V_filt <- SMA001V %>% filter(true_dist >= 382, true_dist <= 10286) %>% na.omit()
make_plot(SMA001V_filt, 'Mn_Ca', color = "brown", 
          'Mn_Ca', span = 0.025, show_In = TRUE, In_style = "points")

SMA001VB_filt <- SMA001VB %>% filter(true_dist <= 9843.7) %>% na.omit() # from measurements
make_plot(SMA001VB_filt, 'Mn_Ca', color = "brown", 
          'Mn_Ca', span = 0.025, show_In = TRUE, In_style = "points")

SMA002TB_filt <- SMA002TB %>% 
  filter(true_dist >= 12810, true_dist <= 45685) %>% na.omit() # gets real yuck towards end
make_plot(SMA002TB_filt, 'Mn_Ca', color = "brown", # just a check
          'Mn_Ca', span = 0.025, show_In = TRUE, In_style = "points") # honestly not usable tho

SMA002T_filt <- SMA002T %>% 
  filter(true_dist >= 180, true_dist <= 43000) %>% na.omit() # he's great tho!
make_plot(SMA002T_filt, 'Mn_Ca', color = "brown", # just a check
          'Mn_Ca', span = 0.025, show_In = TRUE, In_style = "points")

SMA002V_filt <- SMA002V %>% 
  filter(true_dist <= 7952) %>% na.omit() # from measurements
make_plot(SMA002V_filt, 'Mn_Ca', color = "brown", # just a check
          'Mn_Ca', span = 0.025, show_In = TRUE, In_style = "points")

# and get the name for mavs cols!
names(SMA001TB_filt) # they arent nearly as good - lose too much detail with smoothing
mavs1 <- names(SMA001VB_filt)[c(36,38, 41:42)]
mavs2 <- names(SMA001VB_filt)[36:42]
mavsb <- names(SMA001VB_filt)[43:47]
mavsm <- names(SMA001VB_filt)[48:52]

# structure comparison----
## but be systematic----

spans <- c(0.1, 0.05, 0.025, 0.01, 0.0075, 0.005, 0.0025, 0.001)
ns    <- c(500, 1000, 2000)

tooth_comp <- align_grid(SMA001TB_filt, SMA001T_filt, mavs2, spans, ns)
 # enter function for summarising these
result <- select_spans(tooth_comp, n_fixed = 2000, r2_threshold = 0.9) # yes!
span_vec <- setNames(result$chosen$span, result$chosen$element) # for ca ratios
# lets test it now----

range(SMA001TB_filt$true_dist, na.rm = TRUE)
diff(range(SMA001TB_filt$true_dist, na.rm = TRUE)) # for comparison
result$consensus_end - result$consensus_start # better!

align_els(SMA001TB_filt, SMA001T_filt, elements = names(span_vec), 
          spanT = span_vec*2, spanV = span_vec*2, n = 2000) # great

finalt <- align_multi(SMA001TB_filt, SMA001T_filt, 
                      elements = mavs1, 
                      open_begin = TRUE, open_end = TRUE,
                      spanT = 0.005, spanV = span_vec[mavs1], 
                      n = 2000) 

# mm fits the whole thing not great... likely bc of too many elements here 
range(SMA001T_filt$true_dist) #cut ~ 4.4k
range(SMA001TB_filt$true_dist) #cut 4k as well; so end SHOULD be close
# so i actually like the original

# vert----
vert_comp <- align_grid(SMA001VB_filt, SMA001V_filt, mavs2, spans, ns)
result <- select_spans(vert_comp, n_fixed = 2000, r2_threshold = 0.9) # awesome!
span_vecv <- setNames(result$chosen$span, result$chosen$element) # oh wow

# how much wigglier are we talking 
span_vecv/span_vec
# hmm a fair bit

align_els(SMA001VB_filt, SMA001V_filt, elements = names(span_vec), 
          spanT = span_vecv, spanV = span_vecv, n = 2000) 
# plots look EXCELLENT!

finalv <- align_multi(SMA001VB_filt, SMA001V_filt, 
                      elements = mavs1, 
                      open_begin = TRUE, open_end = TRUE,
                      spanT = span_vecv[mavs1]*2, spanV = span_vecv[mavs1]*2, 
                      n = 6000)

# straight up... wow
# so multipling the span 
range(SMA001V_filt$true_dist) # ah! seems to make it better, closer to indivdual alignment
span_verts <-  span_vecv*2

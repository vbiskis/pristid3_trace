# Details ----
#' dtw_run.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 28 Aug 2026
#' Content: 
#' dynamic time warping II
#' application of functions developed previously
#' setting the models built with parameters developed in dtw1
#' -----------

source('helpers/functions/fxn_dtw.R')
source('helpers/functions/fxn_mavs.R')
source('helpers/microHelp.R')

#SMA001----
align_els(SMA001T_filt %>% filter(true_dist <= 43000), 
          SMA001V_filt, # Ba gets confused
          elements = names(span_vec), 
          spanT = 0.08, spanV = 0.02, n = 2000) 

align_els(SMA001TB_filt %>% filter(true_dist <= 46000), 
          SMA001VB_filt, 
          elements = names(span_vec), 
          spanT = 0.08, spanV = 0.02, n = 2000) 

# similar...

## multi----
mult1b2 <- align_multi(SMA001TB_filt %>% filter(true_dist <= 46000), 
                       SMA001VB_filt, 
                     elements = c('Ba_Ca_mav', 'Sr_Ca_mav'),
                     spanT = 0.08, spanV = 0.02,
                     n = 6000, open_end = TRUE, open_begin = TRUE) 

mult1b2.1 <- align_multi(SMA001TB_filt %>% filter(true_dist <= 46000), 
                       SMA001VB_filt, 
                       elements = c('Mg_Ca_mav', 'Sr_Ca_mav'),
                       spanT = 0.08, spanV = 0.02,
                       n = 6000, open_end = TRUE, open_begin = TRUE) 

mult1b3 <- align_multi(SMA001TB_filt %>% filter(true_dist <= 46000), 
                         SMA001VB_filt, 
                         elements = c('Mg_Ca_mav', 'Sr_Ca_mav', "Mn_Ca_mav"),
                         spanT = 0.08, spanV = 0.02,
                         n = 6000, open_end = TRUE, open_begin = TRUE) 

mult1b4 <- align_multi(SMA001TB_filt %>% filter(true_dist <= 46000), 
                       SMA001VB_filt, 
                       elements = mavs1,
                       spanT = 0.08, spanV = 0.02,
                       n = 6000, open_end = TRUE, open_begin = TRUE) 

mult1b6 <- align_multi(SMA001TB_filt %>% filter(true_dist <= 46000),
                       SMA001VB_filt, 
                       elements = mavs2,
                       spanT = 0.08, spanV = 0.02,
                       n = 6000, open_end = TRUE, open_begin = TRUE) 

# ahhh that's interesting!

mult1h2 <- align_multi(SMA001T_filt, SMA001V_filt, 
                      elements = c('Ba_Ca_mav', 'Sr_Ca_mav'),
                      spanT = 0.08, spanV = 0.02,
                      n = 6000, open_begin = TRUE, # lockV_start = 3266.09,
                      open_end = TRUE) # same result for both verts

mult1h2.1 <- align_multi(SMA001T_filt, SMA001V_filt, 
                         elements = c('Mg_Ca_mav', 'Sr_Ca_mav'),
                         spanT = 0.08, spanV = 0.02,
                         n = 6000, open_end = TRUE, open_begin = TRUE) 

mult1h3 <- align_multi(SMA001T_filt, SMA001V_filt, 
                       elements = c('Mg_Ca_mav', 'Sr_Ca_mav', "Mn_Ca_mav"),
                       spanT = 0.08, spanV = 0.02,
                       n = 6000, open_end = TRUE, open_begin = TRUE) 

mult1h4 <- align_multi(SMA001T_filt, SMA001V_filt,
                      elements = mavs1,
                      # spanT = 0.025, spanV = 0.005,
                      spanT = 0.08, spanV = 0.02,
                      n = 6000, open_begin = TRUE, 
                      open_end = TRUE) # but this one suggests longer

mult1h6 <- align_multi(SMA001T_filt, SMA001V_filt,
                      elements = mavs2,
                      # spanT = 0.025, spanV = 0.005,
                      spanT = 0.08, spanV = 0.02,
                      n = 6000, open_begin = TRUE, # lockV_start = 3266.09,
                      open_end = TRUE) 

#SMA002----
SMA002TH_out$plots$Ca #(never gets to the base, tooth incomplete)

align_els(SMA002T_filt, SMA002V_filt, elements = mavs2,
          spanT = 0.08, spanV = 0.02, n = 4000)

mult2 <- align_multi(SMA002T_filt, SMA002V_filt,
                       elements = c('Ba_Ca_mav', 'Sr_Ca_mav'),
                       spanT = 0.08, spanV = 0.02,
                       n = 6000, open_begin = TRUE, # lockV_start = 3266.09,
                       open_end = TRUE) #hmm! 
# oh! it's calling it right before (note tooth incomplete)

mult2.1 <- align_multi(SMA002T_filt, SMA002V_filt,
                       elements = c('Mg_Ca_mav', 'Sr_Ca_mav'),
                       spanT = 0.08, spanV = 0.02,
                       n = 6000, open_begin = TRUE, # lockV_start = 3266.09,
                       open_end = TRUE) 

mult3 <- align_multi(SMA002T_filt, SMA002V_filt,
                       elements = c('Mg_Ca_mav', 'Sr_Ca_mav', 'Mn_Ca_mav'),
                       spanT = 0.08, spanV = 0.02,
                       n = 6000, open_begin = TRUE, # lockV_start = 3266.09,
                       open_end = TRUE) 

mult4 <- align_multi(SMA002T_filt, SMA002V_filt,
            elements =  mavs1,
            spanT = 0.08, spanV = 0.02,
            n = 6000,  open_begin = TRUE,
            open_end = TRUE)

mult6 <- align_multi(SMA002T_filt, SMA002V_filt,
                       elements =  mavs2,
                       spanT = 0.04, spanV = 0.01,
                       n = 6000,  open_begin = TRUE,
                       # open_begin = FALSE, lockV_start = 3800,
                       open_end = TRUE)

# oh woah these all look like they match - not long though and they are def scrunched.
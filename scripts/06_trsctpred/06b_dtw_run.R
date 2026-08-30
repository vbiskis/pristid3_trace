# Details ----
#' dtw.R
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
transsplit <- readRDS('data/output/transsplit.RDS')
names(transsplit)
range(SMA001V_filt$true_dist)
#TB too short; not much use in using.

#SMA001----
align_els(SMA001T_filt, SMA001VB_filt %>% 
            filter(true_dist >= 1300), # option to
          elements = names(span_verts),
          spanT = span_vecv, spanV = span_vec,
          n = 2000) 

align_els(SMA001T_filt, SMA001V_filt %>% 
            filter(true_dist >= 1980), # force remove area prior to band 1
          elements = names(span_vec), 
          spanT = 0.0025, spanV = 0.0100, n = 2000) # and then it behaves!

# can also try to force to start at location identified by diff elements?

## multi----
final <- align_multi(SMA001TB_filt, SMA001VB_filt %>% filter(true_dist >= 1300), 
                     # need to remove; keeps seeing patterns that arent there
                     elements = mavs1,
                     # spanT = 0.075, spanV = 0.025,
                     spanT = span_vecv, spanV = span_vec,
                     n = 6000, open_end = TRUE, open_begin = TRUE) 

# that actually looks cool

final2 <- align_multi(SMA001T_filt, SMA001V_filt%>% filter(true_dist >= 1980),
                      elements = mavs1,
                      spanT = 0.005, spanV = 0.01, 
                      n = 6000, open_begin = TRUE, # lockV_start = 3266.09,
                      open_end = TRUE) #hmm! 

finalm <- align_multi(SMA001T_filt, SMA001V_filt %>% filter(true_dist >= 1980),
                      elements = mavsm,
                      spanT = 0.025, spanV = 0.01, 
                      n = 6000, open_begin = TRUE, # lockV_start = 3266.09,
                      open_end = TRUE) #hmm! 

#SMA002----

align_els(SMA002T_filt, SMA002V_filt, 
          elements = names(span_vec),
          spanT = span_vec, spanV = span_vecv, n = 2000)

finalsma2 <- align_multi(SMA002T_filt, SMA002V_filt,
            elements =  mavs1,
            spanT = 0.005, spanV = 0.01, 
            n = 6000,  open_begin = TRUE,
            # open_begin = FALSE, lockV_start = 3800,
            open_end = TRUE)

# eureka! 
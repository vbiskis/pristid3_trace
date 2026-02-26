# Details ----
#' dtw.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 8 Jul 2025
#' Content: 
#' dynamic time warping strategy
#' -----------

source('helpers/functions/fxn_dtw.R')
source('helpers/microHelp.R')
transsplit <- readRDS('data/output/transsplit.RDS')

plot_els <- c('Sr_Ca', 'Ba_Ca', 'Mg_Ca', 'Li_Ca')

#SMA001----
SMA001T <- transsplit[[55]]
SMA001T <- SMA001T %>% 
  filter(true_dist >= true_dist[which(In_Ca < 5)][1],
         true_dist <= 55402)

SMA001V <- transsplit[[57]]
SMA001V <- SMA001V %>% 
  filter(true_dist >= true_dist[which(In_Ca < 10)][1],
         true_dist <= true_dist[tail(which(In_Ca < 1), 1)])

#test it!----
align_els(SMA001T, SMA001V, plot_els, pattern = asymmetric,
          spanT = 0.015, spanV = 0.03) # WE GOT IT - just had to norm and reframe

elsg <- els[!els %in% c("Cu_Ca", "Pb_Ca", 'Ni_Ca')]

align_multi(SMA001T, SMA001V, elsg, ws = 200, n = 500,
            spanT = 0.025, spanV = 0.05) # funny one, not sure if woth keeping

## bleach----
SMA001TB <- transsplit[[54]]
SMA001TB <- SMA001TB %>% 
  filter(true_dist >= 10000, 
         true_dist <= 54279)

SMA001VB <- transsplit[[56]]
SMA001VB <- SMA001VB %>% 
  filter(true_dist >= true_dist[which(In_Ca < 10)][1],
         true_dist <= true_dist[tail(which(In_Ca < 1), 1)])

align_els(SMA001TB, SMA001VB, plot_els, pattern = asymmetric,
          spanT = 0.015, spanV = 0.03)

align_multi(SMA001TB, SMA001VB, elsg, ws = 200, n = 500,
            spanT = 0.025, spanV = 0.05)

#SMA002----
SMA002T <- transsplit[[59]] 
SMA002T <- SMA002T %>% 
  filter(true_dist >= 130, 
         true_dist <= max(true_dist[which(In_Ca < 1)])) 

SMA002V <- transsplit[[60]]
SMA002V <- SMA002V %>% 
  filter(true_dist >= true_dist[which(In_Ca < 10)][1], 
         true_dist <= 7900)

align_els(SMA002T, SMA002V, plot_els, pattern = asymmetricP0,
          spanT = 0.02, spanV = 0.06)

align_multi(SMA002T, SMA002V, elsg)

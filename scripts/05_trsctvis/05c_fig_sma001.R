# Details ----
#' 05c_fig_sma001.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 16 Jul 2024
#' Content: 
#' vert and teeth pics for aquarium samples, sma001
#' this one's manual, doesn't use the loop, can adjust
#' -----------

source('helpers/microHelp.R')
set_theme(mytheme)
source('helpers/functions/fxn_mavs.R')
source('helpers/functions/fxn_profile.R')
transsplit <- readRDS('data/output/transsplit.RDS')
names(transsplit)[54:60] # and I can easily check this bc I was smart earlier and named things
# want 54-57 for SMA001

# water----
## tooth----
SMA001T <- transsplit[[55]] %>% 
  filter(true_dist >= 165, true_dist <= 55798) %>% na.omit()
# 165 um from tip (start) and 169 from end off the sample

# lets remove In spikes
h1 <- plot_holes(SMA001T, trigger =  3, release = 1.5,
                 elements = c("Sr_Ca", "Ba_Ca", "Li_Ca", "Mg_Ca"))

SMA001T <- h1$data[!h1$data$.hole, ]

# sanity check:
get_spikes(id = SMA001T) %>% 
  round(digits = 2)

## holes----
find_gaps(SMA001T$true_dist) # holes > 50 um (normal spacing ~8.6 um)

# remind myself... what is a window size of 21 (10 pts each side) in um:
samp_intt <- mean(diff(sort(SMA001T$true_dist)), na.rm = TRUE)
samp_intt*5 # ws(21) = 43 um
samp_intt*10 # ws(21) = 86 um
samp_intt*12 # ws(25) = 104 um

SMA001T_out <- profs(SMA001T, el_all, "tooth", # alt_denom = "Ba", # data types incl.
                      span = 0.05, window_size = 11, # stats
                      # img_path = "data/imgs/SMA001/SMA001HT/SMA001-H-sm.jpg", # cool visual
                      show_In = TRUE, In_style = "rect", rect_threshold = 0.75) # contamination check

SMA001T_out$plots$Ca

sma01a <- SMA001T_out$plots$Sr + theme(axis.title.x = element_blank())
sma01b <- SMA001T_out$plots$Ba + theme(axis.title.x = element_blank())
sma01c <- SMA001T_out$plots$Li + theme(axis.title.x = element_blank())
sma01d <- SMA001T_out$plots$Mg + theme(axis.title.x = element_blank())
sma01e <- SMA001T_out$plots$Mn 

## vert----
SMA001V <- transsplit[[57]] %>% 
  filter(# true_dist >= 707, 
    true_dist <= 10659.44) %>% na.omit()

get_spikes(id = SMA001V) %>% 
  # filter(max_In_Ca >= 4) %>% 
  round(digits = 2) # no holes, just In at start

SMA001V_out <- profs(SMA001V, el_all, "vert", # alt_denom = "Ba", 
                      # img_path = "data/imgs/SMA001/SMA001HV/SMA001-H-bw-crop.jpg",
                      span = 0.025, window_size = 31, # thin = 2,
                      show_In = TRUE, In_style = "points")

samp_intv <- mean(diff(sort(SMA001V$true_dist)), na.rm = TRUE)
samp_intv*10 # ws(21) = 29 um! much less, noisier
samp_intv*15 # ws(31) = 43, matching teeth
samp_intv*20 # ws(21) = 60 um
samp_intv*35 # ws(71) = 101 um

SMA001V_out$plots$Ca

sma01va <- SMA001V_out$plots$Sr + theme(axis.title.x = element_blank())
sma01vb <- SMA001V_out$plots$Ba + theme(axis.title.x = element_blank())
sma01vc <- SMA001V_out$plots$Li + theme(axis.title.x = element_blank())
sma01vd <- SMA001V_out$plots$Mg + theme(axis.title.x = element_blank())
sma01ve <- SMA001V_out$plots$Mn

recta <- ggplot2::annotate("rect", xmin = 6100, xmax = 10800, 
                          ymin = 1650, ymax = 2750,
                              fill = "lightgrey", alpha = 0.3) 
sma01va$layers <- c(recta, sma01va$layers)

rectb <- ggplot2::annotate("rect", xmin = 6100, xmax = 10800, 
                           ymin = 5, ymax = 42,
                           fill = "lightgrey", alpha = 0.3) 
sma01vb$layers <- c(rectb, sma01vb$layers)

rectc <- ggplot2::annotate("rect", xmin = 6100, xmax = 10800, 
                              ymin = 15, ymax = 29,
                              fill = "lightgrey", alpha = 0.3) 
sma01vc$layers <- c(rectc, sma01vc$layers)

rectd <- ggplot2::annotate("rect", xmin = 6100, xmax = 10800,
                              ymin = 15450, ymax = 21000,
                   fill = "lightgrey", alpha = 0.3) 
sma01vd$layers <- c(rectd, sma01vd$layers)

sma001a <- (sma01a + sma01va) / (sma01b + sma01vb) / (sma01c + sma01vc) / (sma01d + sma01vd) + 
  plot_layout(widths = c(1, 1)) & theme(plot.margin = margin(5, 5, 5, 5))

ggsave(
  "SMA001VTH.png",
  plot = sma001a,
  device = NULL,
  path = 'figures/fig6',
  scale = 1,
  width = 8,
  height = 8,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
)

# bleach ----
## tooth----
SMA001TB <- transsplit[[54]]  %>% 
  filter(true_dist >= 164, 
         true_dist <= 55369)

b1 <- plot_holes(SMA001TB, trigger =  3, release = 1.5,
                 elements = c("Sr_Ca", "Ba_Ca", "Li_Ca", "Mg_Ca"))
SMA001TB <- b1$data[!b1$data$.hole, ]

get_spikes(SMA001TB)

## holes----
find_gaps(SMA001TB$true_dist) # holes > 50 um (normal spacing ~8.6 um)

SMA001TB_out <- profs(SMA001TB, el_all, "tooth", 
                       # alt_denom = "Ba",
                       span = 0.05, window_size = 11, 
                       show_In = FALSE)

SMA001TB_out$plots$Ca

sma01a2 <- SMA001TB_out$plots$Sr + theme(axis.title.x = element_blank())
sma01b2 <- SMA001TB_out$plots$Ba + theme(axis.title.x = element_blank())
sma01c2 <- SMA001TB_out$plots$Li + theme(axis.title.x = element_blank())
sma01d2 <- SMA001TB_out$plots$Mg + theme(axis.title.x = element_blank())
sma01e2 <- SMA001TB_out$plots$Mn 

SMA001TB_out$plots$Ca

## vert----
SMA001VB <- transsplit[[56]] %>% 
  filter(true_dist> 745, true_dist <= 9968.2)

get_spikes(id = SMA001VB) %>% 
  # filter(max_In_Ca >= 2) %>% 
  round(digits = 2)

SMA001VB_out <- profs(SMA001VB, el_all, "vert", # alt_denom = "Ba",
                       span = 0.025, window_size = 31, 
                       show_In = TRUE, In_style = "points") # wow that is cool
SMA001VB_out$plots$Ca

sma01va2 <- SMA001VB_out$plots$Sr + theme(axis.title.x = element_blank())
sma01vb2 <- SMA001VB_out$plots$Ba + theme(axis.title.x = element_blank())
sma01vc2 <- SMA001VB_out$plots$Li + theme(axis.title.x = element_blank())
sma01vd2 <- SMA001VB_out$plots$Mg + theme(axis.title.x = element_blank())
sma01ve2 <- SMA001VB_out$plots$Mn

recta <- ggplot2::annotate("rect", xmin = 6100, xmax = 10000, 
                           ymin = 1650, ymax = 2650,
                           fill = "lightgrey", alpha = 0.3) 
sma01va2$layers <- c(recta, sma01va2$layers)

rectb <- ggplot2::annotate("rect", xmin = 6100, xmax = 10000, 
                           ymin = 5, ymax = 62,
                           fill = "lightgrey", alpha = 0.3) 
sma01vb2$layers <- c(rectb, sma01vb2$layers)

rectc <- ggplot2::annotate("rect", xmin = 6100, xmax = 10000, 
                           ymin = 15, ymax = 28,
                           fill = "lightgrey", alpha = 0.3) 
sma01vc2$layers <- c(rectc, sma01vc2$layers)

rectd <- ggplot2::annotate("rect", xmin = 6100, xmax = 10000,
                           ymin = 15450, ymax = 20200,
                           fill = "lightgrey", alpha = 0.3) 
sma01vd2$layers <- c(rectd, sma01vd2$layers)

sma001b <- (sma01a2 + sma01va2) / (sma01b2 + sma01vb2) / (sma01c2 + sma01vc2) / (sma01d2 + sma01vd2) & 
  theme(plot.margin = margin(5, 5, 5, 5)) #thats ... identical

ggsave(
  "SMA001VTB.png",
  plot = sma001b,
  device = NULL,
  path = 'figures/fig6',
  scale = 1,
  width = 8,
  height = 8,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
)

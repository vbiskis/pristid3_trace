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
sapply(transsplit, function(df) unique(df$ID))

# want 54-57 for SMA001
#54 tooth bleach
#55 tooth water
#56 vert bleach 
#57 vert water

# water----
## tooth----
SMA001T <- transsplit[[55]]
SMA001T <- SMA001T %>% 
  filter(true_dist >= true_dist[which(In_Ca < 5)][1],
         # true_dist <= true_dist[tail(which(In_Ca < 5), 1)],
         true_dist <= 55002) # manually removing the end due to contamination

SMA001T_filt <- profs(SMA001T, el_all, "tooth", span = 0.03, window_size = 81,
                      # img_path = "data/transects/imgs/SMA001/SMA001-H-sm.jpg", #seems to be half, so double span
                      show_In = FALSE)

sma01a <- a + theme(axis.title.x = element_blank())
sma01b <- b + theme(axis.title.x = element_blank())
sma01c <- cc + theme(axis.title.x = element_blank())
sma01d <- d + theme(axis.title.x = element_blank())
sma01e <- e 

## vert----
SMA001V <- transsplit[[57]]
SMA001V <- SMA001V %>% 
  filter(true_dist >= true_dist[which(In_Ca < 10)][1], # not confident in the effect of In enough to go bigger
         true_dist <= true_dist[tail(which(In_Ca < 1), 1)])

SMA001V_filt <- profs(SMA001V, el_all, "vert", # need image!
                      span = 0.01, window_size = 81, show_In = FALSE)

sma01va <- a + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
sma01vb <- b + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
sma01vc <- cc + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
sma01vd <- d + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
sma01ve <- e
sma01vin <- f

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
SMA001TB <- transsplit[[54]]

SMA001TB <- SMA001TB %>% 
  filter(true_dist >= 10000, # big hole or some sh
         true_dist <= 54279)

SMA001TB_filt <- profs(SMA001TB, el_all, "tooth", 
                       span = 0.03, window_size = 81, show_In = FALSE)

sma01a2 <- a + theme(axis.title.x = element_blank())
sma01b2 <- b + theme(axis.title.x = element_blank())
sma01c2 <- cc + theme(axis.title.x = element_blank())
sma01d2 <- d + theme(axis.title.x = element_blank())
sma01e2 <- e + theme(axis.title.x = element_blank())
sma01f2 <- f + theme(axis.title.x = element_blank())

## vert----
SMA001VB <- transsplit[[56]]
SMA001VB <- SMA001VB %>% 
  filter(true_dist >= true_dist[which(In_Ca < 10)][1], # adjust to see window move (2, 5, 10)
         true_dist <= true_dist[tail(which(In_Ca < 1), 1)])

SMA001VB_filt <- profs(SMA001VB, el_all, "vert", span = 0.01, window_size = 81, show_In = FALSE)

sma01va2 <- a + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
sma01vb2 <- b + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
sma01vc2 <- cc + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
sma01vd2 <- d + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
sma01ve2 <- e 
sma01vf2 <- f

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

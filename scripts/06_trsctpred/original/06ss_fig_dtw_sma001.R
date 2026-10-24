#dual plot# Details ----
#' dtw_sma001.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 8 Jul 2025
#' Content: 
#' graphing - automated ggplot
#' # mav datasets made in the dtw script
#' -----------

source('helpers/microHelp.R')
source('helpers/functions/fxn_mavs.R')
transsplit <- readRDS('data/output/transsplit.RDS')
sapply(transsplit, function(df) unique(df$ID))

# sma001----
SMA001T <- transsplit[[55]]
SMA001T <- SMA001T %>% 
  filter(true_dist >= true_dist[which(In_Ca < 5)][1],
         true_dist <= 55402)

SMA001V <- transsplit[[57]]
SMA001V <- SMA001V %>% 
  filter(true_dist >= true_dist[which(In_Ca < 10)][1],
         true_dist <= true_dist[tail(which(In_Ca < 1), 1)])

for(element in els) {
  smooth_col <- paste0(element, "_mav")
  SMA001T[[smooth_col]] <- predict(loess(SMA001T[[element]] ~ SMA001T$true_dist, span = 0.1))
  SMA001V[[smooth_col]] <- predict(loess(SMA001V[[element]] ~ SMA001V$true_dist, span = 0.05))
}

#from script 06a
slow_range <- c(6430, 10650) # vert range
fast_range <- c(717, 55393) # tooth range

scale_factor <- diff(slow_range) / diff(fast_range)  
offset <- slow_range[1] - fast_range[1] * scale_factor

SMA001T$distance_scaled <- SMA001T$true_dist * scale_factor + offset

comp_plot <- function(tooth_data, vert_data, y_var, y_lab, 
                                 tooth_span = 0.15, vert_span = 0.05,
                                 show_xlab = FALSE){
  df1 <- data.frame(distance = tooth_data$distance_scaled,
                    value = tooth_data[[y_var]],
                    series = "T_scaled")
  
  df2 <- data.frame(distance = vert_data$true_dist,
                    value = vert_data[[paste0(y_var, "_mav")]],
                    series = "V")
  
  # scale tooth to vert range
  y_scale  <- diff(range(df2$value, na.rm = TRUE)) / diff(range(df1$value, na.rm = TRUE))
  y_offset <- min(df2$value, na.rm = TRUE) - min(df1$value, na.rm = TRUE) * y_scale
  df1$value <- df1$value * y_scale + y_offset
  
  p <- ggplot() +
    geom_smooth(data = df2, aes(distance, value), method = "loess", span = vert_span, color = "darkgrey") +
    geom_smooth(data = df1, aes(distance, value), method = "loess", span = tooth_span, 
                color = "brown", se = TRUE) +
    scale_y_continuous(name = y_lab) +
    labs(x = if(show_xlab) "Distance (um)" else NULL) +
    theme_bw() +
    theme(axis.title.x = if(show_xlab) element_text() else element_blank())
  
  p
}

# then just:
sr <- comp_plot(SMA001T, SMA001V, "Sr_Ca", "[Sr:Ca]")
ba <- comp_plot(SMA001T, SMA001V, "Ba_Ca", "[Ba:Ca]")
li <- comp_plot(SMA001T, SMA001V, "Li_Ca", "[Li:Ca]")
mg <- comp_plot(SMA001T, SMA001V, "Mg_Ca", "[Mg:Ca]")

SMA001 <- (sr + ba) / (li + mg)

ggsave("SMA001-dtw.png", 
       plot = SMA001, path = "figures/fig7",
       width = 6, height = 4, units = "in", dpi = 300)

# bleach----

SMA001TB <- transsplit[[54]]
SMA001TB <- SMA001TB %>% 
  filter(true_dist >= 10000, 
         true_dist <= 54279)

SMA001VB <- transsplit[[56]]
SMA001VB <- SMA001VB %>% 
  filter(true_dist >= true_dist[which(In_Ca < 10)][1],
         true_dist <= true_dist[tail(which(In_Ca < 1), 1)])

for(element in els) {
  smooth_col <- paste0(element, "_mav")
  SMA001TB[[smooth_col]] <- predict(loess(SMA001TB[[element]] ~ SMA001TB$true_dist, span = 0.1))
  SMA001VB[[smooth_col]] <- predict(loess(SMA001VB[[element]] ~ SMA001VB$true_dist, span = 0.05))
}

slow_range <- c(5838.065, 9991.241) # vert range
fast_range <- c(10005.66, 54270.67) # tooth range

scale_factor <- diff(slow_range) / diff(fast_range)  
offset <- slow_range[1] - fast_range[1] * scale_factor

SMA001TB$distance_scaled <- SMA001TB$true_dist * scale_factor + offset

sr <- comp_plot(SMA001TB, SMA001VB, "Sr_Ca", "[Sr:Ca]")
ba <- comp_plot(SMA001TB, SMA001VB, "Ba_Ca", "[Ba:Ca]")
li <- comp_plot(SMA001TB, SMA001VB, "Li_Ca", "[Li:Ca]")
mg <- comp_plot(SMA001TB, SMA001VB, "Mg_Ca", "[Mg:Ca]")

SMA001B <- (sr + ba) / (li + mg)
ggsave("SMA001b-dtw.png", 
       plot = SMA001B, path = "figures/fig7",
       width = 6, height = 4, units = "in", dpi = 300)

# sma002----

SMA002T <- transsplit[[59]] 
SMA002T <- SMA002T %>% 
  filter(true_dist >= 130, 
         true_dist <= max(true_dist[which(In_Ca < 1)])) 

SMA002V <- transsplit[[60]]
SMA002V <- SMA002V %>% 
  filter(true_dist >= true_dist[which(In_Ca < 10)][1], 
         true_dist <= 7900)

for(element in els) {
  smooth_col <- paste0(element, "_mav")
  SMA002T[[smooth_col]] <- predict(loess(SMA002T[[element]] ~ SMA002T$true_dist, span = 0.1))
  SMA002V[[smooth_col]] <- predict(loess(SMA002V[[element]] ~ SMA002V$true_dist, span = 0.05))
}

slow_range <- c(5486.686, 7520.072) # vert range
fast_range <- c(138.246, 43055.37) # tooth range

scale_factor <- diff(slow_range) / diff(fast_range)  
offset <- slow_range[1] - fast_range[1] * scale_factor

SMA002T$distance_scaled <- SMA002T$true_dist * scale_factor + offset

sr <- comp_plot(SMA002T, SMA002V, "Sr_Ca", "[Sr:Ca]")
ba <- comp_plot(SMA002T, SMA002V, "Ba_Ca", "[Ba:Ca]")
li <- comp_plot(SMA002T, SMA002V, "Li_Ca", "[Li:Ca]")
mg <- comp_plot(SMA002T, SMA002V, "Mg_Ca", "[Mg:Ca]")

SMA002 <- (sr + ba) / (li + mg)

ggsave("SMA002-dtw.png", 
       plot = SMA002, path = "figures/fig7",
       width = 6, height = 4, units = "in", dpi = 300)

# Details ----
#' dtw.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 30 Aug 2026
#' Content: 
#' dynamic time warping III
#' lets check out that warp rate
#' 
#' x3 functions - depending on the question:
#' + first one is the grouped warp rate
#' + second is warp reate per element if they are nto forced aligned
#' + 3rd is now the warp rate but lookng at elements that are forced to the same area along a tooth
#' -----------

warp_rate <- function(alignment, datT, datV, n = 2000){
  
  tooth_idx <- round(seq(1, nrow(datT), length.out = n))
  vert_idx  <- round(seq(1, nrow(datV), length.out = n))
  
  tooth_dist <- datT$true_dist[tooth_idx[alignment$index1]]
  vert_dist  <- datV$true_dist[vert_idx[alignment$index2]]
  
  avg_slope <- (vert_dist[length(vert_dist)] - vert_dist[1]) / (tooth_dist[length(tooth_dist)] - tooth_dist[1])
  predicted <- vert_dist[1] + avg_slope * (tooth_dist - tooth_dist[1])
  residual  <- vert_dist - predicted   # + = ahead of average pace, - = behind
  
  df <- data.frame(tooth_dist = tooth_dist, vert_dist = vert_dist,
                   predicted = predicted, residual = residual)
  
  insert_zero_crossings <- function(df, x_col = "tooth_dist", y_col = "residual"){
    x <- df[[x_col]]
    y <- df[[y_col]]
    
    out_x <- x[1]; out_y <- y[1]
    for(i in 2:length(x)){
      if(sign(y[i-1]) != sign(y[i]) && y[i-1] != 0 && y[i] != 0){
        frac <- y[i-1] / (y[i-1] - y[i])
        cross_x <- x[i-1] + frac * (x[i] - x[i-1])
        out_x <- c(out_x, cross_x)
        out_y <- c(out_y, 0)
      }
      out_x <- c(out_x, x[i])
      out_y <- c(out_y, y[i])
    }
    
    data.frame(x = out_x, y = out_y)
  }
  
  p1 <- ggplot(df, aes(tooth_dist, vert_dist)) +
    geom_line(aes(color = residual), linewidth = 1) +
    geom_line(aes(y = predicted), linetype = "dashed", color = "grey40") +
    scale_color_gradient2(low = "#2166ac", mid = "grey85", high = "#b2182b", midpoint = 0,
                          name = "Ahead / behind\naverage pace (um)") +
    theme_bw() +
    labs(x = "Position along tooth (um)", y = "Corresponding position on vert (um)",
         title = "Warping path vs. average deposition rate")
  
  df_clean <- insert_zero_crossings(df, "tooth_dist", "residual")
  
  p2 <- ggplot(df_clean, aes(x)) +
    geom_ribbon(aes(ymin = 0, ymax = pmax(y, 0)), fill = "#b2182b", alpha = 0.4) +
    geom_ribbon(aes(ymin = pmin(y, 0), ymax = 0), fill = "#2166ac", alpha = 0.4) +
    geom_line(aes(y = y), color = "black", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    labs(x = "Position along tooth (um)", y = "Deviation from\naverage rate (um)")
  
  p1 / p2
}

warp_rate(final, SMA001TB_filt, SMA001VB_filt %>% filter(true_dist >= 1300), n = 6000)
warp_rate(final2, SMA001T_filt, SMA001V_filt %>% filter(true_dist >= 1980), n = 6000)
warp_rate(finalm, SMA001T_filt, SMA001V_filt %>% filter(true_dist >= 1980), n = 6000)
warp_rate(finalsma2, SMA002T_filt, SMA002V_filt, n = 6000) # hmm!

warp_rate(finalt, SMA001TB_filt, SMA001T_filt, n = 2000) # just a comparison across the teeth
warp_rate(finalv, SMA001VB_filt, SMA001V_filt, n = 6000) # and lets look at verts

# individual
warps_solo <- function(datT, datV, elements, spanT, spanV, n = 2000, pattern = asymmetricP0){
  
  datT <- datT[order(datT$true_dist), ]
  datV <- datV[order(datV$true_dist), ]
  
  if(length(spanT) == 1) spanT <- setNames(rep(spanT, length(elements)), elements)
  if(length(spanV) == 1) spanV <- setNames(rep(spanV, length(elements)), elements)
  
  plots <- list()
  for(element in elements){
    smooth_col <- paste0(element, "_smooth")
    datT[[smooth_col]] <- predict(loess(datT[[element]] ~ datT$true_dist, span = spanT[[element]]))
    datV[[smooth_col]] <- predict(loess(datV[[element]] ~ datV$true_dist, span = spanV[[element]]))
    
    tSmooth <- approx(seq_along(datT[[smooth_col]]), datT[[smooth_col]], n = n)$y
    vSmooth <- approx(seq_along(datV[[smooth_col]]), datV[[smooth_col]], n = n)$y
    tSmooth <- (tSmooth - mean(tSmooth)) / sd(tSmooth)
    vSmooth <- (vSmooth - mean(vSmooth)) / sd(vSmooth)
    
    align <- dtw(tSmooth, vSmooth, keep = TRUE, step = pattern, open.end = TRUE, open.begin = TRUE)
    
    tooth_idx <- round(seq(1, nrow(datT), length.out = n))
    vert_idx  <- round(seq(1, nrow(datV), length.out = n))
    
    tooth_dist <- datT$true_dist[tooth_idx[align$index1]]
    vert_dist  <- datV$true_dist[vert_idx[align$index2]]
    
    avg_slope <- (vert_dist[length(vert_dist)] - vert_dist[1]) / (tooth_dist[length(tooth_dist)] - tooth_dist[1])
    predicted <- vert_dist[1] + avg_slope * (tooth_dist - tooth_dist[1])
    residual  <- vert_dist - predicted
    
    df <- data.frame(tooth_dist = tooth_dist, vert_dist = vert_dist, predicted = predicted, residual = residual)
    
    plots[[element]] <- ggplot(df, aes(tooth_dist, vert_dist)) +
      geom_line(aes(color = residual), linewidth = 1) +
      geom_line(aes(y = predicted), linetype = "dashed", color = "grey40") +
      scale_color_gradient2(low = "#2166ac", mid = "grey85", high = "#b2182b", midpoint = 0,
                            name = "Deviation\n(um)") +
      theme_bw() + labs(x = "Tooth (um)", y = "Vert (um)", title = element)
  }
  
  wrap_plots(plots, ncol = 2)
}

warps_solo(SMA001T_filt, SMA001V_filt %>% filter(true_dist >= 1980),
           elements = mavs1, spanT = span_vecv, spanV = span_vec, n = 6000)

warps_solo(SMA001TB_filt, SMA001VB_filt %>% filter(true_dist >= 1300), # and then li flips? hmmmm
           elements = mavs1, spanT = span_vecv, spanV = span_vec, n = 6000) 

warps_solo(SMA001T_filt, SMA001V_filt %>% filter(true_dist >= 1980), # that makes... more sense?
           elements = mavsb, spanT = 0.01, spanV = 0.05, n = 6000)

warps_solo(SMA002T_filt, SMA002V_filt, elements = mavs1, spanT = 0.005, spanV = 0.01, n = 6000)

warp_team <- function(alignment, datT, datV, elements, spanT, spanV, n = 2000){
  
  tooth_idx <- round(seq(1, nrow(datT), length.out = n))
  vert_idx  <- round(seq(1, nrow(datV), length.out = n))
  
  tooth_dist <- datT$true_dist[tooth_idx[alignment$index1]]
  vert_dist  <- datV$true_dist[vert_idx[alignment$index2]]
  
  if(length(spanT) == 1) spanT <- setNames(rep(spanT, length(elements)), elements)
  if(length(spanV) == 1) spanV <- setNames(rep(spanV, length(elements)), elements)
  
  plots <- list()
  for(element in elements){
    smooth_col <- paste0(element, "_smooth")
    datT[[smooth_col]] <- predict(loess(datT[[element]] ~ datT$true_dist, span = spanT[[element]]))
    datV[[smooth_col]] <- predict(loess(datV[[element]] ~ datV$true_dist, span = spanV[[element]]))
    
    q <- scale(datT[[smooth_col]][tooth_idx[alignment$index1]])[,1]
    r <- scale(datV[[smooth_col]][vert_idx[alignment$index2]])[,1]
    
    value_residual <- q - r   # how far tooth's value is from vert's value at each matched point
    
    df <- data.frame(tooth_dist = tooth_dist, vert_dist = vert_dist, value_residual = value_residual)
    
    plots[[element]] <- ggplot(df, aes(tooth_dist, vert_dist)) +
      geom_line(aes(color = value_residual), linewidth = 1) +
      scale_color_gradient2(low = "#2166ac", mid = "grey85", high = "#b2182b", midpoint = 0,
                            name = "Tooth\u2212Vert\n(z-score)") +
      theme_bw() + labs(x = "Position on tooth (um)", y = "Position on vert (um)", title = element)
  }
  
  wrap_plots(plots, ncol = 2)
}

warp_team(final, SMA001TB_filt, SMA001VB_filt %>% filter(true_dist >= 1300),
                            elements = mavs1, spanT = span_vecv, spanV = span_vec, n = 6000)

warp_team(final2, SMA001T_filt, SMA001V_filt %>% filter(true_dist >= 1980),
           elements = mavs1, spanT = span_vecv, spanV = span_vec, n = 6000) # yikes completely opposite li/ba results!

warp_team(finalm, SMA001T_filt, SMA001V_filt %>% filter(true_dist >= 1980),
           elements = mavsb, spanT = 0.01, spanV = 0.05, n = 6000)

warp_team(finalsma2, SMA002T_filt, SMA002V_filt, elements = mavs1, spanT = 0.005, spanV = 0.01, n = 6000)


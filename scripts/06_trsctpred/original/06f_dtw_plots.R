# Details ----
#' dtwV.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 30 Aug 2026
#' Content: 
#' dynamic time warping IV
#' output - plots
#' -----------

source('helpers/functions/fxn_gaps.R') # smoothing that stops at holes

plot_overlap <- function(alignment, datT, datV, elements, spanT, spanV, n = 2000,
                         gap_um = 50){ # holes wider than this (um): smoothed either side, not across
  
  datT <- datT[order(datT$true_dist), ]
  datV <- datV[order(datV$true_dist), ]
  
  if(length(spanT) == 1) spanT <- setNames(rep(spanT, length(elements)), elements)
  if(length(spanV) == 1) spanV <- setNames(rep(spanV, length(elements)), elements)
  
  for(element in elements) {
    smooth_col <- paste0(element, "_smooth")
    datT[[smooth_col]] <- smooth_gaps(datT$true_dist, datT[[element]], spanT[[element]], gap_um)
    datV[[smooth_col]] <- smooth_gaps(datV$true_dist, datV[[element]], spanV[[element]], Inf) # verts: never split
  }
  
  tooth_idx <- round(seq(1, nrow(datT), length.out = n))
  vert_idx  <- round(seq(1, nrow(datV), length.out = n))
  
  vert_start <- datV$true_dist[vert_idx[alignment$index2[1]]]
  vert_end   <- datV$true_dist[vert_idx[alignment$index2[length(alignment$index2)]]]
  
  plots <- list()
  for(element in elements){
    smooth_col <- paste0(element, "_smooth")
    
    df_vert <- data.frame(distance = datV$true_dist, value = datV[[smooth_col]])
    df_tooth <- data.frame(
      distance = datV$true_dist[vert_idx[alignment$index2]],
      value = datT[[smooth_col]][tooth_idx[alignment$index1]]
    )
    
    vert_region <- df_vert[df_vert$distance >= vert_start & df_vert$distance <= vert_end, ]
    y_scale <- diff(range(vert_region$value, na.rm = TRUE)) / diff(range(df_tooth$value, na.rm = TRUE))
    y_offset <- min(vert_region$value, na.rm = TRUE) - min(df_tooth$value, na.rm = TRUE) * y_scale
    df_tooth$value <- df_tooth$value * y_scale + y_offset
    
    cat(element, "y_scale:", y_scale, "y_offset:", y_offset,
        "| vert_region range:", range(vert_region$value, na.rm = TRUE),
        "| df_tooth range (post-scale):", range(df_tooth$value, na.rm = TRUE), "\n")
    
    plots[[element]] <- ggplot() +
      ggplot2::annotate("rect", xmin = vert_start, xmax = vert_end, ymin = -Inf, ymax = Inf,
               fill = "steelblue", alpha = 0.08) +
      geom_vline(xintercept = c(vert_start, vert_end), linetype = "dashed", color = "grey40") +
      geom_line(data = df_vert, aes(distance, value), color = "darkgrey") +
      geom_line(data = df_tooth, aes(distance, value), color = "brown") +
      theme_bw() +
      labs(x = "Distance along vert (um)", y = element)
  }
  
  wrap_plots(plots, ncol = 2) +
    plot_annotation(
      caption = paste(strwrap(sprintf(paste(
        "Smoothed vertebra (grey) with the smoothed tooth (brown) laid along it by the alignment, rescaled to the",
        "vertebra's range over the matched region (shaded, between dashed lines): tooth overlaps vertebra",
        "%.0f to %.0f um (%.0f um)."), vert_start, vert_end, vert_end - vert_start), 120), collapse = "\n"),
      theme = theme(plot.caption = element_text(hjust = 0))
    )
}

plot_overlap(mult1h6, SMA001T_filt, SMA001V_filt, 
             mavs2, spanT = 0.04, spanV = 0.01, n = 6000)

plot(mult1h6, type = "twoway", offset = 1)

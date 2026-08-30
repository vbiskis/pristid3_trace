# Details ----
#' dtw.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 30 Aug 2026
#' Content: 
#' dynamic time warping IV
#' output - plots
#' -----------

plot_overlap <- function(alignment, datT, datV, elements, spanT, spanV, n = 2000){
  
  datT <- datT[order(datT$true_dist), ]
  datV <- datV[order(datV$true_dist), ]
  
  if(length(spanT) == 1) spanT <- setNames(rep(spanT, length(elements)), elements)
  if(length(spanV) == 1) spanV <- setNames(rep(spanV, length(elements)), elements)
  
  for(element in elements) {
    smooth_col <- paste0(element, "_smooth")
    datT[[smooth_col]] <- predict(loess(datT[[element]] ~ datT$true_dist, span = spanT[[element]]))
    datV[[smooth_col]] <- predict(loess(datV[[element]] ~ datV$true_dist, span = spanV[[element]]))
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
      labs(x = "Distance along vert (um)", y = element, title = element)
  }
  
  wrap_plots(plots, ncol = 2) +
    plot_annotation(
      title = sprintf("Tooth overlap region: %.0f to %.0f um (span: %.0f um)",
                      vert_start, vert_end, vert_end - vert_start)
    )
}

plot_overlap(final, SMA001TB_filt, SMA001VB_filt %>% filter(true_dist >= 1300), 
             names(span_vec[mavs1]),
             spanT = span_vecv, spanV = span_vec,
             n = 6000)
plot(final, type = "twoway", off = 1)

plot_overlap(final2, SMA001T_filt, SMA001V_filt %>% filter(true_dist >= 1980), 
             names(span_vec[mavs1]),
             spanT = span_vecv, spanV = span_vec,
             n = 6000)

plot(final2, type = "twoway", off = 1) # its cool that they have the same pattern basically 
# for both teeth-vert combo

# sma002
plot_overlap(finalsma2, SMA002T_filt, SMA002V_filt, 
             names(span_vec[mavs1]),
             spanT = span_vecv, spanV = span_vec,
             n = 6000)

plot(finalsma2, type = "twoway", off = 1) # its cool that they have the same pattern basically 
# for both teeth-vert combo


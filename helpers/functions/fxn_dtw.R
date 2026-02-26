# Details ----
#' dtw.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 8 Jul 2025
#' Updated: 25 Feb 2026
#' ++ turned original code into a programming for looping easier
#' Content: 
#' dynamic time warping functions x 2
#' ++ one basic visualisation for each element
#' ++ one for clustering
#' -----------

library(dtw) # fxn 1
library(dtwclust) # fxn 2

align_els <- function(datT, datV, elements, n = 1000,
                      pattern = asymmetric, spanT = 0.05, spanV = 0.15) {
  
  # fit smooths (from my first script version)
  for(element in elements) {
    smooth_col <- paste0(element, "_smooth")
    datT[[smooth_col]] <- predict(loess(datT[[element]] ~ datT$true_dist, span = spanT))
    datV[[smooth_col]] <- predict(loess(datV[[element]] ~ datV$true_dist, span = spanV))
  } # could modify with other smoothers
  
  results <- list()
  
  for(element in elements) {
    smooth_col <- paste0(element, "_smooth")
    
    # interpolate to same length (as in Hegg 2021)
    tSmooth <- approx(seq_along(datT[[smooth_col]]), datT[[smooth_col]], n = n)$y
    vSmooth <- approx(seq_along(datV[[smooth_col]]), datV[[smooth_col]], n = n)$y
    
    # then normalize (as in Likenday 2025)
    tSmooth <- (tSmooth - mean(tSmooth)) / sd(tSmooth)
    vSmooth <- (vSmooth - mean(vSmooth)) / sd(vSmooth)
    
    alignment <- dtw(tSmooth, vSmooth,
                     keep = TRUE, step = pattern, #changeable
                     open.end = TRUE, open.begin = TRUE)
    
    datT <- datT[order(datT$true_dist), ]
    tooth_idx <- round(seq(1, nrow(datT), length.out = n)) # and bc those dist were scaled arbitrary
    vert_idx  <- round(seq(1, nrow(datV), length.out = n)) # we transform them back
    
    tooth_start <- datT$true_dist[tooth_idx[alignment$index1[1]]]
    tooth_end   <- datT$true_dist[tooth_idx[alignment$index1[length(alignment$index1)]]]
    vert_start  <- datV$true_dist[vert_idx[alignment$index2[1]]]
    vert_end    <- datV$true_dist[vert_idx[alignment$index2[length(alignment$index2)]]]
    
    cat(element, "\n",
        "Tooth:", tooth_start, "to", tooth_end, "\n",
        "Vert: ", vert_start, "to", vert_end, "\n",
        "Normalized distance:", alignment$normalizedDistance, "\n\n")
    
    results[[element]] <- alignment
      
      df1 <- data.frame(
        distance = datV$true_dist[vert_idx[alignment$index2]],
        value    = datT[[smooth_col]][tooth_idx[alignment$index1]],
        series   = "Tooth"
      )
      
      df2 <- data.frame(
        distance = datV$true_dist,
        value    = datV[[smooth_col]],
        series   = "Vert"
      )
      
      vert_region <- df2[df2$distance >= vert_start & df2$distance <= vert_end, ]
      
      y_scale  <- diff(range(vert_region$value, na.rm = TRUE)) / diff(range(df1$value, na.rm = TRUE))
      y_offset <- min(vert_region$value, na.rm = TRUE) - min(df1$value, na.rm = TRUE) * y_scale
      df1$value <- df1$value * y_scale + y_offset
      
      p <- ggplot() +
        geom_line(data = df2, aes(distance, value), color = "darkgrey") +
        geom_line(data = df1, aes(distance, value), color = "brown") +
        theme_bw() +
        labs(x = "Distance (um)", y = element, title = element)
      
      print(p)
    
    tryCatch({ # for when one of the silly plots fails
      plot(alignment, type = "two", off = 1, main = element)
      # contour(alignment$costMatrix,col=terrain.colors(100)) # contour is slow! at your peril
      plot(alignment, type = "threeway", off = 1, main = element)
      
    }, error = function(err) {
      cat("Skipping", element, "plots - error:", conditionMessage(err), "\n")
    })
  }
  
  return(results)
}

align_multi <- function(datT, datV, elements, n = 1000, 
                        spanT = 0.05, spanV = 0.15, ws = 200) {
  
  # smooth
  for(element in elements) {
    smooth_col <- paste0(element, "_smooth")
    datT[[smooth_col]] <- predict(loess(datT[[element]] ~ datT$true_dist, span = spanT))
    datV[[smooth_col]] <- predict(loess(datV[[element]] ~ datV$true_dist, span = spanV))
  }
  
  smooth_cols <- paste0(elements, "_smooth")
  
  # interpolate to same length
  matT <- sapply(smooth_cols, function(col) approx(seq_along(datT[[col]]), datT[[col]], n = n)$y)
  matV <- sapply(smooth_cols, function(col) approx(seq_along(datV[[col]]), datV[[col]], n = n)$y)
  
  # z-normalize each element
  matT <- apply(matT, 2, function(x) (x - mean(x)) / sd(x))
  matV <- apply(matV, 2, function(x) (x - mean(x)) / sd(x))
  
  # one warping path across all elements
  alignment <- dtw(matT, matV, keep = TRUE,
                   step = asymmetricP0,
                   open.end = FALSE,
                   open.begin = TRUE,
                   window.type = "sakoechiba",
                   window.size = ws)
  
  # back-calculate real distances
  tooth_idx <- round(seq(1, nrow(datT), length.out = n))
  vert_idx  <- round(seq(1, nrow(datV), length.out = n))
  
  tooth_start <- datT$true_dist[tooth_idx[alignment$index1[1]]]
  tooth_end   <- datT$true_dist[tooth_idx[alignment$index1[length(alignment$index1)]]]
  vert_start  <- datV$true_dist[vert_idx[alignment$index2[1]]]
  vert_end    <- datV$true_dist[vert_idx[alignment$index2[length(alignment$index2)]]]
  
  cat("Tooth:", tooth_start, "to", tooth_end, "\n",
      "Vert: ", vert_start, "to", vert_end, "\n",
      "Normalized distance:", alignment$normalizedDistance, "\n")
  
  # plot each element with single warping path
  par(mfrow = c(ceiling(length(elements)/2), 2))
  for(element in elements) {
    smooth_col <- paste0(element, "_smooth")
    plot(datV[[smooth_col]] ~ datV$true_dist, type = "l", main = element, ylab = element)
    lines(datT[[smooth_col]][tooth_idx[alignment$index1]] ~ 
            datV$true_dist[vert_idx[alignment$index2]], col = "red")
  }
  par(mfrow = c(1,1))
  
  return(alignment)
}

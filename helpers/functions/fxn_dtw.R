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

library(dtw) # fxn 1-3
library(dtwclust) # fxn 4

# compare vert & tooth, free element behaviour
align_els <- function(datT, datV, elements, n = 1000, # same steps as in thesis but a fxn now
                      pattern = asymmetric, spanT = 0.02, spanV = 0.01) {
  
  datT <- datT[order(datT$true_dist), ]
  datV <- datV[order(datV$true_dist), ]   # this sort fixes everything.
  
  # allow a single shared span (old behavior) OR a named vector per element
  if(length(spanT) == 1) spanT <- setNames(rep(spanT, length(elements)), elements)
  if(length(spanV) == 1) spanV <- setNames(rep(spanV, length(elements)), elements)
  
  for(element in elements) {
    smooth_col <- paste0(element, "_smooth")
    datT[[smooth_col]] <- predict(loess(datT[[element]] ~ datT$true_dist, span = spanT[[element]]))
    datV[[smooth_col]] <- predict(loess(datV[[element]] ~ datV$true_dist, span = spanV[[element]]))
  }
  
  results <- list()
  
  for(element in elements) {
    smooth_col <- paste0(element, "_smooth")
    
    tSmooth <- approx(seq_along(datT[[smooth_col]]), datT[[smooth_col]], n = n)$y
    vSmooth <- approx(seq_along(datV[[smooth_col]]), datV[[smooth_col]], n = n)$y
    
    tSmooth <- (tSmooth - mean(tSmooth)) / sd(tSmooth)
    vSmooth <- (vSmooth - mean(vSmooth)) / sd(vSmooth)
    
    alignment <- dtw(tSmooth, vSmooth,
                     keep = TRUE, step = pattern,
                     open.end = TRUE, open.begin = TRUE)
    
    tooth_idx <- round(seq(1, nrow(datT), length.out = n))
    vert_idx  <- round(seq(1, nrow(datV), length.out = n))
    
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
    
    tryCatch({
      plot(alignment, type = "two", offset = 1, main = element)
      plot(alignment, type = "threeway", main = element)
    }, error = function(err) {
      cat("Skipping", element, "plots - error:", conditionMessage(err), "\n")
    })
  }
  
  return(results)
}

# compare spans across elements to identify fixed location
align_grid <- function(datA, datB, 
                       elements, spans, ns, 
                       pattern = asymmetric){
  
  grid <- expand.grid(element = elements, span = spans, n = ns, stringsAsFactors = FALSE)
  grid$norm_dist  <- NA_real_
  grid$startA     <- NA_real_
  grid$endA       <- NA_real_
  grid$startB     <- NA_real_
  grid$endB       <- NA_real_
  grid$warp_slope <- NA_real_
  grid$warp_r2    <- NA_real_
  
  datA <- datA[order(datA$true_dist), ]
  datB <- datB[order(datB$true_dist), ]
  
  for(i in seq_len(nrow(grid))){
    el <- grid$element[i]
    sp <- grid$span[i]
    nn <- grid$n[i]
    
    smooth_col <- paste0(el, "_smooth")
    datA[[smooth_col]] <- predict(loess(datA[[el]] ~ datA$true_dist, span = sp))
    datB[[smooth_col]] <- predict(loess(datB[[el]] ~ datB$true_dist, span = sp))
    
    aSmooth <- approx(seq_along(datA[[smooth_col]]), datA[[smooth_col]], n = nn)$y
    bSmooth <- approx(seq_along(datB[[smooth_col]]), datB[[smooth_col]], n = nn)$y
    
    aSmooth <- (aSmooth - mean(aSmooth)) / sd(aSmooth)
    bSmooth <- (bSmooth - mean(bSmooth)) / sd(bSmooth)
    
    alignment <- tryCatch(
      dtw(aSmooth, bSmooth, keep = TRUE, step = pattern, open.end = TRUE, open.begin = TRUE),
      error = function(e) NULL
    )
    
    if(is.null(alignment)) next  # skip combos that fail (e.g. span too tight for n)
    
    idxA <- round(seq(1, nrow(datA), length.out = nn))
    idxB <- round(seq(1, nrow(datB), length.out = nn))
    
    grid$startA[i] <- datA$true_dist[idxA[alignment$index1[1]]]
    grid$endA[i]   <- datA$true_dist[idxA[alignment$index1[length(alignment$index1)]]]
    grid$startB[i] <- datB$true_dist[idxB[alignment$index2[1]]]
    grid$endB[i]   <- datB$true_dist[idxB[alignment$index2[length(alignment$index2)]]]
    grid$norm_dist[i] <- alignment$normalizedDistance
    
    fit <- lm(alignment$index2 ~ alignment$index1)
    grid$warp_slope[i] <- coef(fit)[2]
    grid$warp_r2[i]    <- summary(fit)$r.squared
  }
  
  grid
}

# helper to summarise the above output (bc I got tired of trying to visualise with ggplot)
select_spans <- function(grid, n_fixed, r2_threshold = 0.9, max_iter = 5){
  
  # filter out any groupings that result in poor correlation (r2, norm_dist)
  candidates <- grid %>% filter(n == n_fixed, warp_r2 >= r2_threshold)
  
  # flag elements where NO span works (looking at you Zn)
  all_elements <- unique(grid$element)
  qualifying   <- unique(candidates$element)
  dropped      <- setdiff(all_elements, qualifying)
  if(length(dropped) > 0){
    warning("Elements with no span reaching r2_threshold = ", r2_threshold,
            " (excluded from consensus): ", paste(dropped, collapse = ", "))
  }
  
  # per element - the span that leads you to the mean agreement 
  consensus_start <- median(candidates$startB, na.rm = TRUE)
  consensus_end   <- median(candidates$endB, na.rm = TRUE)
  
  chosen <- data.frame()
  for(iter in seq_len(max_iter)){
    chosen <- candidates %>%
      mutate(dist_to_consensus = abs(startB - consensus_start) + abs(endB - consensus_end)) %>%
      group_by(element) %>%
      slice_min(dist_to_consensus, n = 1, with_ties = FALSE) %>%
      ungroup()
    
    new_start <- median(chosen$startB, na.rm = TRUE)
    new_end   <- median(chosen$endB, na.rm = TRUE)
    
    if(abs(new_start - consensus_start) < 1 && abs(new_end - consensus_end) < 1) break
    consensus_start <- new_start
    consensus_end   <- new_end
  }
  
  list(
    chosen = chosen %>% select(element, span, warp_r2, norm_dist, startB, endB),
    consensus_start = consensus_start,
    consensus_end = consensus_end,
    spread_start = sd(chosen$startB, na.rm = TRUE),
    spread_end   = sd(chosen$endB, na.rm = TRUE),
    dropped_elements = dropped
  )
}

# (option to apply spans ID'd from fxn above)
# lock all elements in when identifying matching region
align_multi <- function(datT, datV, elements, n = 1000,
         spanT = 0.15, spanV = 0.05, # around 3x the size
         open_begin = TRUE, open_end = TRUE,
         lockV_start = NULL, pattern = asymmetric){   
  
  datT <- datT[order(datT$true_dist), ]
  datV <- datV[order(datV$true_dist), ]
  
  if(length(spanT) == 1) spanT <- setNames(rep(spanT, length(elements)), elements)
  if(length(spanV) == 1) spanV <- setNames(rep(spanV, length(elements)), elements)
  
  # smooth on FULL datV first, always — before any lock-cropping
  for(element in elements) {
    smooth_col <- paste0(element, "_smooth")
    datT[[smooth_col]] <- predict(loess(datT[[element]] ~ datT$true_dist, span = spanT[[element]]))
    datV[[smooth_col]] <- predict(loess(datV[[element]] ~ datV$true_dist, span = spanV[[element]]))
  }
  
  # NEW: crop datV to start at lockV_start, and force open.begin = FALSE on that side
  if(!is.null(lockV_start)){
    datV <- datV[datV$true_dist >= lockV_start, ]
    open_begin <- FALSE
  }
  
  smooth_cols <- paste0(elements, "_smooth")
  
  matT <- sapply(smooth_cols, function(col) approx(seq_along(datT[[col]]), datT[[col]], n = n)$y)
  matV <- sapply(smooth_cols, function(col) approx(seq_along(datV[[col]]), datV[[col]], n = n)$y)
  
  matT <- apply(matT, 2, function(x) (x - mean(x)) / sd(x))
  matV <- apply(matV, 2, function(x) (x - mean(x)) / sd(x))
  
  alignment <- dtw(matT, matV, keep = TRUE,
                   step = pattern,
                   open.end = open_end,
                   open.begin = open_begin)
  
  tooth_idx <- round(seq(1, nrow(datT), length.out = n))
  vert_idx  <- round(seq(1, nrow(datV), length.out = n))
  
  tooth_start <- datT$true_dist[tooth_idx[alignment$index1[1]]]
  tooth_end   <- datT$true_dist[tooth_idx[alignment$index1[length(alignment$index1)]]]
  vert_start  <- datV$true_dist[vert_idx[alignment$index2[1]]]
  vert_end    <- datV$true_dist[vert_idx[alignment$index2[length(alignment$index2)]]]
  
  cat("Tooth:", tooth_start, "to", tooth_end, "\n",
      "Vert: ", vert_start, "to", vert_end, "\n",
      "Normalized distance:", alignment$normalizedDistance/sqrt(length(elements)), "\n")
  
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

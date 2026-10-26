# Details ----
#' fxn_mavs.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 6 Sep 2024
#' Content: moving averages fnx for profiles
#' + separated from original script 23 Feb 2026
#' -----------

library(zoo)
add_mavs <- function(df, ratio_columns, window_size = 21,
                     gap_um = 50) { # holes wider than this (um): the median doesn't reach across
  # 10 points before and after
  if (window_size %% 2 == 0) {
    stop("Window size must be odd.")
  } #remind that can't do that silly

  # rolling median within each piece between holes (gaps in true_dist), so it doesn't mix
  # in values from across a hole; at a hole edge it uses the points on that side only.
  # the record's own start and end stay NA, as before
  o    <- if(is.null(df$true_dist)) seq_len(nrow(df)) else order(df$true_dist)
  pc   <- if(is.null(df$true_dist)) rep(1, nrow(df)) else cumsum(c(1, diff(df$true_dist[o]) > gap_um))
  half <- (window_size - 1) / 2
  n    <- nrow(df)

  for (ratio in ratio_columns) {
    y <- df[[ratio]][o]
    m <- rep(NA_real_, n)
    for(p in unique(pc)){
      i <- which(pc == p)
      m[i] <- rollapply(y[i], width = window_size, FUN = median, partial = TRUE, align = "center")
    }
    m[seq_len(min(half, n))] <- NA                  # record start
    m[seq.int(max(1, n - half + 1), n)] <- NA       # record end
    out <- rep(NA_real_, n); out[o] <- m
    df[[paste0(ratio, "_mav")]] <- out
  } #this is with rolling median, removes spikes

  return(df)
}
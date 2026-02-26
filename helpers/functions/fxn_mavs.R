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
add_mavs <- function(df, ratio_columns, window_size = 21) {
  # 10 points before and after
  if (window_size %% 2 == 0) {
    stop("Window size must be odd.")
  } #remind that can't do that silly
  
  for (ratio in ratio_columns) {
    df[[paste0(ratio, "_mav")]] <- rollmedian(df[[ratio]], 
                                              k = window_size, fill = NA, align = "center")
  } #this is with rolling median, removes spikes
  
  return(df)
}
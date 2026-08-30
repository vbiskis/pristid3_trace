# Details ----
#' fxn_profile.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 8 Jul 2025
#' Content: plots profile data for chosen elements
#'  + no changes Feb 2026
#' -----------

library(magick)

# simple and useful----
scale_In <- function(in_vals, primary_vals){
  range_p <- range(primary_vals, na.rm = TRUE)
  range_i <- range(in_vals, na.rm = TRUE)
  scale <- diff(range_p) / diff(range_i)
  shift <- range_p[1] - range_i[1] * scale
  list(scale = scale, shift = shift)
} #to get on the plot with a secondary scale

get_spikes <- function(id, threshold = 0.5){
  with(id, {
    spike <- In_Ca > threshold
    r <- rle(spike)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1
    
    keep <- r$values & r$lengths > 2
    s <- starts[keep]
    e <- ends[keep]
    
    data.frame(
      start_dist = true_dist[s],
      end_dist   = true_dist[e],
      n_points   = r$lengths[keep],
      max_In_Ca  = mapply(function(a, b) max(In_Ca[a:b]), s, e)
    )
  })
}

alt_ratio <- function(id, el, denom, use_mav = FALSE){
  in_suf  <- if(use_mav) "_Ca_mav" else "_Ca"
  out_suf <- if(use_mav) "_mav" else ""
  
  el_col    <- paste0(el, in_suf)
  denom_col <- paste0(denom, in_suf)
  
  if(!el_col %in% names(id))    stop("Column '", el_col, "' not found in id")
  if(!denom_col %in% names(id)) stop("Column '", denom_col, "' not found in id")
  
  id[[paste0(el, "_", denom, out_suf)]] <- id[[el_col]] / id[[denom_col]]
  id
}

add_alt <- function(id, elements, denoms){
  for(denom in denoms){
    els_to_run <- setdiff(elements, denom)
    for(el in els_to_run){
      id <- alt_ratio(id, el, denom = denom)
    }
  }
  id
}

# full plots----
make_plot <- function(id, y_var, color, y_lab, smooth_var = NULL, span = 0.25, 
                      show_In = TRUE, In_style = c("points", "rect"), 
                      thin = 1, rect_threshold = 0.5){
  
  In_style <- match.arg(In_style)
  smooth_var <- if(is.null(smooth_var)) y_var else smooth_var
  s <- scale_In(id$In_Ca_mav, id[[y_var]])
  id_thin <- id[order(id$true_dist), ][seq(1, nrow(id), by = thin), ]
  
  p <- ggplot(id) + theme_bw() + labs(x = "Distance (um)") +
    geom_point(data = id_thin, aes(x = true_dist, y = .data[[y_var]]), color = color, alpha = 0.4) +
    geom_smooth(aes(x = true_dist, y = .data[[smooth_var]]), method = "loess", 
                span = span, color = "black", se = TRUE, n = 300)
  
  if(show_In && In_style == "points"){
    p <- p +
      geom_point(aes(x = true_dist, y = In_Ca_mav * s$scale + s$shift), color = "grey", alpha = 0.3) +
      scale_y_continuous(name = y_lab, sec.axis = sec_axis(~ (. - s$shift) / s$scale, name = "In:Ca")) +
      theme(axis.title.y.right = element_text(color = "grey"),
            axis.text.y.right = element_text(color = "grey"))
  } else {
    p <- p + scale_y_continuous(name = y_lab)
  }
  
  p
} # now we can toggle

profs <- function(id, elrats, hardtype,
                  span = 0.25, window_size = 21, thin = 1,
                  img_path = NULL, show_In = TRUE,
                  In_style = c("points", "rect"),
                  rect_threshold = 0.5,
                  alt_denom = NULL){
  
  In_style <- match.arg(In_style)
  spikes <- get_spikes(id, threshold = rect_threshold)
  
  if(hardtype == "tooth"){
    id <- id %>% filter(In_Ca <= 1)
  } else {
    id <- id %>% filter(In_Ca <= 10)
  }
  
  all_els <- list(
    Sr = "darkgreen", Ba = "darkblue", Mn = "lightblue",
    Li = "orange",    Mg = "purple",   Zn = "brown"
  )
  
  # --- compute alt ratios FIRST, from raw _Ca columns, before any smoothing ---
  ratio_els <- NULL
  alt_cols  <- character(0)
  if(!is.null(alt_denom) && alt_denom != "Ca"){
    ratio_els <- all_els
    ratio_els[[alt_denom]] <- NULL  # don't ratio the denominator against itself
    
    for(el in names(ratio_els)){
      id <- alt_ratio(id, el, alt_denom, use_mav = FALSE) # for computing - get ratio first
    }
    alt_cols <- paste0(names(ratio_els), "_", alt_denom)
  }
  
  # smooth the original elrats AND the new raw alt-ratio columns together
  elrats_all <- unique(c(elrats, alt_cols))
  id <- id %>% group_modify(~ add_mavs(.x, elrats_all, window_size = window_size))
  
  if(hardtype == "tooth") id <- id %>% filter(In_Ca_mav <= 5)
  
  a  <<- make_plot(id, "Sr_Ca_mav", "darkgreen", "Sr:Ca", "Sr_Ca", span, show_In, In_style, thin, rect_threshold)
  b  <<- make_plot(id, "Ba_Ca_mav", "darkblue",  "Ba:Ca", "Ba_Ca", span, show_In, In_style, thin, rect_threshold)
  e  <<- make_plot(id, "Mn_Ca_mav", "lightblue", "Mn:Ca", "Mn_Ca", span, show_In, In_style, thin, rect_threshold)
  cc <<- make_plot(id, "Li_Ca_mav", "orange",    "Li:Ca", "Li_Ca", span, show_In, In_style, thin, rect_threshold)
  d  <<- make_plot(id, "Mg_Ca_mav", "purple",    "Mg:Ca", "Mg_Ca", span, show_In, In_style, thin, rect_threshold)
  f  <<- make_plot(id, "Zn_Ca_mav", "brown",      "Zn:Ca", "Zn_Ca", span, show_In, In_style, thin, rect_threshold)
  g  <<- make_plot(id, "Ca43_CPS", "pink", "Ca", "Ca43_CPS", span, show_In, In_style, thin, rect_threshold)
  
  # --- alt-ratio panels, now using the smoothed mav column for points/rects
  #     and the raw (pre-smoothing) column for the loess smooth line ---
  alt_plots <- list()
  if(length(alt_cols) > 0){
    for(el in names(ratio_els)){
      alt_col_mav <- paste0(el, "_", alt_denom, "_mav")
      alt_col_raw <- paste0(el, "_", alt_denom)
      alt_lab <- paste0(el, ":", alt_denom)
      alt_plots[[el]] <- make_plot(id, alt_col_mav, all_els[[el]], alt_lab, alt_col_raw, span,
                                   show_In, In_style, thin, rect_threshold)
    }
  }
  
  if(show_In && In_style == "rect" && nrow(spikes) > 0){
    rect_layer <- list(
      geom_rect(data = spikes, aes(xmin = start_dist, xmax = end_dist, ymin = -Inf, ymax = Inf,
                                   alpha = max_In_Ca), fill = "grey30", inherit.aes = FALSE),
      scale_alpha_binned(range = c(0.1, 0.7), breaks = c(1, 2, 4), labels = c("<1", "2", ">4"),
                         limits = c(0, 10), oob = scales::squish, name = "In:Ca peak")
    )
    a <<- a + rect_layer; b <<- b + rect_layer; e <<- e + rect_layer
    cc <<- cc + rect_layer; d <<- d + rect_layer; f <<- f + rect_layer; g <<- g + rect_layer
    if(length(alt_plots) > 0){
      for(el in names(alt_plots)) alt_plots[[el]] <- alt_plots[[el]] + rect_layer
    }
  }
  
  if(!is.null(img_path)){
    img <- image_read(img_path)
    img_plot <- ggdraw() + draw_image(img)
    
    plot(img_plot / plot_spacer() / a / b + plot_layout(heights = c(1, 0.05, 1, 1)))
    plot(img_plot / plot_spacer() / cc / d + plot_layout(heights = c(1, 0.05, 1, 1)))
    plot(img_plot / plot_spacer() / e / f + plot_layout(heights = c(1, 0.05, 1, 1)))
  } else {
    plot(a / b)
    plot(cc / d)
    plot(e / f)
  }
  
  if(length(alt_plots) > 0){
    default_map <- list(Sr = a, Ba = b, Mn = e, Li = cc, Mg = d, Zn = f)
    for(el in names(alt_plots)){
      plot(default_map[[el]] / alt_plots[[el]])
    }
  }
  
  return(id)
}

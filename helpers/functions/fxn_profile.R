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

scale_In <- function(in_vals, primary_vals){
  range_p <- range(primary_vals, na.rm = TRUE)
  range_i <- range(in_vals, na.rm = TRUE)
  scale <- diff(range_p) / diff(range_i)
  shift <- range_p[1] - range_i[1] * scale
  list(scale = scale, shift = shift)
} #to get on the plot with a secondary scale

make_plot <- function(id, y_var, color, y_lab, smooth_var = NULL, span = 0.25, show_In = TRUE, thin = 1){
  
  smooth_var <- if(is.null(smooth_var)) y_var else smooth_var
  s <- scale_In(id$In_Ca_mav, id[[y_var]])
  id_thin <- id[order(id$true_dist), ][seq(1, nrow(id), by = thin), ]
  
  p <- ggplot(id) + # repetitive ggplotting
    geom_point(data = id_thin, aes(x = true_dist, y = .data[[y_var]]), color = color, alpha = 0.4) +
    geom_smooth(aes(x = true_dist, y = .data[[smooth_var]]), method = "loess", 
                span = span, color = "black", se = TRUE) +
    theme_bw() +
    labs(x = "Distance (um)")
  
  if(show_In){
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
                  img_path = NULL, show_In = TRUE){
  
  if(hardtype == "tooth"){
    id <- id %>% filter(In_Ca <= 2)
  } else {
    id <- id %>% filter(In_Ca <= 10)
  }
  
  id <- id %>% group_modify(~ add_mavs(.x, elrats, window_size = window_size))
  
  if(hardtype == "tooth") id <- id %>% filter(In_Ca_mav <= 5)
  
  sp <- if(hardtype == "tooth") span else 0.1
  
  # for tooth, smooth var = mav; for else, smooth var = raw
  sv <- if(hardtype == "tooth") list(sr = "Sr_Ca_mav", ba = "Ba_Ca_mav", mn = "Mn_Ca_mav") else
    list(sr = "Sr_Ca",     ba = "Ba_Ca",     mn = "Mn_Ca")
  
  a  <<- make_plot(id, "Sr_Ca_mav", "darkgreen", "Sr:Ca", sv$sr, sp, show_In, thin)
  b  <<- make_plot(id, "Ba_Ca_mav", "darkblue",  "Ba:Ca", sv$ba, sp, show_In, thin)
  e  <<- make_plot(id, "Mn_Ca_mav", "lightblue", "Mn:Ca", sv$mn, sp, show_In, thin)
  cc <<- make_plot(id, "Li_Ca_mav", "orange",    "Li:Ca", NULL, sp, show_In, thin)
  d  <<- make_plot(id, "Mg_Ca_mav", "purple",    "Mg:Ca", NULL, sp, show_In, thin)
  f  <<- make_plot(id, "Zn_Ca_mav", "brown",      "Zn:Ca", NULL, sp, show_In, thin)
  g  <<- make_plot(id, "Ca43_CPS", "pink", "Ca", sv$sr, sp, show_In, thin)
  
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
    #plot(g)
  }
  
  return(id)
}

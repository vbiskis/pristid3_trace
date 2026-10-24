# Details ----
#' dtw.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 14 Aug 2026
#' Content: 
#' dynamic time warping IV
#' instead of building up sections, trialling cutting down the tooth 
#' + until the same elements are aligned in a section 
#' + and then growing it
#' -----------

## try again----

stack_alignments <- function(tooth, vert, elements,
                             spanT = 0.04, spanV = 0.01, n = 500,
                             n_steps = 5,            # cut positions tried at each end
                             max_trim = 0.35,        # trim up to this share of the tooth from each end
                             t_starts = NULL, t_ends = NULL,   # or give cut positions directly
                             tol_um = 300,           # elements agree if within this (vert um)
                             ref_step = 2500,        # tooth reference points every ref_step um
                             bands = NULL,           # optional vert band positions (um)
                             plot = c("segments", "heat", "both"),
                             ...){                   # passed to the heat map: tooth_bin, vert_bin, by_element, flip
  
  tooth <- tooth[order(tooth$true_dist), ]
  vert  <- vert[order(vert$true_dist), ]
  if(is.null(tooth$.hole)) tooth$.hole <- FALSE
  
  # loess on kept points -> grid -> NA at holes -> z-score
  sm_grid <- function(dat, el, span, grid, hole){
    d  <- dat[!dat$.hole & is.finite(dat[[el]]), ]
    y  <- approx(d$true_dist, predict(loess(d[[el]] ~ d$true_dist, span = span)),
                 xout = grid, ties = mean, rule = 2)$y
    y[hole] <- NA
    (y - mean(y, na.rm = TRUE)) / sd(y, na.rm = TRUE)
  }
  
  # vert: smoothed once
  vert$.hole <- FALSE
  vg <- seq(min(vert$true_dist), max(vert$true_dist), length.out = n)
  Vz <- sapply(elements, function(el) sm_grid(vert, el, spanV, vg, rep(FALSE, n)))
  
  # tooth cuts
  tr <- range(tooth$true_dist); L <- diff(tr)
  if(is.null(t_starts)) t_starts <- seq(tr[1], tr[1] + max_trim * L, length.out = n_steps)
  if(is.null(t_ends))   t_ends   <- seq(tr[2] - max_trim * L, tr[2], length.out = n_steps)
  cuts <- expand.grid(t_start = t_starts, t_end = t_ends)
  cuts <- cuts[cuts$t_end > cuts$t_start, ]
  cuts <- cuts[order(cuts$t_start, cuts$t_end), ]
  cuts$cut <- seq_len(nrow(cuts))
  
  runs <- list(); refs <- list(); paths <- list()
  for(k in seq_len(nrow(cuts))){
    message("Cut ", k, " of ", nrow(cuts))
    tc <- tooth[tooth$true_dist >= cuts$t_start[k] & tooth$true_dist <= cuts$t_end[k], ]
    if(nrow(tc) < 30) next
    tg <- seq(min(tc$true_dist), max(tc$true_dist), length.out = n)
    hg <- approx(tc$true_dist, as.numeric(tc$.hole), xout = tg, method = "constant",
                 rule = 2, ties = max)$y > 0
    rp <- seq(ceiling(min(tg) / ref_step) * ref_step, max(tg), by = ref_step)
    
    for(el in elements){
      res <- tryCatch({
        Tz <- sm_grid(tc, el, spanT, tg, hg)
        C  <- abs(outer(Tz, Vz[, el], "-")); C[is.na(Tz), ] <- 0
        a  <- dtw(C, step.pattern = asymmetric, open.begin = TRUE, open.end = TRUE)
        list(a = a, C = C, gap = is.na(Tz))
      }, error = function(e) NULL)
      if(is.null(res)) next
      a <- res$a
      v <- tapply(vg[a$index2], a$index1, mean)          # vert position per tooth point
      ti <- tg[as.integer(names(v))]
      paths[[length(paths) + 1]] <- data.frame(cut = k, element = el, tooth = ti, vert = as.numeric(v))
      kept <- !res$gap[a$index1]
      runs[[length(runs) + 1]] <- data.frame(
        cut = k, t_start = cuts$t_start[k], t_end = cuts$t_end[k], element = el,
        v_start = min(v), v_end = max(v),
        fit = mean(res$C[cbind(a$index1, a$index2)][kept]))   # mean |dz| over kept tooth points
      if(length(rp) > 0)
        refs[[length(refs) + 1]] <- data.frame(cut = k, element = el, tooth_ref = rp,
                                               vert = approx(ti, as.numeric(v), xout = rp)$y)
    }
  }
  runs <- do.call(rbind, runs)
  refs <- do.call(rbind, refs)
  paths <- do.call(rbind, paths)
  
  # --- groups within each cut ---
  runs$group <- NA_integer_
  groups <- do.call(rbind, lapply(split(runs, runs$cut), function(d){
    if(nrow(d) == 1) g <- 1L else {
      D <- outer(seq_len(nrow(d)), seq_len(nrow(d)), Vectorize(function(i, j)
        max(abs(d$v_start[i] - d$v_start[j]), abs(d$v_end[i] - d$v_end[j]))))
      g <- cutree(hclust(as.dist(D), "complete"), h = tol_um)
    }
    runs$group[runs$cut == d$cut[1]] <<- g
    do.call(rbind, lapply(split(seq_len(nrow(d)), g), function(i) data.frame(
      cut = d$cut[1], t_start = d$t_start[1], t_end = d$t_end[1],
      n_elements = length(i),
      members = paste(sort(sub("_Ca.*$", "", d$element[i])), collapse = "+"),
      v_start = mean(d$v_start[i]), v_end = mean(d$v_end[i]),
      spread_start = diff(range(d$v_start[i])), spread_end = diff(range(d$v_end[i])))))
  }))
  rownames(groups) <- NULL
  
  # --- uncertainty at tooth reference points ---
  ref_tab <- NULL
  if(!is.null(refs) && nrow(refs) > 0){
    ref_tab <- do.call(rbind, lapply(split(refs, refs$tooth_ref), function(d){
      d <- d[!is.na(d$vert), ]
      if(nrow(d) == 0) return(NULL)
      med <- tapply(d$vert, d$element, median)
      if(length(med) > 1){
        D <- abs(outer(as.numeric(med), as.numeric(med), "-")); g <- cutree(hclust(as.dist(D), "complete"), h = tol_um)
        best <- names(which.max(table(g)))
        mem <- names(med)[g == as.integer(best)]
        if(length(mem) < 2) return(data.frame(tooth_ref = d$tooth_ref[1], n_agree = 0,
                                              n_elements = length(med), members = "none", n_cuts = NA,
                                              vert_median = NA, vert_lo = NA, vert_hi = NA, uncert_um = NA))
      } else mem <- names(med)
      v <- d$vert[d$element %in% mem]
      data.frame(tooth_ref = d$tooth_ref[1], n_agree = length(mem), n_elements = length(med),
                 members = paste(sort(sub("_Ca.*$", "", mem)), collapse = "+"),
                 n_cuts = length(unique(d$cut[d$element %in% mem])),
                 vert_median = median(v), vert_lo = min(v), vert_hi = max(v),
                 uncert_um = max(v) - min(v))
    }))
    rownames(ref_tab) <- NULL
  }
  
  # --- plot: elements stacked within row ---
  ne  <- length(elements)
  off <- setNames(seq(-0.3, 0.3, length.out = ne), elements)
  runs$y <- runs$cut + off[runs$element]
  boxes <- groups[groups$n_elements >= 2, ]
  labs_y <- sprintf("%.0f–%.0f", cuts$t_start, cuts$t_end)
  
  p_seg <- ggplot() +
    { if(!is.null(bands)) geom_vline(xintercept = bands, colour = "grey75", linewidth = 0.3) } +
    { if(nrow(boxes) > 0) geom_rect(data = boxes,
                                    aes(xmin = v_start, xmax = v_end, ymin = cut - 0.42, ymax = cut + 0.42,
                                        alpha = n_elements), fill = "grey40") } +
    geom_segment(data = runs, aes(x = v_start, xend = v_end, y = y, yend = y, colour = element),
                 linewidth = 1.2) +
    scale_alpha_continuous(range = c(0.1, 0.35), breaks = seq(2, ne), name = "Elements\nagreeing") +
    scale_y_continuous(breaks = cuts$cut, labels = labs_y, expand = c(0.01, 0.01)) +
    theme(panel.grid.minor = element_blank(), plot.caption = element_text(hjust = 0)) +
    labs(x = "Vertebra (µm)", y = "Tooth kept (µm)", colour = NULL,
         caption = sprintf(paste0("Grey = elements whose matched vert start and end are each within %d µm ",
                                  "(complete linkage). %d tooth cuts × %d elements."),
                           tol_um, nrow(cuts), ne))
  
  # --- plot option: heat map----
  p_heat <- function(s, tooth_bin = 250, vert_bin = 50, bands = NULL, by_element = FALSE, flip = FALSE){
    d <- s$paths
    d$tb <- (floor(d$tooth / tooth_bin) + 0.5) * tooth_bin
    d$vb <- (floor(d$vert  / vert_bin)  + 0.5) * vert_bin
    hit   <- unique(d[, c("element", "cut", "tb", "vb")])
    cover <- aggregate(cut ~ element + tb, data = unique(d[, c("element", "cut", "tb")]), FUN = length)
    names(cover)[3] <- "n_cover"                                  # cuts containing this tooth bin
    agg <- aggregate(cut ~ element + tb + vb, data = hit, FUN = length)
    agg <- merge(agg, cover)
    agg$share <- agg$cut / agg$n_cover                            # share of those cuts mapping here
    ne <- length(unique(d$element))
    
    if(by_element){
      df <- agg; lim <- c(0, 1); lab <- "Share of cuts"
    } else {
      df <- aggregate(share ~ tb + vb, data = agg, FUN = sum); lim <- c(0, ne); lab <- "Elements"
    }
    ggplot(df, aes(vb, tb, fill = share)) +
      geom_tile(height = tooth_bin, width = vert_bin) +
      { if(!is.null(bands)) geom_vline(xintercept = bands, colour = "brown", linewidth = 0.5) } +
      { if(by_element) facet_wrap(~ element) } +
      scale_fill_viridis_c(limits = lim, name = lab) +
      { if(flip) coord_flip() } +
      theme(plot.caption = element_text(hjust = 0)) +
      labs(x = "Vertebra (\u00b5m)", y = "Tooth (\u00b5m)",
           caption = sprintf(paste0("S(t, v) = \u03a3\u1d62 c\u1d62(t, v) / C\u1d62(t);  c\u1d62 = cuts in which element i maps tooth bin t ",
                                    "to vert bin v, C\u1d62 = cuts containing t. Bins %d \u00d7 %d \u00b5m; %d tooth cuts."),
                             tooth_bin, vert_bin, nrow(s$cuts)))
  }
  
  out <- list(runs = runs, groups = groups, ref_points = ref_tab, ref_raw = refs,
              paths = paths, cuts = cuts, plot = p_seg)
  out$heat <- function(...) p_heat(out, bands = bands, ...)
  plot <- match.arg(plot)
  if(plot %in% c("segments", "both")) print(p_seg)
  if(plot %in% c("heat", "both"))     print(p_heat(out, bands = bands, ...))
  if(!is.null(ref_tab)) print(ref_tab, digits = 4, row.names = FALSE)
  invisible(out)
}

s <- stack_alignments(SMA001T_filt, SMA001V_filt, # sawfish 1
                      elements = mavs2, # or mavs1 for just the main 4
                      spanT = 0.05, spanV = 0.01, # adjust res again, path ~ same
                      n = 1000, n_steps = 25, # how high res the map is
                      tol_um = 800, ref_step = 1000, # how (detailed) the table prints
                      plot = "heat", # for continuous scale
                      bands = c(710, # a visual of how path relates to banding
                                1715, 2900, 5900, 
                                8500, 9175, 9600))

s$groups       # stats per tooth cut: which elements & where
s$plot
s$heat(by_element = TRUE)

t <- stack_alignments(SMA001TB_filt, SMA001VB_filt, # sawfish 1 - bleach
                      elements = mavs2, # or mavs1 for just the main 4
                      spanT = 0.05, spanV = 0.01, # adjust res again, path ~ same
                      n = 1000, n_steps = 5, # how high res the map is
                      tol_um = 800, ref_step = 1000, # how (detailed) the table prints
                      plot = "heat", # for continuous scale
                      bands = c(630, # birth
                                1210, 2400, 5600, 
                                8000, 8800, 9300))

t$plot
t$heat(by_element = TRUE)

u <- stack_alignments(SMA002T_filt, SMA002V_filt, # sawfish 1 - bleach
                      elements = mavs2, # or mavs1 for just the main 4
                      spanT = 0.05, spanV = 0.01, # adjust res again, path ~ same
                      n = 1000, n_steps = 5, # how high res the map is
                      tol_um = 800, ref_step = 1000, # how (detailed) the table prints
                      plot = "heat", # for continuous scale
                      bands = c(770, 1025, # birth and/or 1yr
                                1650, 2700, 3700,
                                4550, 7550))

u$plot
u$heat(by_element = TRUE) 
#oh Ba is useful for this one - makes sense - it lived longer in the wild.

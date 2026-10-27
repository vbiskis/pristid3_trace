# Details ----
#' match_segment()
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' Content:
#' Sequentially calc. which part of the tooth best matches a given section of vert.
#' Idea comes from the dtw package setup where:
#' + vert section = query (matched in full or pieces as called
#' + whole tooth = reference (open begin/end)
#' Under the assumption that
#' + chosen elements accumulate at the same time together (thus one shared path)
#' + cost profile along the tooth: one clear dip = well-identified match
#' **each element also matched on its own: do they pick the same spot?
#' Notes:
#' + resamples by DISTANCE (stepT / stepV um), not by point count, so the slope
#'   constraint has a physical meaning: matched tooth length can be 0 to
#'   2 * stepT / stepV times the vert section length (asymmetric pattern)
#' ** z-scores use the WHOLE record (whole tooth, whole vert), so the section keeps
#'   its level relative to the rest of the record
#' -----------

library(dtw)
library(ggplot2)
library(patchwork)
source('helpers/functions/fxn_gaps.R') # smoothing that stops at holes

# shared prep: sort, crop, smooth vert, and a tooth smoother that matches the vert
prep_records <- function(datT, datV, elements, spanV, tooth_max, min_pts = 10, gap_um = 50){
  datT <- datT[order(datT$true_dist), ]
  datT <- datT[datT$true_dist <= tooth_max, ]
  datV <- datV[order(datV$true_dist), ]
  if(length(spanV) == 1) spanV <- setNames(rep(spanV, length(elements)), elements)
  
  xT <- datT$true_dist; xV <- datV$true_dist
  LT <- diff(range(xT)); LV <- diff(range(xV))
  span_min <- min_pts / nrow(datT)    # loess needs enough points in each window
  
  zV <- scale(sapply(elements, function(el) smooth_gaps(xV, datV[[el]], spanV[[el]], Inf))) # verts: never split
  colnames(zV) <- elements
  
  cache <- list()
  # tooth smoothed for a given ratio (cached so each ratio is only fitted once)
  tooth_at_ratio <- function(r){
    key <- format(r)
    if(!is.null(cache[[key]])) return(cache[[key]])
    sp <- setNames(pmin(1, pmax(span_min, spanV * r * LV / LT)), elements)  # pmin/pmax drop name
    z  <- scale(sapply(elements, function(el) smooth_gaps(xT, datT[[el]], sp[[el]], gap_um)))
    colnames(z) <- elements
    out <- list(z = z, spanT = sp, clamped = any(spanV * r * LV / LT < span_min))
    cache[[key]] <<- out
    out
  }
  
  list(xT = xT, xV = xV, zV = zV, tooth_at_ratio = tooth_at_ratio, spanV = spanV,
       gaps = find_gaps(xT, gap_um))
}

at_pts <- function(z, x, xout, elements){
  m <- sapply(seq_len(ncol(z)), function(k) approx(x, z[, k], xout = xout, ties = mean)$y)
  matrix(m, ncol = ncol(z), dimnames = list(NULL, elements))
}

match_window <- function(datT, datV, elements, v_from, win = 1000,
                         spanV, ratio_range = c(1, 15),
                         n_ratios     = 8,
                         tooth_search = c(-Inf, Inf), # only search this part of the tooth
                         start_search = c(-Inf, Inf), # and only pieces starting in this range
                         start_step   = 100,          # um between candidate tooth starts
                         n_pts        = 100,          # points each piece is resampled to
                         window_frac  = 0.3,          # how freely each element may warp
                         keep_best    = 0.1,          # share of best-fitting candidates kept
                         tooth_max    = Inf,
                         gap_um       = 50,           # holes wider than this (um) in the tooth are skipped
                         plot         = FALSE,
                         zoom         = Inf,          # plot: windows of life shown either side (Inf = whole vert)
                         offset       = 3,
                         prep         = NULL,         # optional: reuse prep_records() output in loops
                         bands        = NULL,         # optional: vert positions (um) of the dark bands (centres)
                         band_width   = NULL,         # optional: draw bands as blocks this wide (um); NULL = lines
                         band_alpha   = 0.3,
                         show_paths   = TRUE){        # FALSE = two-way plot only
  
  if(is.null(prep)) prep <- prep_records(datT, datV, elements, spanV, tooth_max, gap_um = gap_um)
  xT <- prep$xT; xV <- prep$xV; zV <- prep$zV
  
  v_to <- v_from + win
  na_row <- function(sig){
    out <- data.frame(v_from = v_from, v_to = v_to, tooth_start = NA_real_, tooth_end = NA_real_,
                      ratio = NA_real_, spanT = NA_real_, disagree_pct = NA_real_,
                      disagree_um = NA_real_, cost = NA_real_, signal = sig)
    out[paste0("cost_", elements)] <- NA_real_
    out[paste0("warp_", elements)] <- NA_real_
    out
  }
  
  gv <- seq(v_from, v_to, length.out = n_pts)
  Vw <- at_pts(zV, xV, gv, elements)
  if(anyNA(Vw)){
    message(sprintf("Vert window %.0f-%.0f um is outside the vert data (%.0f-%.0f um) - skipped",
                    v_from, v_to, min(xV), max(xV)))
    return(na_row(NA))
  }
  signal <- mean(apply(Vw, 2, sd))   # low = flat piece, hard to place
  w      <- max(1, round(window_frac * n_pts))
  align1 <- function(q, r) dtw(q, r, step.pattern = symmetric2,
                               window.type = "sakoechiba", window.size = w)
  
  # tooth piece at n_pts even positions, minus any that fall in a hole;
  # NULL if more of it is missing than the warping window can absorb
  tooth_piece <- function(zT, s, tlen){
    pos  <- seq(s, s + tlen, length.out = n_pts)
    keep <- which(!in_gaps(pos, prep$gaps))
    if(n_pts - length(keep) > w) return(NULL)
    list(pos = pos[keep], keep = keep, z = at_pts(zT, xT, pos[keep], elements))
  }
  
  lo <- max(min(xT), tooth_search[1])
  hi <- min(max(xT), tooth_search[2])
  ratios <- seq(ratio_range[1], ratio_range[2], length.out = n_ratios)
  
  rows <- list()
  for(r in ratios){
    tlen <- r * win
    if(tlen > hi - lo) next
    zT <- prep$tooth_at_ratio(r)$z
    s_from <- max(lo, start_search[1]); s_to <- min(hi - tlen, start_search[2])
    if(s_from > s_to) next
    for(s in seq(s_from, s_to, by = start_step)){
      tp   <- tooth_piece(zT, s, tlen)
      if(is.null(tp)) next
      Tw   <- tp$z
      idx  <- matrix(NA_real_, n_pts, length(elements))
      cost <- setNames(numeric(length(elements)), paste0("cost_", elements))
      for(e in seq_along(elements)){
        a <- align1(Vw[, e, drop = FALSE], Tw[, e, drop = FALSE])
        idx[, e] <- approx(seq_along(tp$keep), tp$keep,   # tooth point for each vert point,
                           tapply(a$index2, a$index1, mean))$y  # numbered as if no hole
        cost[e]  <- a$normalizedDistance
      }
      # disagreement: at each vert point, gap between the elements whose matched
      # tooth points are furthest apart; averaged; as a share of the piece
      spread <- mean(apply(idx, 1, function(x) max(x) - min(x))) / n_pts
      # warp per element: average distance of its path from an even pace
      # (straight diagonal), as % of the piece. 0 = no warping
      warp <- setNames(100 * colMeans(abs(idx - seq_len(n_pts))) / n_pts,
                       paste0("warp_", elements))
      rows[[length(rows) + 1]] <- data.frame(ratio = r, tooth_start = s, tooth_end = s + tlen,
                                             disagree = spread, cost = mean(cost),
                                             t(cost), t(warp))
    }
  }
  if(length(rows) == 0) return(na_row(signal))
  res <- do.call(rbind, rows)
  
  # keep the best-fitting candidates at each ratio, then the one where elements warp most alike
  thr  <- ave(res$cost, res$ratio, FUN = function(x) quantile(x, keep_best))
  cand <- which(res$cost <= thr)
  b    <- res[cand[which.min(res$disagree[cand])], ]
  tr   <- prep$tooth_at_ratio(b$ratio)
  
  out <- data.frame(v_from = v_from, v_to = v_to,
                    tooth_start = b$tooth_start, tooth_end = b$tooth_end, ratio = b$ratio,
                    spanT = mean(tr$spanT),
                    disagree_pct = 100 * b$disagree,
                    disagree_um  = b$disagree * (b$tooth_end - b$tooth_start),
                    cost = b$cost, signal = signal,
                    b[grep("^cost_", names(b))], b[grep("^warp_", names(b))])
  
  # bands on a vert x-axis: lines, or blocks band_width um wide
  band_layer <- function(bb){
    if(is.null(bb) || length(bb) == 0) return(NULL)
    if(is.null(band_width))
      return(geom_vline(xintercept = bb, colour = "darkorange", linewidth = 0.4, alpha = 0.7))
    geom_rect(data = data.frame(xmin = bb - band_width / 2, xmax = bb + band_width / 2,
                                ymin = -Inf, ymax = Inf),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "darkorange", alpha = band_alpha, inherit.aes = FALSE)
  }
  
  # two-way: whole vert along the bottom, tooth piece above, grey lines = warp
  if(plot){
    tp <- tooth_piece(tr$z, b$tooth_start, b$tooth_end - b$tooth_start)
    gt <- tp$pos
    Tw <- tp$z
    tx <- v_from + (gt - b$tooth_start) / (b$tooth_end - b$tooth_start) * win
    
    vert_df <- do.call(rbind, lapply(elements, function(el)
      data.frame(element = el, x = xV, value = zV[, el])))
    tooth_df <- do.call(rbind, lapply(elements, function(el)
      data.frame(element = el, x = tx, value = Tw[, el] + offset)))
    links <- do.call(rbind, lapply(seq_along(elements), function(e){
      a <- align1(Vw[, e, drop = FALSE], Tw[, e, drop = FALSE])
      data.frame(element = elements[e],
                 x = tx[a$index2], y = Tw[a$index2, e] + offset,
                 xend = gv[a$index1], yend = Vw[a$index1, e])
    }))
    
    p <- ggplot() +
      ggplot2::annotate("rect", xmin = v_from, xmax = v_to, ymin = -Inf, ymax = Inf,
                        fill = "red", alpha = 0.06) +
      geom_segment(data = links, aes(x, y, xend = xend, yend = yend),
                   colour = "grey60", linetype = "dashed", linewidth = 0.3) +
      band_layer(bands) +
      geom_line(data = vert_df,  aes(x, value), colour = "black") +
      geom_line(data = tooth_df, aes(x, value), colour = "brown") +
      facet_wrap(~ element, ncol = 1, scales = "free_y") +
      { if(is.finite(zoom)) coord_cartesian(xlim = c(v_from - zoom * win, v_to + zoom * win)) } +
      theme_bw() +
      labs(x = "Vertebra position (um)", y = "z-score",
           caption = sprintf("Tooth %.0f-%.0f um (red) relative to vertebra %.0f-%.0f um (black)\nGrey = element warp. Ratio = %.1f (tooth span %.4f). Fit %.2f SD",
                             b$tooth_start, b$tooth_end, v_from, v_to,
                             b$ratio, mean(tr$spanT), b$cost))
    # each element's warp path vs an even pace
    paths <- do.call(rbind, lapply(seq_along(elements), function(e){
      a <- align1(Vw[, e, drop = FALSE], Tw[, e, drop = FALSE])
      data.frame(element = elements[e], vert = gv[a$index1], tooth = gt[a$index2])
    }))
    p_warp <- ggplot(paths, aes(vert, tooth, colour = element)) +
      ggplot2::annotate("segment", x = v_from, xend = v_to,
                        y = b$tooth_start, yend = b$tooth_end,
                        linetype = "dashed", colour = "grey40") +
      band_layer(bands[bands >= v_from & bands <= v_to]) +
      geom_line(linewidth = 0.8) +
      labs(x = "Vertebra position (um)", y = "Tooth position (um)", colour = NULL,
           caption = "Dashed = 0 warping. Distance from dashed line = warping value
           Warp path aligned = element agreement. Orange = vertebra banding")
    
    if(show_paths) print(p / p_warp + plot_layout(heights = c(length(elements), 1.5))) else print(p)
  }
  
  out
}

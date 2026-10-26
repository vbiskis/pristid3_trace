# Details ----
#' fxn_stack.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' Content:
#' stack_alignments(): cut the tooth down from either end, align each cut to the vert per element,
#' and see where the elements agree (moved here from 06d so the .R and .Rmd share one copy)
#' + holes in the tooth (gaps > gap_um, e.g. from removed In:Ca spikes) are smoothed either side,
#'   not across, and the grid points inside them cost nothing, so the path passes through freely
#' -----------

library(dtw)
library(ggplot2)
library(patchwork)
source('helpers/functions/fxn_gaps.R') # smoothing that stops at holes

stack_alignments <- function(tooth, vert, elements,
                             spanT = 0.04, spanV = 0.01, n = 500,
                             n_steps = 5,            # cut positions tried at each end
                             max_trim = 0.35,        # trim up to this share of the tooth from each end
                             t_starts = NULL, t_ends = NULL,   # or give cut positions directly
                             tol_um = 300,           # elements agree if within this (vert um)
                             ref_step = 2500,        # tooth reference points every ref_step um
                             bands = NULL,           # optional vert band positions (um)
                             plot = c("segments", "heat", "both", "none"),
                             gap_um = 50,            # tooth holes wider than this (um): smoothed either side, cost 0 inside
                             ...){                   # passed to the heat map: tooth_bin, vert_bin, by_element, flip
  
  tooth <- tooth[order(tooth$true_dist), ]
  vert  <- vert[order(vert$true_dist), ]
  if(is.null(tooth$.hole)) tooth$.hole <- FALSE
  
  # loess on kept points (each piece between holes on its own) -> grid -> NA at holes -> z-score
  sm_grid <- function(dat, el, span, grid, hole, gap = gap_um){
    d  <- dat[!dat$.hole & is.finite(dat[[el]]), ]
    y  <- approx(d$true_dist, smooth_gaps(d$true_dist, d[[el]], span, gap),
                 xout = grid, ties = mean, rule = 2)$y
    y[hole] <- NA
    (y - mean(y, na.rm = TRUE)) / sd(y, na.rm = TRUE)
  }
  
  # vert: smoothed once
  vert$.hole <- FALSE
  vg <- seq(min(vert$true_dist), max(vert$true_dist), length.out = n)
  Vz <- sapply(elements, function(el) sm_grid(vert, el, spanV, vg, rep(FALSE, n), gap = Inf)) # verts: never split
  
  # tooth cuts
  tr <- range(tooth$true_dist); L <- diff(tr)
  if(length(n_steps) == 1)  n_steps  <- c(n_steps, n_steps)
  if(length(max_trim) == 1) max_trim <- c(max_trim, max_trim)
  if(is.null(t_starts)) t_starts <- seq(tr[1], tr[1] + max_trim[1] * L, length.out = n_steps[1])
  if(is.null(t_ends))   t_ends   <- sort(seq(tr[2], tr[2] - max_trim[2] * L, length.out = n_steps[2]))
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
    hg <- hg | in_gaps(tg, find_gaps(tc$true_dist[!tc$.hole], gap_um))   # and grid points in a gap
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
  
  p_seg <- function(by_element = FALSE){
  r <- runs
  if(by_element) r$y <- r$cut          # one element per panel, so no stacking within the row
  ggplot() +
    { if(!is.null(bands)) geom_vline(xintercept = bands, colour = "grey75", linewidth = 0.3) } +
    { if(!by_element && nrow(boxes) > 0) geom_rect(data = boxes,
                                    aes(xmin = v_start, xmax = v_end, ymin = cut - 0.42, ymax = cut + 0.42,
                                        alpha = n_elements), fill = "grey40") } +
    geom_segment(data = r, aes(x = v_start, xend = v_end, y = y, yend = y, colour = element),
                 linewidth = 1.2) +
    { if(by_element) facet_wrap(~ element) } +
    { if(by_element) guides(colour = "none") } +
    scale_alpha_continuous(range = c(0.1, 0.35), breaks = seq(2, ne), name = "Elements\nagreeing") +
    scale_y_continuous(breaks = cuts$cut, labels = labs_y, expand = c(0.01, 0.01)) +
    theme(panel.grid.minor = element_blank(), plot.caption = element_text(hjust = 0)) +
    labs(x = "Vertebra (µm)", y = "Tooth kept (µm)", colour = NULL,
         caption = if(by_element) sprintf("%d tooth cuts × %d elements.", nrow(cuts), ne) else
           sprintf(paste0("Grey = elements whose matched vert start and end are each within %d µm ",
                          "(complete linkage). %d tooth cuts × %d elements."),
                   tol_um, nrow(cuts), ne))
  }
  
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
  
  out <- list(runs = runs, groups = groups, ref_points = ref_tab, ref_raw = refs, paths = paths, cuts = cuts)
  out$plot <- function(...) p_seg(...)
  out$heat <- function(...) p_heat(out, bands = bands, ...)
  plot <- match.arg(plot)
  if(plot %in% c("segments", "both")) print(p_seg())
  if(plot %in% c("heat", "both"))     print(p_heat(out, bands = bands, ...))
  if(!is.null(ref_tab) && plot != "none") print(ref_tab, digits = 4, row.names = FALSE)
  invisible(out)
}

# cut_trends(): does cutting the tooth move each element's match on the vertebra the way it should? ----
#' + start: for each end cut, run through the start cuts: cutting the tooth start later should start
#'   the match later on the vertebra
#' + end:   for each start cut, run through the end cuts: cutting the tooth end earlier should end the
#'   match earlier on the vertebra
#' within each, the trend shows if it holds over at least min_steps consecutive cuts (its longest run).
#' per element: for how many of the other end's cuts it shows (_n of _of), which of those cuts
#' (_at, tooth um), and the tooth -> vertebra range its runs cover (_cuts); a trend holds overall
#' (start / end) if it shows for at least half of the other end's cuts
cut_trends <- function(runs, min_steps = 3){
  run_of <- function(t, v){                       # longest run where v keeps rising with t
    o <- order(t); t <- t[o]; v <- v[o]
    up <- rle(diff(v) > 0)
    if(!any(up$values)) return(NULL)
    k <- which(up$values)[which.max(up$lengths[up$values])]
    if(up$lengths[k] + 1 < min_steps) return(NULL)
    i <- sum(up$lengths[seq_len(k - 1)]) + 1
    j <- i + up$lengths[k]
    c(t[i], t[j], v[i], v[j])
  }
  trend <- function(d, by, t, v){                  # within each cut of the other end
    g <- Filter(function(x) nrow(x) >= min_steps, split(d, d[[by]]))
    r <- lapply(g, function(x) run_of(x[[t]], x[[v]]))
    ok <- !vapply(r, is.null, logical(1))
    m <- do.call(rbind, r[ok])
    list(n = sum(ok), of = length(r),
         at = if(any(ok)) sprintf("%.0f-%.0f", min(as.numeric(names(r)[ok])), max(as.numeric(names(r)[ok]))) else "-",
         cuts = if(any(ok)) sprintf("tooth %.0f-%.0f -> vert %.0f-%.0f", min(m[, 1]), max(m[, 2]), min(m[, 3]), max(m[, 4])) else "-")
  }
  out <- do.call(rbind, lapply(split(runs, runs$element), function(d){
    s <- trend(d, "t_end", "t_start", "v_start")
    e <- trend(d, "t_start", "t_end", "v_end")
    data.frame(element = sub("_Ca.*$", "", d$element[1]),
               start = s$n > 0 && s$n >= s$of / 2, start_n = s$n, start_of = s$of, start_at = s$at, start_cuts = s$cuts,
               end   = e$n > 0 && e$n >= e$of / 2, end_n = e$n,   end_of = e$of,   end_at = e$at,   end_cuts = e$cuts)
  }))
  rownames(out) <- NULL
  out
}

# sweep_spans(): cut_trends() for every spanT x spanV combination ----
#' everything else passed to stack_alignments() through ... (n, n_steps, gap_um, ...)
sweep_spans <- function(tooth, vert, elements, spansT, spansV, min_steps = 3, workers = 1, ...){
  grid <- expand.grid(spanT = spansT, spanV = spansV)
  one  <- function(k) cbind(grid[k, ], cut_trends(suppressMessages(stack_alignments(
    tooth, vert, elements, spanT = grid$spanT[k], spanV = grid$spanV[k], plot = "none", ...))$runs, min_steps),
    row.names = NULL)
  if(workers > 1){
    old <- future::plan(future::multisession, workers = workers)
    on.exit(future::plan(old), add = TRUE)
    out <- furrr::future_map(seq_len(nrow(grid)), one, .options = furrr::furrr_options(packages = c("dtw", "ggplot2")))
  } else out <- lapply(seq_len(nrow(grid)), one)
  do.call(rbind, out)
}

# plot_trends(): which trends hold, per element and span combination ----
plot_trends <- function(tab, min_steps = 3){
  tab$trend <- factor(ifelse(tab$start & tab$end, "start and end", ifelse(tab$start, "start only",
                      ifelse(tab$end, "end only", "neither"))),
                      levels = c("start and end", "start only", "end only", "neither"))
  ggplot(tab, aes(factor(spanV), factor(spanT), fill = trend)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%d/%d\n%d/%d", start_n, start_of, end_n, end_of)), size = 3.2) +
    facet_wrap(~ element) +
    scale_fill_manual(values = c("start and end" = "#1b9e77", "start only" = "#a6d96a",
                                 "end only" = "#fdae61", "neither" = "grey85"), drop = FALSE) +
    theme_bw() +
    theme(plot.caption = element_text(hjust = 0)) +
    labs(x = "spanV (loess span on the vertebra)", y = "spanT (loess span on the tooth)", fill = NULL,
         caption = paste(strwrap(paste(
           "Each tile is one span combination for one element; labels = start n/of, end n/of. Start: for each",
           "end cut, cutting the tooth start later made the match start later on the vertebra over at least %d",
           "consecutive cuts - n of the end cuts. End: for each start cut, cutting the tooth end earlier made the",
           "match end earlier over at least %d consecutive cuts - n of the start cuts. A trend holds (colour) if",
           "it shows for at least half of the other end's cuts."), 115), collapse = "\n") |> sprintf(min_steps, min_steps))
}

# Details ----
#' fxn_gaps.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' Content:
#' smoothing that respects holes (gaps in true_dist left by removing In:Ca spikes)
#' + one loess over the whole record smooths ACROSS a hole, blending the data either side
#' + resampling onto an even distance grid then draws a straight line through the hole
#' these fit each continuous piece on its own, and find grid points that fall in a hole
#' + gaps: normal point spacing is ~8.6 um, so anything over gap_um = 50 um is a hole
#' -----------

library(ggplot2)

# holes: gaps between neighbouring points wider than gap_um
find_gaps <- function(x, gap_um = 50){
  x <- sort(x[is.finite(x)])
  i <- which(diff(x) > gap_um)
  data.frame(start = x[i], end = x[i + 1])
}

# TRUE for positions inside a hole
in_gaps <- function(x, gaps){
  out <- rep(FALSE, length(x))
  for(k in seq_len(nrow(gaps))) out <- out | (x > gaps$start[k] & x < gaps$end[k])
  out
}

# loess fitted separately on each piece between holes
# span is a share of the WHOLE record, so each window holds the same number of points
# as one loess over everything (identical result when there are no holes)
# at a hole edge the usual quadratic loess curls past the last points, so within half a
# window of a hole it eases into a local linear fit; a piece shorter than one window
# next to a hole is linear throughout (the record's own start and end are left as they were)
# se = TRUE returns a data frame (fit, se, piece) for plotting ribbons
smooth_gaps <- function(x, y, span, gap_um = 50, min_pts = 10, se = FALSE){
  ok <- which(is.finite(x) & is.finite(y))
  o  <- ok[order(x[ok])]; xs <- x[o]; ys <- y[o]; N <- length(xs)
  pc <- cumsum(c(1, diff(xs) > gap_um))
  h  <- span * N / 2                          # half a smoothing window, in points

  fit <- sef <- rep(NA_real_, N)
  for(p in unique(pc)){
    i <- which(pc == p); n <- length(i)
    d <- data.frame(x = xs[i], y = ys[i])
    if(n < min_pts){ fit[i] <- mean(d$y); next }
    sp <- min(1, max(span * N, min_pts) / n)   # at least min_pts points per window

    # weight on the quadratic: 1 in the middle, 0 at an edge that's a hole
    from_hole <- pmin(if(p > 1) seq_len(n) else Inf, if(p < max(pc)) n:1 else Inf)
    w <- if(n < 2 * h && max(pc) > 1) rep(0, n) else pmin(1, from_hole / h)

    q <- predict(loess(y ~ x, data = d, span = sp, degree = 2), se = se)
    if(all(w == 1)){
      if(se){ fit[i] <- q$fit; sef[i] <- q$se.fit } else fit[i] <- q
      next
    }
    l <- predict(loess(y ~ x, data = d, span = sp, degree = 1), se = se)
    if(se){
      fit[i] <- w * q$fit + (1 - w) * l$fit
      sef[i] <- w * q$se.fit + (1 - w) * l$se.fit
    } else fit[i] <- w * q + (1 - w) * l
  }
  out <- rep(NA_real_, length(x)); out[o] <- fit
  if(!se) return(out)
  s <- p_ <- rep(NA_real_, length(x)); s[o] <- sef; p_[o] <- pc
  data.frame(fit = out, se = s, piece = p_)
}

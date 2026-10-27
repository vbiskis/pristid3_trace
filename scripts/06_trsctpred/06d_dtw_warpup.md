06d · DTW warp up
================
Nikki Biskis
Rendered 27 Oct 2026 01:54

- [Settings](#settings)
- [Searches (all sets)](#searches-all-sets)
- [SMA001 (water)](#sma001-water)
- [SMA001 (bleach)](#sma001-bleach)
- [SMA002 (water)](#sma002-water)

- **Paper:** The potential of elemental analysis in determining origin
  of sawfish rostra
- **DOI:** tbd
- **Written:** 28 Aug 2026 (all three sets, all six elements, Oct 2026)
- **Content:** dynamic time warping III: build the match up along the
  vertebra, window by window

The original version is in `original/06c_dtw_warpup.R`. For each set,
with all six elements together in `match_window()`:

1.  **Mg’s window**: one window over the stretch of vertebra where Mg
    places the whole tooth (06c), to see how every element fits there.
2.  **Sequential**: 500 µm vertebra windows every 200 µm from the start
    of the vertebra, each forced to match later on the tooth than the
    last.
3.  **Random order**: the same windows, each on its own, free to match
    anywhere on the tooth.
4.  **Anchor and build-out**: the random-order windows inside the anchor
    (the vertebra range where the set’s anchor pair maps the same tooth
    to the same place, from 06c), then built out forward (later vertebra
    must match later tooth) and backward (earlier vertebra, earlier
    tooth).
5.  **Final window**: one `match_window()` over the whole built-out
    stretch, the final output.

Starts from 06a and 06c (`data/output/dtw_cache/06a.RData`,
`06c.RData`); saves `06d.RData` for 06e and 06f. Vertebra spans are each
set’s spans in 06c.

``` r
source('helpers/microHelp.R')
source('helpers/functions/fxn_mavs.R')
source('helpers/functions/fxn_dtw.R')
source('helpers/functions/fxn_warps.R')
library(furrr)
load('data/output/dtw_cache/06a.RData') # objects from 06a
load('data/output/dtw_cache/06c.RData') # sets, gap_um, down (runs, spans), anchors from 06c
```

## Settings

Mg’s window and the anchor for each set (vertebra µm), both from 06c:
Mg’s window is where Mg places the whole, uncut tooth; the anchor is
where the set’s anchor pair agrees. Either can be set by hand below.

``` r
# Mg's window: where Mg places the whole (uncut) tooth in 06c's run; anchor: where the set's
# anchor pair agrees (06c)
mg_window <- function(n){
  r <- down[[n]]$runs
  r <- r[r$element == "Mg_Ca_mav" & r$t_start == min(r$t_start) & r$t_end == max(r$t_end), ]
  c(r$v_start, r$v_end)
}
windows <- lapply(setNames(names(sets), names(sets)), function(n)
  list(mg = mg_window(n), anchor = anchors[[n]]))
# to set one by hand, e.g.: windows$SMA002W$anchor <- c(5000, 6000)
missing <- names(windows)[sapply(windows, function(w) is.null(w$anchor) || is.null(w$mg))]
if(length(missing)) stop("no Mg window or anchor from 06c for ", paste(missing, collapse = ", "),
                         ": set it by hand above")

up_settings <- list(elements = mavs2, win = 500, step = 200, ratio_range = c(3, 16), n_ratios = 8)

# one window of one set, all elements (prep reused so each record is only smoothed once)
preps <- lapply(setNames(names(sets), names(sets)), function(n)
  prep_records(sets[[n]]$tooth, sets[[n]]$vert, up_settings$elements, down[[n]]$spanV, Inf, gap_um = gap_um))
mw <- function(n, v, search, win = up_settings$win, ...)
  cbind(set = n, match_window(sets[[n]]$tooth, sets[[n]]$vert, up_settings$elements, v_from = v, win = win,
                              spanV = down[[n]]$spanV, tooth_search = search, gap_um = gap_um,
                              prep = preps[[n]], ...))
starts <- function(n) seq(min(sets[[n]]$vert$true_dist), max(sets[[n]]$vert$true_dist) - up_settings$win,
                          by = up_settings$step)

# matched tooth against vertebra, as in the original 06c
plot_res <- function(d) ggplot(d, aes((v_from + v_to) / 2, (tooth_start + tooth_end) / 2,
                                      colour = disagree_pct, size = signal)) +
  geom_point() + theme_bw() +
  theme(plot.caption = element_text(hjust = 0)) +
  labs(x = "Vert (um)", y = "Matched tooth (um)", colour = "Elements\ndisagree (%)", size = "Signal",
       caption = paste("Each point is one 500 um vertebra window: x = window centre, y = centre of the matched tooth piece.",
                       "\nDisagree = how far apart the elements place the piece; signal = how much the vertebra window varies."))
```

## Searches (all sets)

The random-order windows of all three sets run in parallel; then each
set’s sequential search and build-out (which have to go in order) run in
parallel across sets.

``` r
jobs <- do.call(rbind, lapply(names(sets), function(n) data.frame(set = n, v = starts(n))))
plan(multisession, workers = future::availableCores() - 1)
rand <- future_map2_dfr(jobs$set, jobs$v, function(n, v)
  mw(n, v, c(-Inf, Inf), ratio_range = up_settings$ratio_range, n_ratios = up_settings$n_ratios),
  .options = furrr_options(packages = c("dtw", "ggplot2")))
```

``` r
seq_build <- function(n){
  args <- list(ratio_range = up_settings$ratio_range, n_ratios = up_settings$n_ratios)
  # sequential: each window later on the tooth than the last
  sq <- list(); search <- c(-Inf, Inf)
  for(v in starts(n)){
    r <- do.call(mw, c(list(n, v, search), args)); sq[[length(sq) + 1]] <- r
    if(!is.na(r$tooth_start)) search <- c(r$tooth_start, Inf)
  }
  # anchor: random-order windows inside the anchor range; build out from its ends
  a <- rand[rand$set == n & rand$v_from >= windows[[n]]$anchor[1] & rand$v_from <= windows[[n]]$anchor[2] &
              !is.na(rand$tooth_start), ]
  a <- a[order(a$v_from), ]
  a1 <- a[1, ]; aN <- a[nrow(a), ]
  fwd <- list(); search <- c(aN$tooth_start, Inf)          # forward: later vert must match later tooth
  for(v in starts(n)[starts(n) > aN$v_from]){
    r <- do.call(mw, c(list(n, v, search), args)); fwd[[length(fwd) + 1]] <- r
    if(!is.na(r$tooth_start)) search <- c(r$tooth_start, Inf)
  }
  bwd <- list(); search <- c(-Inf, a1$tooth_end)            # backward: earlier vert must match earlier tooth
  for(v in rev(starts(n)[starts(n) < a1$v_from])){
    r <- do.call(mw, c(list(n, v, search), args)); bwd[[length(bwd) + 1]] <- r
    if(!is.na(r$tooth_end)) search <- c(-Inf, r$tooth_end)
  }
  list(seq = do.call(rbind, sq), anchor = a, built = rbind(do.call(rbind, bwd), a, do.call(rbind, fwd)))
}
plan(multisession, workers = length(sets))
up <- future_map(setNames(names(sets), names(sets)), seq_build,
                 .options = furrr_options(packages = c("dtw", "ggplot2")))
plan(sequential)
```

## SMA001 (water)

``` r
mg_check <- function(n) mw(n, windows[[n]]$mg[1], c(-Inf, Inf), win = diff(windows[[n]]$mg),
                           ratio_range = c(3, 16), n_ratios = 24, plot = TRUE, bands = sets[[n]]$bands)
invisible(mg_check("SMA001W"))
```

![](06d_dtw_warpup_files/figure-gfm/w-mg-1.png)<!-- -->

``` r
plot_res(up$SMA001W$seq)
```

![](06d_dtw_warpup_files/figure-gfm/w-seq-1.png)<!-- -->

``` r
plot_res(rand[rand$set == "SMA001W", ])
```

![](06d_dtw_warpup_files/figure-gfm/w-rand-1.png)<!-- -->

``` r
plot_res(up$SMA001W$built)
```

![](06d_dtw_warpup_files/figure-gfm/w-built-1.png)<!-- -->

``` r
final_window <- function(n){
  b <- up[[n]]$built[!is.na(up[[n]]$built$tooth_start), ]
  mw(n, min(b$v_from), c(-Inf, Inf), win = max(b$v_to) - min(b$v_from),
     ratio_range = c(3, 24), n_ratios = 24, plot = TRUE, bands = sets[[n]]$bands)
}
final <- list()
final$SMA001W <- final_window("SMA001W")
```

![](06d_dtw_warpup_files/figure-gfm/w-final-1.png)<!-- -->

## SMA001 (bleach)

``` r
invisible(mg_check("SMA001B"))
```

![](06d_dtw_warpup_files/figure-gfm/b-mg-1.png)<!-- -->

``` r
plot_res(up$SMA001B$seq)
```

![](06d_dtw_warpup_files/figure-gfm/b-seq-1.png)<!-- -->

``` r
plot_res(rand[rand$set == "SMA001B", ])
```

![](06d_dtw_warpup_files/figure-gfm/b-rand-1.png)<!-- -->

``` r
plot_res(up$SMA001B$built)
```

![](06d_dtw_warpup_files/figure-gfm/b-built-1.png)<!-- -->

``` r
final$SMA001B <- final_window("SMA001B")
```

![](06d_dtw_warpup_files/figure-gfm/b-final-1.png)<!-- -->

## SMA002 (water)

``` r
invisible(mg_check("SMA002W"))
```

![](06d_dtw_warpup_files/figure-gfm/s2-mg-1.png)<!-- -->

``` r
plot_res(up$SMA002W$seq)
```

![](06d_dtw_warpup_files/figure-gfm/s2-seq-1.png)<!-- -->

``` r
plot_res(rand[rand$set == "SMA002W", ])
```

![](06d_dtw_warpup_files/figure-gfm/s2-rand-1.png)<!-- -->

``` r
plot_res(up$SMA002W$built)
```

![](06d_dtw_warpup_files/figure-gfm/s2-built-1.png)<!-- -->

``` r
final$SMA002W <- final_window("SMA002W")
```

![](06d_dtw_warpup_files/figure-gfm/s2-final-1.png)<!-- -->

``` r
final <- do.call(rbind, final)
save(windows, up_settings, rand, up, final, file = 'data/output/dtw_cache/06d.RData')
```

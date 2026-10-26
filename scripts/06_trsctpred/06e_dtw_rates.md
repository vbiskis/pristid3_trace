06e · DTW final matches
================
Nikki Biskis
Rendered 26 Oct 2026 20:43

- **Paper:** The potential of elemental analysis in determining origin
  of sawfish rostra
- **DOI:** tbd
- **Written:** 30 Aug 2026 (final table Oct 2026)
- **Content:** dynamic time warping V: where the tooth sits on the
  vertebra, four ways, for each set

Starts from 06a, 06c and 06d (`data/output/dtw_cache/`). The earlier
version is in `original/06e_dtw_rates.R`. For each set:

1.  **Mg, warp down** (06c): where cutting the tooth moves Mg’s match
    the way it should.
2.  **Mg, warp up** (06d): Mg’s longest run of vertebra windows where
    the random-order and sequential matches agree.
3.  **Forced, warp-down mix**: `align_multi()` on the whole tooth and
    vertebra (open start and end), all elements of the warp-down mix
    locked to one path.
4.  **Forced, warp-up mix**: the same with the warp-up mix.

Spans for the forced alignments are each set’s spans in 06c.

``` r
source('helpers/microHelp.R')
source('helpers/functions/fxn_mavs.R')
source('helpers/functions/fxn_dtw.R')
source('helpers/functions/fxn_stack.R')
load('data/output/dtw_cache/06a.RData') # objects from 06a
load('data/output/dtw_cache/06c.RData') # sets, gap_um, down, down_mix, down_best from 06c
load('data/output/dtw_cache/06d.RData') # up_windows, up_best, up_mix from 06d

# align_multi() on the whole records with one mix of elements; where it starts and ends on each
force_align <- function(name, mix, n = 1000){
  set <- sets[[name]]; d <- down[[name]]
  al <- align_multi(set$tooth, set$vert, mix, n = n, spanT = d$spanT, spanV = d$spanV,
                    open_begin = TRUE, open_end = TRUE, gap_um = gap_um)
  tt <- sort(set$tooth$true_dist)[round(seq(1, nrow(set$tooth), length.out = n))]
  vv <- sort(set$vert$true_dist)[round(seq(1, nrow(set$vert), length.out = n))]
  data.frame(tooth_from = tt[al$index1[1]], tooth_to = tt[al$index1[length(al$index1)]],
             vert_from  = vv[al$index2[1]], vert_to  = vv[al$index2[length(al$index2)]],
             norm_dist  = al$normalizedDistance)
}
```

``` r
final_tab <- do.call(rbind, lapply(names(sets), function(name){
  mg_down <- down_best[down_best$set == name, ]
  mg_up   <- up_best[up_best$set == name & up_best$element == "Mg", ]
  invisible(capture.output({
    f_down <- force_align(name, down_mix[[name]])
    f_up   <- force_align(name, up_mix[[name]])
  }))
  els <- function(m) paste(sub("_Ca.*$", "", m), collapse = "+")
  rbind(
    data.frame(set = name, match = "Mg, warp down", elements = "Mg",
               mg_down[c("tooth_from", "tooth_to", "vert_from", "vert_to")], norm_dist = NA, agree_pct = NA),
    data.frame(set = name, match = "Mg, warp up", elements = "Mg",
               mg_up[c("tooth_from", "tooth_to", "vert_from", "vert_to")], norm_dist = NA, agree_pct = mg_up$agree_pct),
    data.frame(set = name, match = "Forced, warp-down mix", elements = els(down_mix[[name]]),
               f_down, agree_pct = NA),
    data.frame(set = name, match = "Forced, warp-up mix", elements = els(up_mix[[name]]),
               f_up, agree_pct = NA))
}))
```

``` r
final_tab$R <- (final_tab$tooth_to - final_tab$tooth_from) / (final_tab$vert_to - final_tab$vert_from)
final_tab[c("tooth_from", "tooth_to", "vert_from", "vert_to")] <- lapply(final_tab[c("tooth_from", "tooth_to", "vert_from", "vert_to")], round)
final_tab$R <- round(final_tab$R, 1); final_tab$norm_dist <- round(final_tab$norm_dist, 2)
rownames(final_tab) <- NULL
```

``` r
knitr::kable(final_tab)
```

| set | match | elements | tooth_from | tooth_to | vert_from | vert_to | norm_dist | agree_pct | R |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|
| SMA001W | Mg, warp down | Mg | 12860 | 43524 | 3101 | 7033 | NA | NA | 7.8 |
| SMA001W | Mg, warp up | Mg | 10521 | 12821 | 711 | 1611 | NA | 6 | 2.6 |
| SMA001W | Forced, warp-down mix | Mg+Zn | 821 | 54875 | 3597 | 10268 | 0.81 | NA | 8.1 |
| SMA001W | Forced, warp-up mix | Mg+Sr | 821 | 54875 | 711 | 1915 | 1.05 | NA | 44.9 |
| SMA001B | Mg, warp down | Mg | 44295 | 50018 | 8563 | 9827 | NA | NA | 4.5 |
| SMA001B | Mg, warp up | Mg | 35502 | 37002 | 746 | 1246 | NA | 2 | 3.0 |
| SMA001B | Forced, warp-down mix | Li+Mn | 12002 | 52880 | 5878 | 9965 | 0.72 | NA | 10.0 |
| SMA001B | Forced, warp-up mix | Mn+Zn | 12002 | 52880 | 4510 | 9513 | 0.63 | NA | 8.2 |
| SMA002W | Mg, warp down | Mg | NA | NA | NA | NA | NA | NA | NA |
| SMA002W | Mg, warp up | Mg | 20812 | 22312 | 334 | 834 | NA | 3 | 3.0 |
| SMA002W | Forced, warp-down mix | Zn+Ba | 112 | 43012 | 2177 | 5806 | 0.59 | NA | 11.8 |
| SMA002W | Forced, warp-up mix | Zn+Ba | 112 | 43012 | 2177 | 5806 | 0.59 | NA | 11.8 |

- **match**: how the tooth was placed on the vertebra (see the list at
  the top).
- **elements**: elements used.
- **tooth_from, tooth_to, vert_from, vert_to**: where the match starts
  and ends on the tooth and on the vertebra (µm).
- **R**: growth ratio over the match, µm of tooth per µm of vertebra.
- **norm_dist**: fit of the forced alignments (DTW normalised distance,
  combined over the elements; lower = closer match).
- **agree_pct**: for Mg warp up, the share of vertebra windows where the
  random-order and sequential matches agree.

``` r
dir.create('data/output/dtw_final', showWarnings = FALSE)
write.csv(final_tab, 'data/output/dtw_final/final.csv', row.names = FALSE)
write.csv(up_best, 'data/output/dtw_final/warpup_best.csv', row.names = FALSE)
save(final_tab, force_align, file = 'data/output/dtw_cache/06e.RData')
```

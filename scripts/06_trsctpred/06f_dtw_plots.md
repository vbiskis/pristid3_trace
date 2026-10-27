06e · DTW summary
================
Nikki Biskis
Rendered 27 Oct 2026 02:11

- [Summary table](#summary-table)
- [Figure](#figure)

- **Paper:** The potential of elemental analysis in determining origin
  of sawfish rostra
- **DOI:** tbd
- **Written:** summary Oct 2026
- **Content:** plots: by window build-out, final window and forced
  alignment

Starts from 06a - 06d (`data/output/dtw_cache/`). The earlier version is
in `original/06e_dtw_rates.R`. For each set, three results are compared,
all with the six elements together and each set’s spans from 06c:

- **Build-out**: 06d’s 500 µm vertebra windows, from the anchor out
  (forward and backward).
- **Final window**: 06d’s one `match_window()` over the whole built-out
  stretch.
- **align_multi**: one DTW path over the whole tooth and vertebra (open
  start and end).

``` r
source('helpers/microHelp.R')
source('helpers/functions/fxn_mavs.R')
source('helpers/functions/fxn_dtw.R')
load('data/output/dtw_cache/06a.RData') # objects from 06a
load('data/output/dtw_cache/06c.RData') # sets, gap_um, down from 06c
load('data/output/dtw_cache/06d.RData') # windows, up_settings, rand, up, final from 06d
```

``` r
# align_multi() on the whole records, all elements; its path in um
force_align <- function(n, k = 1000){
  set <- sets[[n]]; d <- down[[n]]
  invisible(capture.output(al <- align_multi(set$tooth, set$vert, up_settings$elements, n = k,
                                             spanT = d$spanT, spanV = d$spanV, open_begin = TRUE,
                                             open_end = TRUE, gap_um = gap_um)))
  tt <- sort(set$tooth$true_dist)[round(seq(1, nrow(set$tooth), length.out = k))]
  vv <- sort(set$vert$true_dist)[round(seq(1, nrow(set$vert), length.out = k))]
  list(path = data.frame(set = n, vert = vv[al$index2], tooth = tt[al$index1]),
       norm_dist = al$normalizedDistance)
}
fa <- lapply(setNames(names(sets), names(sets)), force_align)
```

## Summary table

``` r
built <- do.call(rbind, lapply(names(up), function(n) cbind(set = n, up[[n]]$built[-1])))
built <- built[!is.na(built$tooth_start), ]
built$v_mid <- (built$v_from + built$v_to) / 2
built$t_mid <- (built$tooth_start + built$tooth_end) / 2

summary_tab <- do.call(rbind, lapply(names(sets), function(n){
  b <- built[built$set == n, ]; f <- final[final$set == n, ]; p <- fa[[n]]$path
  fit <- lm(t_mid ~ v_mid, data = b)
  rbind(
    data.frame(set = n, result = "Build-out", windows = sprintf("%d / %d", nrow(b), nrow(up[[n]]$built)),
               vert_from = min(b$v_from), vert_to = max(b$v_to), tooth_from = min(b$tooth_start), tooth_to = max(b$tooth_end),
               R = coef(fit)[2], R_r2 = summary(fit)$r.squared, fit = median(b$cost), disagree_pct = median(b$disagree_pct)),
    data.frame(set = n, result = "Final window", windows = "1",
               vert_from = f$v_from, vert_to = f$v_to, tooth_from = f$tooth_start, tooth_to = f$tooth_end,
               R = f$ratio, R_r2 = NA, fit = f$cost, disagree_pct = f$disagree_pct),
    data.frame(set = n, result = "align_multi", windows = "",
               vert_from = min(p$vert), vert_to = max(p$vert), tooth_from = min(p$tooth), tooth_to = max(p$tooth),
               R = (max(p$tooth) - min(p$tooth)) / (max(p$vert) - min(p$vert)), R_r2 = NA,
               fit = fa[[n]]$norm_dist, disagree_pct = NA))
}))
num <- c("vert_from", "vert_to", "tooth_from", "tooth_to")
summary_tab[num] <- lapply(summary_tab[num], round)
summary_tab[c("R", "R_r2", "fit", "disagree_pct")] <- lapply(summary_tab[c("R", "R_r2", "fit", "disagree_pct")], round, 2)
rownames(summary_tab) <- NULL
knitr::kable(summary_tab)
```

| set | result | windows | vert_from | vert_to | tooth_from | tooth_to | R | R_r2 | fit | disagree_pct |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| SMA001W | Build-out | 48 / 48 | 711 | 10611 | 1021 | 53692 | 6.23 | 0.93 | 0.89 | 25.78 |
| SMA001W | Final window | 1 | 711 | 10611 | 3221 | 50999 | 4.83 | NA | 0.42 | 28.56 |
| SMA001W | align_multi |  | 3597 | 5072 | 821 | 54875 | 36.66 | NA | 2.15 | NA |
| SMA001B | Build-out | 44 / 44 | 746 | 9846 | 12002 | 52502 | 4.99 | 0.92 | 0.53 | 33.00 |
| SMA001B | Final window | 1 | 746 | 9846 | 15902 | 43202 | 3.00 | NA | 0.30 | 26.93 |
| SMA001B | align_multi |  | 4116 | 9476 | 12002 | 52880 | 7.63 | NA | 2.04 | NA |
| SMA002W | Build-out | 36 / 36 | 334 | 7834 | 912 | 41312 | 4.54 | 0.66 | 0.54 | 31.75 |
| SMA002W | Final window | 1 | 334 | 7834 | 15012 | 37512 | 3.00 | NA | 0.33 | 21.58 |
| SMA002W | align_multi |  | 3908 | 4185 | 112 | 43012 | 155.16 | NA | 2.30 | NA |

- **windows**: build-out windows that matched, of all windows.
- **vert_from, vert_to, tooth_from, tooth_to**: where the result starts
  and ends on the vertebra and the tooth (µm).
- **R**: growth ratio, µm of tooth per µm of vertebra. For the
  build-out, the slope of matched tooth position on vertebra position
  over its windows, and **R_r2** how straight that is (1 = constant
  growth ratio).
- **fit**: DTW fit (normalised distance, z-score units; lower = closer
  match). For the build-out, the median over its windows.
- **disagree_pct**: how far apart the elements place the matched tooth
  piece, as a share of the piece (build-out: median over its windows).

## Figure

``` r
bands <- do.call(rbind, lapply(names(sets), function(n) data.frame(set = n, band = sets[[n]]$bands)))
anch  <- do.call(rbind, lapply(names(windows), function(n) data.frame(set = n, from = windows[[n]]$anchor[1], to = windows[[n]]$anchor[2])))
paths <- do.call(rbind, lapply(fa, `[[`, "path"))
p <- ggplot() +
  geom_rect(data = anch, aes(xmin = from, xmax = to, ymin = -Inf, ymax = Inf), fill = "grey88") +
  geom_vline(data = bands, aes(xintercept = band), colour = "brown", linewidth = 0.3) +
  geom_line(data = paths, aes(vert, tooth), colour = "#1b9e77", linewidth = 0.8) +
  geom_segment(data = final, aes(x = v_from, xend = v_to, y = tooth_start, yend = tooth_end),
               colour = "#d95f02", linewidth = 1.2) +
  geom_point(data = built, aes(v_mid, t_mid, size = signal, colour = disagree_pct)) +
  scale_colour_viridis_c(name = "Elements\ndisagree (%)") +
  scale_size_continuous(range = c(0.5, 3), name = "Signal") +
  facet_wrap(~ set, ncol = 1, scales = "free") +
  theme_bw() +
  theme(plot.caption = element_text(hjust = 0)) +
  labs(x = "Vertebra (um)", y = "Tooth (um)",
       caption = paste(strwrap(paste(
         "Points = build-out windows (x = vertebra window centre, y = centre of the matched tooth piece; colour = how",
         "far apart the elements place it; size = how much the vertebra window varies). Orange line = the final",
         "match_window() over the built-out stretch (start to end; slope = growth ratio). Green line = align_multi()",
         "path over the whole records. Grey band = anchor; brown lines = vertebral bands."), 115), collapse = "\n"))
p
```

![](06e_dtw_rates_files/figure-gfm/figure-1.png)<!-- -->

``` r
dir.create('data/output/dtw_final', showWarnings = FALSE)
dir.create('figures/dtw', showWarnings = FALSE)
write.csv(summary_tab, 'data/output/dtw_final/summary.csv', row.names = FALSE)
ggsave('figures/dtw/summary.png', p, width = 10, height = 9, dpi = 300)
save(summary_tab, built, fa, file = 'data/output/dtw_cache/06e.RData')
```

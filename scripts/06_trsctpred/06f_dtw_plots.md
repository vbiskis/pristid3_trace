06f · DTW plots
================
Nikki Biskis
Rendered 27 Oct 2026 21:19

- [Checks](#checks)
  - [SMA001W1](#sma001w1)
  - [SMA001W3](#sma001w3)
  - [SMA001B1](#sma001b1)
  - [SMA001B3](#sma001b3)
  - [SMA002W1](#sma002w1)
  - [SMA002W2](#sma002w2)
  - [SMA002W3](#sma002w3)
- [Summary table](#summary-table)
- [Figure](#figure)

- **Paper:** The potential of elemental analysis in determining origin
  of sawfish rostra
- **DOI:** tbd
- **Written:** summary Oct 2026
- **Content:** plots: match_window over the anchor, Mg’s window, the
  build-out and the align_multi range, plus the summary table

Starts from 06a - 06d (`data/output/dtw_cache/`). The earlier version is
in `original/06e_dtw_rates.R`. For each set, `match_window()` with all
six elements is run over four stretches of the vertebra, each showing
the overlay and the warp path:

1.  **Anchor**: where the set’s anchor pair maps the same tooth to the
    same place (06c2).
2.  **Mg window**: where Mg places the whole, uncut tooth (06c); Sr’s if
    Mg isn’t in the run (see 06d).
3.  **Build-out**: what 06d’s build aligned: the whole stretch of
    vertebra its windows built out to, held to the stretch of tooth they
    matched (so the average ratio is fixed, too). The warp path lets the
    growth ratio change along it, as it did window to window in the
    build; the warping allowance is widened (`fixed_warp`) to give it
    room.
4.  **align_multi range**: the stretch `align_multi()` matched the tooth
    to, run on the whole records with all six elements.

Checks 1, 2 and 4 place the tooth themselves, anywhere on it, letting it
stretch by 2–8 µm of tooth per µm of vertebra: independent checks on the
build. The table at the end gives all four, with the build-out windows’
own statistics and `align_multi()` itself.

``` r
source('helpers/microHelp.R')
source('helpers/functions/fxn_mavs.R')
source('helpers/functions/fxn_dtw.R')
source('helpers/functions/fxn_warps.R')
load('data/output/dtw_cache/06a.RData') # objects from 06a
load('data/output/dtw_cache/06c2.RData') # sets, gap_um, down, anchors from 06c2
load('data/output/dtw_cache/06d.RData') # up, built (build-out windows), up_settings from 06d
```

``` r
check_settings <- list(ratio_range = c(2, 8), n_ratios = 24,
                       fixed_warp = 0.4)   # check 3's warping allowance (match_window's default is 0.3)
checks <- list()

# align_multi() on the whole records, all elements; its path in um
force_align <- function(n, k = 1000){
  set <- sets[[n]]
  invisible(capture.output(al <- align_multi(set$tooth, set$vert, mavs2, n = k, open_begin = TRUE,
                                             open_end = TRUE, gap_um = gap_um)))
  tt <- sort(set$tooth$true_dist)[round(seq(1, nrow(set$tooth), length.out = k))]
  vv <- sort(set$vert$true_dist)[round(seq(1, nrow(set$vert), length.out = k))]
  list(path = data.frame(set = n, vert = vv[al$index2], tooth = tt[al$index1]),
       norm_dist = al$normalizedDistance)
}
fa <- lapply(setNames(names(up), names(up)), force_align)   # the versions 06d built
```

``` r
# build-out windows (all versions, from 06d)
built$v_mid <- (built$v_from + built$v_to) / 2
built$t_mid <- (built$tooth_start + built$tooth_end) / 2
built_range <- function(n) c(min(built$v_from[built$set == n]), max(built$v_to[built$set == n]))
built_tooth <- function(n) c(min(built$tooth_start[built$set == n]), max(built$tooth_end[built$set == n]))

# Mg's window: where Mg places the whole (uncut) tooth in 06c's run (Sr's if Mg isn't in it)
mg_window <- function(n){
  r <- down[[n]]$runs
  for(el in c("Mg_Ca_mav", "Sr_Ca_mav")){
    x <- r[r$element == el & r$t_start == min(r$t_start) & r$t_end == max(r$t_end), ]
    if(nrow(x)) return(c(x$v_start, x$v_end))
  }
  NULL
}

# one match_window() over a stretch of vertebra, all elements, with its figure
preps <- lapply(setNames(names(up), names(up)), function(n)
  prep_records(sets[[n]]$tooth, sets[[n]]$vert, mavs2, gap_um = gap_um))
# tooth = a stretch of tooth to hold it to (one ratio, the tooth's length / the vert's; the warp
# path lets it change along the way); NULL = free, anywhere on the tooth
win_check <- function(n, label, v, tooth = NULL){
  fixed <- !is.null(tooth)
  cbind(set = n, check = label,
        match_window(sets[[n]]$tooth, sets[[n]]$vert, mavs2, v_from = v[1], win = diff(v), prep = preps[[n]], gap_um = gap_um,
                     ratio_range = if(fixed) rep(diff(tooth) / diff(v), 2) else check_settings$ratio_range,
                     n_ratios = if(fixed) 1 else check_settings$n_ratios,
                     tooth_search = if(fixed) tooth + c(-1, 1) else c(-Inf, Inf),   # +-1 um: rounding room
                     window_frac = if(fixed) check_settings$fixed_warp else 0.3,
                     plot = TRUE, bands = sets[[n]]$bands))
}
```

## Checks

``` r
for(n in names(up)){
  cat("\n\n###", n, "\n\n")
  ranges <- list("1 Anchor" = anchors[[n]], "2 Mg window" = mg_window(n),
                 "3 Build-out" = built_range(n), "4 align_multi range" = range(fa[[n]]$path$vert))
  for(k in names(ranges)){
    if(length(ranges[[k]]) != 2) next   # e.g. neither Mg nor Sr in the run
    cat("\n\n**", k, "** (vertebra ", round(ranges[[k]][1]), "-", round(ranges[[k]][2]), " um)\n\n", sep = "")
    checks[[length(checks) + 1]] <- win_check(n, k, ranges[[k]], if(k == "3 Build-out") built_tooth(n))
  }
}
```

### SMA001W1

**1 Anchor** (vertebra 5948-10517 um)

![](06f_dtw_plots_files/figure-gfm/checks-1.png)<!-- -->

**2 Mg window** (vertebra 2185-10258 um)

![](06f_dtw_plots_files/figure-gfm/checks-2.png)<!-- -->

**3 Build-out** (vertebra 5911-9061 um)

![](06f_dtw_plots_files/figure-gfm/checks-3.png)<!-- -->

**4 align_multi range** (vertebra 3658-5081 um)

![](06f_dtw_plots_files/figure-gfm/checks-4.png)<!-- -->

### SMA001W3

**1 Anchor** (vertebra 10079-10617 um)

![](06f_dtw_plots_files/figure-gfm/checks-5.png)<!-- -->

**2 Mg window** (vertebra 2185-10258 um)

![](06f_dtw_plots_files/figure-gfm/checks-6.png)<!-- -->

**3 Build-out** (vertebra 9711-10461 um)

![](06f_dtw_plots_files/figure-gfm/checks-7.png)<!-- -->

**4 align_multi range** (vertebra 3658-5081 um)

![](06f_dtw_plots_files/figure-gfm/checks-8.png)<!-- -->

### SMA001B1

**1 Anchor** (vertebra 5656-9264 um)

![](06f_dtw_plots_files/figure-gfm/checks-9.png)<!-- -->

**2 Mg window** (vertebra 3865-9799 um)

![](06f_dtw_plots_files/figure-gfm/checks-10.png)<!-- -->

**3 Build-out** (vertebra 5546-8296 um)

![](06f_dtw_plots_files/figure-gfm/checks-11.png)<!-- -->

**4 align_multi range** (vertebra 4087-9504 um)

![](06f_dtw_plots_files/figure-gfm/checks-12.png)<!-- -->

### SMA001B3

**1 Anchor** (vertebra 5877-9476 um)

![](06f_dtw_plots_files/figure-gfm/checks-13.png)<!-- -->

**2 Mg window** (vertebra 3865-9799 um)

![](06f_dtw_plots_files/figure-gfm/checks-14.png)<!-- -->

**3 Build-out** (vertebra 5746-9896 um)

![](06f_dtw_plots_files/figure-gfm/checks-15.png)<!-- -->

**4 align_multi range** (vertebra 4087-9504 um)

![](06f_dtw_plots_files/figure-gfm/checks-16.png)<!-- -->

### SMA002W1

**1 Anchor** (vertebra 3993-4572 um)

![](06f_dtw_plots_files/figure-gfm/checks-17.png)<!-- -->

**2 Mg window** (vertebra 3635-5899 um)

![](06f_dtw_plots_files/figure-gfm/checks-18.png)<!-- -->

**3 Build-out** (vertebra 1534-7884 um)

![](06f_dtw_plots_files/figure-gfm/checks-19.png)<!-- -->

**4 align_multi range** (vertebra 3931-4251 um)

![](06f_dtw_plots_files/figure-gfm/checks-20.png)<!-- -->

### SMA002W2

**1 Anchor** (vertebra 4168-5479 um)

![](06f_dtw_plots_files/figure-gfm/checks-21.png)<!-- -->

**2 Mg window** (vertebra 3635-5899 um)

![](06f_dtw_plots_files/figure-gfm/checks-22.png)<!-- -->

**3 Build-out** (vertebra 4134-6284 um)

![](06f_dtw_plots_files/figure-gfm/checks-23.png)<!-- -->

**4 align_multi range** (vertebra 3931-4251 um)

![](06f_dtw_plots_files/figure-gfm/checks-24.png)<!-- -->

### SMA002W3

**1 Anchor** (vertebra 1363-4275 um)

![](06f_dtw_plots_files/figure-gfm/checks-25.png)<!-- -->

**2 Mg window** (vertebra 3635-5899 um)

![](06f_dtw_plots_files/figure-gfm/checks-26.png)<!-- -->

**3 Build-out** (vertebra 334-2084 um)

![](06f_dtw_plots_files/figure-gfm/checks-27.png)<!-- -->

**4 align_multi range** (vertebra 3931-4251 um)

![](06f_dtw_plots_files/figure-gfm/checks-28.png)<!-- -->

``` r
checks <- do.call(rbind, checks)
```

## Summary table

``` r
summary_tab <- do.call(rbind, lapply(names(up), function(n){
  b <- built[built$set == n, ]; p <- fa[[n]]$path; ck <- checks[checks$set == n, ]
  fit <- lm(t_mid ~ v_mid, data = b)
  rbind(
    data.frame(set = n, result = ck$check, windows = "1",
               vert_from = ck$v_from, vert_to = ck$v_to, tooth_from = ck$tooth_start, tooth_to = ck$tooth_end,
               R = ck$ratio, R_r2 = NA, fit = ck$cost, disagree_pct = ck$disagree_pct),
    data.frame(set = n, result = "Build-out windows", windows = sprintf("%d / %d", nrow(b), nrow(up[[n]])),
               vert_from = min(b$v_from), vert_to = max(b$v_to), tooth_from = min(b$tooth_start), tooth_to = max(b$tooth_end),
               R = coef(fit)[2], R_r2 = summary(fit)$r.squared, fit = median(b$cost), disagree_pct = median(b$disagree_pct)),
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
| SMA001W1 | 1 Anchor | 1 | 5948 | 10517 | 4364 | 26613 | 4.87 | NA | 0.36 | 23.72 |
| SMA001W1 | 2 Mg window | 1 | 2185 | 10258 | 10464 | 35036 | 3.04 | NA | 0.33 | 21.64 |
| SMA001W1 | 3 Build-out | 1 | 5911 | 9061 | 40299 | 50199 | 3.14 | NA | 0.60 | 44.62 |
| SMA001W1 | 4 align_multi range | 1 | 3658 | 5081 | 24064 | 32848 | 6.17 | NA | 0.26 | 21.78 |
| SMA001W1 | Build-out windows | 13 / 13 | 5911 | 9061 | 40300 | 50200 | 2.76 | 0.74 | 0.38 | 7.32 |
| SMA001W1 | align_multi |  | 3658 | 5081 | 864 | 54832 | 37.93 | NA | 2.10 | NA |
| SMA001W3 | 1 Anchor | 1 | 10079 | 10617 | 8364 | 9440 | 2.00 | NA | 0.51 | 27.27 |
| SMA001W3 | 2 Mg window | 1 | 2185 | 10258 | 10464 | 35036 | 3.04 | NA | 0.33 | 21.64 |
| SMA001W3 | 3 Build-out | 1 | 9711 | 10461 | 46899 | 51149 | 5.67 | NA | 0.86 | 44.46 |
| SMA001W3 | 4 align_multi range | 1 | 3658 | 5081 | 24064 | 32848 | 6.17 | NA | 0.26 | 21.78 |
| SMA001W3 | Build-out windows | 1 / 1 | 9711 | 10461 | 46900 | 51150 | NA | 0.00 | 0.33 | 7.70 |
| SMA001W3 | align_multi |  | 3658 | 5081 | 864 | 54832 | 37.93 | NA | 2.10 | NA |
| SMA001B1 | 1 Anchor | 1 | 5656 | 9264 | 17502 | 45424 | 7.74 | NA | 0.36 | 29.41 |
| SMA001B1 | 2 Mg window | 1 | 3865 | 9799 | 12302 | 52033 | 6.70 | NA | 0.31 | 28.36 |
| SMA001B1 | 3 Build-out | 1 | 5546 | 8296 | 34601 | 48449 | 5.04 | NA | 0.47 | 29.17 |
| SMA001B1 | 4 align_multi range | 1 | 4087 | 9504 | 12102 | 52616 | 7.48 | NA | 0.32 | 27.83 |
| SMA001B1 | Build-out windows | 11 / 11 | 5546 | 8296 | 34602 | 48450 | 5.62 | 0.93 | 0.52 | 8.80 |
| SMA001B1 | align_multi |  | 4087 | 9504 | 12002 | 52836 | 7.54 | NA | 2.03 | NA |
| SMA001B3 | 1 Anchor | 1 | 5877 | 9476 | 19702 | 45677 | 7.22 | NA | 0.36 | 28.64 |
| SMA001B3 | 2 Mg window | 1 | 3865 | 9799 | 12302 | 52033 | 6.70 | NA | 0.31 | 28.36 |
| SMA001B3 | 3 Build-out | 1 | 5746 | 9896 | 34601 | 46249 | 2.81 | NA | 0.30 | 35.93 |
| SMA001B3 | 4 align_multi range | 1 | 4087 | 9504 | 12102 | 52616 | 7.48 | NA | 0.32 | 27.83 |
| SMA001B3 | Build-out windows | 18 / 18 | 5746 | 9896 | 34602 | 46250 | 1.86 | 0.70 | 0.38 | 7.09 |
| SMA001B3 | align_multi |  | 4087 | 9504 | 12002 | 52836 | 7.54 | NA | 2.03 | NA |
| SMA002W1 | 1 Anchor | 1 | 3993 | 4572 | 14156 | 18183 | 6.96 | NA | 0.28 | 26.79 |
| SMA002W1 | 2 Mg window | 1 | 3635 | 5899 | 14056 | 20946 | 3.04 | NA | 0.32 | 21.67 |
| SMA002W1 | 3 Build-out | 1 | 1534 | 7884 | 4955 | 38949 | 5.35 | NA | 0.36 | 43.66 |
| SMA002W1 | 4 align_multi range | 1 | 3931 | 4251 | 19856 | 21412 | 4.87 | NA | 0.24 | 28.08 |
| SMA002W1 | Build-out windows | 29 / 29 | 1534 | 7884 | 4956 | 38950 | 5.92 | 0.93 | 0.44 | 6.50 |
| SMA002W1 | align_multi |  | 3931 | 4251 | 156 | 42969 | 133.92 | NA | 2.20 | NA |
| SMA002W2 | 1 Anchor | 1 | 4168 | 5479 | 27756 | 31062 | 2.52 | NA | 0.38 | 25.54 |
| SMA002W2 | 2 Mg window | 1 | 3635 | 5899 | 14056 | 20946 | 3.04 | NA | 0.32 | 21.67 |
| SMA002W2 | 3 Build-out | 1 | 4134 | 6284 | 26299 | 41149 | 6.91 | NA | 0.48 | 43.41 |
| SMA002W2 | 4 align_multi range | 1 | 3931 | 4251 | 19856 | 21412 | 4.87 | NA | 0.24 | 28.08 |
| SMA002W2 | Build-out windows | 8 / 8 | 4134 | 6284 | 26300 | 41150 | 4.24 | 0.50 | 0.17 | 7.07 |
| SMA002W2 | align_multi |  | 3931 | 4251 | 156 | 42969 | 133.92 | NA | 2.20 | NA |
| SMA002W3 | 1 Anchor | 1 | 1363 | 4275 | 16156 | 26537 | 3.57 | NA | 0.40 | 23.17 |
| SMA002W3 | 2 Mg window | 1 | 3635 | 5899 | 14056 | 20946 | 3.04 | NA | 0.32 | 21.67 |
| SMA002W3 | 3 Build-out | 1 | 334 | 2084 | 1555 | 16899 | 8.77 | NA | 0.62 | 53.86 |
| SMA002W3 | 4 align_multi range | 1 | 3931 | 4251 | 19856 | 21412 | 4.87 | NA | 0.24 | 28.08 |
| SMA002W3 | Build-out windows | 6 / 6 | 334 | 2084 | 1556 | 16900 | 14.10 | 0.81 | 0.74 | 12.47 |
| SMA002W3 | align_multi |  | 3931 | 4251 | 156 | 42969 | 133.92 | NA | 2.20 | NA |

- **result**: 1–4 = the `match_window()` checks above; **Build-out
  windows** = 06d’s 500 µm windows taken together; **align_multi** = the
  DTW path over the whole records.
- **windows**: for the build-out, windows that matched, of all windows.
- **vert_from, vert_to, tooth_from, tooth_to**: where the result starts
  and ends on the vertebra and the tooth (µm).
- **R**: growth ratio, µm of tooth per µm of vertebra. For the build-out
  windows, the slope of matched tooth position on vertebra position, and
  **R_r2** how straight that is (1 = constant growth ratio).
- **fit**: DTW fit (normalised distance, z-score units; lower = closer
  match). For the build-out windows, the median over them.
- **disagree_pct**: how far apart the elements place the matched tooth
  piece, as a share of the piece.

## Figure

``` r
bands <- do.call(rbind, lapply(names(up), function(n) data.frame(set = n, band = sets[[n]]$bands)))
anch  <- do.call(rbind, lapply(names(up), function(n) data.frame(set = n, from = anchors[[n]][1], to = anchors[[n]][2])))
paths <- do.call(rbind, lapply(fa, `[[`, "path"))
p <- ggplot() +
  geom_rect(data = anch, aes(xmin = from, xmax = to, ymin = -Inf, ymax = Inf), fill = "grey88") +
  geom_vline(data = bands, aes(xintercept = band), colour = "brown", linewidth = 0.3) +
  geom_line(data = paths, aes(vert, tooth), colour = "#1b9e77", linewidth = 0.8) +
  geom_segment(data = checks[checks$check == "3 Build-out", ], aes(x = v_from, xend = v_to, y = tooth_start, yend = tooth_end),
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
         "far apart the elements place it; size = how much the vertebra window varies). Orange line = check 3: the",
         "built-out stretch of vertebra held to the tooth the build matched (slope = average growth ratio). Green line = align_multi()",
         "path over the whole records. Grey band = anchor; brown lines = vertebral bands."), 115), collapse = "\n"))
p
```

![](06f_dtw_plots_files/figure-gfm/figure-1.png)<!-- -->

``` r
dir.create('data/output/dtw_final', showWarnings = FALSE)
dir.create('figures/dtw', showWarnings = FALSE)
write.csv(summary_tab, 'data/output/dtw_final/summary.csv', row.names = FALSE)
ggsave('figures/dtw/summary.png', p, width = 10, height = 9, dpi = 300)
save(summary_tab, checks, built, fa, file = 'data/output/dtw_cache/06f.RData')
```

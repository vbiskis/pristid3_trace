06f · DTW warping and misfit
================
Nikki Biskis
Rendered 27 Oct 2026 02:11

- [Growth ratio along the vertebra](#growth-ratio-along-the-vertebra)
- [Misfit and warp of each element](#misfit-and-warp-of-each-element)
- [What goes with poor alignment?](#what-goes-with-poor-alignment)

- **Paper:** The potential of elemental analysis in determining origin
  of sawfish rostra
- **DOI:** tbd
- **Written:** 30 Aug 2026 (warping and misfit Oct 2026)
- **Content:** dynamic time warping VI: how the warping and the misfit
  vary along the vertebra, and what goes with poor alignment

Starts from 06a, 06c, 06d and 06e (`data/output/dtw_cache/`); uses 06d’s
build-out windows. The earlier version is in `original/06f_dtw_plots.R`.
For each window:

- **Growth ratio**: µm of tooth matched per µm of vertebra (`ratio`).
- **Warp** of each element: how far its path through the window departs
  from an even pace, as % of the window (0 = no warping).
- **Misfit** of each element: its DTW fit in the window (normalised
  distance; higher = poorer match), and how far apart the elements place
  the tooth piece (`disagree_pct`).
- **Possible causes**: how much the vertebra window varies (`signal`;
  flat stretches are hard to place), the share of the matched tooth
  piece inside holes, median In:Ca in the tooth piece and in the
  vertebra window (contamination), and the distance from the window to
  the nearest vertebral band.

``` r
source('helpers/microHelp.R')
source('helpers/functions/fxn_mavs.R')
source('helpers/functions/fxn_dtw.R')    # also find_gaps(), in_gaps()
load('data/output/dtw_cache/06a.RData')  # objects from 06a
load('data/output/dtw_cache/06c.RData')  # sets, gap_um from 06c
load('data/output/dtw_cache/06d.RData')  # up_settings from 06d
load('data/output/dtw_cache/06e.RData')  # built (matched build-out windows) from 06e
```

``` r
# candidate causes of misfit, per window
built <- do.call(rbind, lapply(split(built, built$set), function(b){
  set   <- sets[[b$set[1]]]
  holes <- find_gaps(set$tooth$true_dist, gap_um)
  b$hole_pct <- sapply(seq_len(nrow(b)), function(i){
    x <- seq(b$tooth_start[i], b$tooth_end[i], length.out = 200)
    100 * mean(in_gaps(x, holes))
  })
  b$In_tooth <- sapply(seq_len(nrow(b)), function(i)
    median(set$tooth$In_Ca[set$tooth$true_dist >= b$tooth_start[i] & set$tooth$true_dist <= b$tooth_end[i]]))
  b$In_vert <- sapply(seq_len(nrow(b)), function(i)
    median(set$vert$In_Ca[set$vert$true_dist >= b$v_from[i] & set$vert$true_dist <= b$v_to[i]]))
  b$band_dist <- sapply(b$v_mid, function(v) min(abs(v - set$bands)))
  b
}))

els  <- sub("_Ca.*$", "", up_settings$elements)
long <- do.call(rbind, lapply(up_settings$elements, function(el)
  data.frame(set = built$set, v_mid = built$v_mid, element = sub("_Ca.*$", "", el),
             cost = built[[paste0("cost_", el)]], warp = built[[paste0("warp_", el)]])))
bands <- do.call(rbind, lapply(names(sets), function(n) data.frame(set = n, band = sets[[n]]$bands)))
```

## Growth ratio along the vertebra

``` r
R_all <- do.call(rbind, lapply(split(built, built$set), function(b)
  data.frame(set = b$set[1], R = coef(lm(t_mid ~ v_mid, data = b))[2])))
ggplot(built, aes(v_mid, ratio)) +
  geom_vline(data = bands, aes(xintercept = band), colour = "brown", linewidth = 0.3) +
  geom_hline(data = R_all, aes(yintercept = R), linetype = "dashed") +
  geom_line(colour = "grey50") + geom_point(aes(colour = disagree_pct)) +
  scale_colour_viridis_c(name = "Elements\ndisagree (%)") +
  facet_wrap(~ set, ncol = 1, scales = "free_x") +
  theme_bw() + theme(plot.caption = element_text(hjust = 0)) +
  labs(x = "Vertebra (um)", y = "Growth ratio (um tooth per um vertebra)",
       caption = paste("Growth ratio chosen by each build-out window. Dashed = the build-out's overall growth ratio",
                       "\n(slope of matched tooth on vertebra position). Brown lines = vertebral bands."))
```

![](06f_dtw_plots_files/figure-gfm/ratio-1.png)<!-- -->

## Misfit and warp of each element

``` r
ggplot(long, aes(v_mid, cost, colour = element)) +
  geom_vline(data = bands, aes(xintercept = band), colour = "brown", linewidth = 0.3) +
  geom_line() +
  facet_wrap(~ set, ncol = 1, scales = "free") +
  theme_bw() + theme(plot.caption = element_text(hjust = 0)) +
  labs(x = "Vertebra (um)", y = "Misfit (DTW normalised distance)", colour = NULL,
       caption = paste("Each element's fit between the vertebra window and its matched tooth piece; higher = poorer match.",
                       "\nBrown lines = vertebral bands."))
```

![](06f_dtw_plots_files/figure-gfm/misfit-1.png)<!-- -->

``` r
ggplot(long, aes(v_mid, warp, colour = element)) +
  geom_vline(data = bands, aes(xintercept = band), colour = "brown", linewidth = 0.3) +
  geom_line() +
  facet_wrap(~ set, ncol = 1, scales = "free") +
  theme_bw() + theme(plot.caption = element_text(hjust = 0)) +
  labs(x = "Vertebra (um)", y = "Warp (% of the window)", colour = NULL,
       caption = paste("How far each element's path through the window departs from an even pace (straight diagonal),",
                       "\nas % of the window; 0 = no warping. Brown lines = vertebral bands."))
```

![](06f_dtw_plots_files/figure-gfm/warp-1.png)<!-- -->

## What goes with poor alignment?

Spearman correlation of each window’s misfit (median over the elements)
and of how far apart the elements place the tooth piece, with each
possible cause.

``` r
built$misfit <- apply(built[paste0("cost_", up_settings$elements)], 1, median)
causes <- c(signal = "Vertebra signal", hole_pct = "Holes in tooth piece (%)", In_tooth = "In:Ca, tooth piece",
            In_vert = "In:Ca, vertebra window", band_dist = "Distance to nearest band (um)")
cause_tab <- do.call(rbind, lapply(split(built, built$set), function(b) do.call(rbind, lapply(names(causes), function(cz){
  ct <- function(y) suppressWarnings(cor.test(b[[cz]], b[[y]], method = "spearman"))
  m <- ct("misfit"); d <- ct("disagree_pct")
  data.frame(set = b$set[1], cause = causes[[cz]],
             rho_misfit = round(m$estimate, 2), p_misfit = signif(m$p.value, 2),
             rho_disagree = round(d$estimate, 2), p_disagree = signif(d$p.value, 2))
}))))
rownames(cause_tab) <- NULL
knitr::kable(cause_tab)
```

| set | cause | rho_misfit | p_misfit | rho_disagree | p_disagree |
|:---|:---|---:|---:|---:|---:|
| SMA001B | Vertebra signal | 0.57 | 7.8e-05 | 0.39 | 9.0e-03 |
| SMA001B | Holes in tooth piece (%) | 0.25 | 1.1e-01 | 0.00 | 9.8e-01 |
| SMA001B | In:Ca, tooth piece | 0.56 | 8.0e-05 | 0.66 | 1.0e-06 |
| SMA001B | In:Ca, vertebra window | 0.06 | 6.9e-01 | 0.70 | 4.0e-07 |
| SMA001B | Distance to nearest band (um) | -0.53 | 2.3e-04 | 0.04 | 7.9e-01 |
| SMA001W | Vertebra signal | 0.32 | 2.8e-02 | 0.40 | 5.3e-03 |
| SMA001W | Holes in tooth piece (%) | 0.48 | 5.0e-04 | 0.38 | 7.0e-03 |
| SMA001W | In:Ca, tooth piece | 0.69 | 0.0e+00 | 0.42 | 2.7e-03 |
| SMA001W | In:Ca, vertebra window | 0.14 | 3.5e-01 | 0.12 | 4.1e-01 |
| SMA001W | Distance to nearest band (um) | -0.34 | 2.0e-02 | -0.27 | 6.9e-02 |
| SMA002W | Vertebra signal | 0.37 | 2.8e-02 | 0.57 | 4.0e-04 |
| SMA002W | Holes in tooth piece (%) | -0.03 | 8.7e-01 | -0.34 | 4.4e-02 |
| SMA002W | In:Ca, tooth piece | 0.83 | 0.0e+00 | 0.53 | 9.2e-04 |
| SMA002W | In:Ca, vertebra window | -0.36 | 3.3e-02 | -0.48 | 3.6e-03 |
| SMA002W | Distance to nearest band (um) | -0.29 | 8.6e-02 | -0.36 | 3.3e-02 |

- **rho_misfit, p_misfit**: Spearman correlation (and p-value) between
  the cause and the window’s misfit (median over elements).
- **rho_disagree, p_disagree**: the same with how far apart the elements
  place the tooth piece.
- Reading it: a positive rho with **In:Ca** points to contamination;
  with **holes**, to missing tooth data; a negative rho with **vertebra
  signal** means flat stretches are hard to place; a negative rho with
  **distance to nearest band** means poorer matches close to bands (slow
  growth compresses the vertebra there). NA = the cause doesn’t vary in
  that set (e.g. no holes).

``` r
cl <- do.call(rbind, lapply(names(causes), function(cz)
  data.frame(set = built$set, cause = causes[[cz]], value = built[[cz]], misfit = built$misfit)))
ggplot(cl, aes(value, misfit, colour = set)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
  facet_wrap(~ cause, scales = "free_x", ncol = 2) +
  theme_bw() + theme(plot.caption = element_text(hjust = 0)) +
  labs(x = NULL, y = "Misfit (median over elements)", colour = NULL,
       caption = "Each point is one build-out window; lines = linear trend per set.")
```

![](06f_dtw_plots_files/figure-gfm/causes-plot-1.png)<!-- -->

``` r
dir.create('data/output/dtw_final', showWarnings = FALSE)
write.csv(cause_tab, 'data/output/dtw_final/misfit_causes.csv', row.names = FALSE)
write.csv(built, 'data/output/dtw_final/buildout_windows.csv', row.names = FALSE)
```

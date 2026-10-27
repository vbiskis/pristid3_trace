06a2 · Smoothing test
================
Nikki Biskis
Rendered 27 Oct 2026 21:03

- [Results](#results)

- **Paper:** The potential of elemental analysis in determining origin
  of sawfish rostra
- **DOI:** tbd
- **Written:** Oct 2026
- **Content:** how wide the rolling median should be, tested before the
  DTW

The DTW (06b–06f) uses only 06a’s rolling medians (`add_mavs()`); there
is no further smoothing. This checks the medians’ width \[a la Hegg et
al. (2018)\]: run a breakpoint analysis at several smoothing widths,
from the raw data up, and see whether the result depends on the
smoothing.

Points are every 0.288 s: 2.9 µm apart on the vertebrae (10 µm/s) and
8.6 µm on the teeth (30 µm/s). The spot is 50 µm, so ~17 vertebra points
and ~6 tooth points fall under one spot. The widths tried are matched in
spot widths across the two records:

\| raw \| ~1 spot \| ~2 spot \| ~2.6 spot \| ~3.5 spots \|

\|—\|—\| - - - \| - - - \| - - - \| - - - \| \| V \| 1 \| 17 \| 31 \| 45
\| 61 \| \| T \| 1 \| 7 \| 11 \| 15 \| 21 \|

For each record, element and width:

1.  the raw element:Ca is smoothed with `add_mavs()` (as in 06a);
2.  one value is taken per 50 µm (one per spot width), every width uses
    the same stretch (the widest median’s ends dropped);
3.  breakpoints = changes in mean and variance; found with
    `changepoint::cpt.meanvar()` (PELT, MBIC penalty, at least 4 values
    = 200 µm between breakpoints)ß

**agree_next** is the share of a width’s breakpoints that are still
found (within 2 spots, 100 µm) at the next wider width. When it stays
high, the smoothing isn’t moving the breakpoints.

``` r
source('helpers/microHelp.R')
source('helpers/functions/fxn_mavs.R')
library(changepoint)
load('data/output/dtw_cache/06a.RData') # objects from 06a
```

``` r
els    <- c("Sr_Ca", "Ba_Ca", "Li_Ca", "Mg_Ca", "Mn_Ca", "Zn_Ca")
spot   <- 50                                                     # laser spot (um)
widths <- list(vert = c(1, 17, 31, 45, 61), tooth = c(1, 7, 11, 15, 21))
gaps   <- c(vert = Inf, tooth = 350)                             # holes, as in 06a
recs <- list(SMA001W = list(tooth = SMA001T_filt,  vert = SMA001V_filt),
             SMA001B = list(tooth = SMA001TB_filt, vert = SMA001VB_filt),
             SMA002W = list(tooth = SMA002T_filt,  vert = SMA002V_filt))

# breakpoints (um) in one element of one record, smoothed with a w-point rolling median
breaks_at <- function(d, el, w, gap){
  d <- d[order(d$true_dist), c("true_dist", el)]
  if(w > 1) d[[el]] <- add_mavs(d, el, window_size = w, gap_um = gap)[[paste0(el, "_mav")]]
  half <- (max(unlist(widths)) - 1) / 2                          # same stretch for every width
  d <- d[(half + 1):(nrow(d) - half), ]
  g <- seq(min(d$true_dist), max(d$true_dist), by = spot)        # one value per spot width
  i <- findInterval(g, d$true_dist, all.inside = TRUE)
  i <- ifelse(abs(d$true_dist[i + 1] - g) < abs(d$true_dist[i] - g), i + 1, i)
  keep <- abs(d$true_dist[i] - g) <= spot / 2 & is.finite(d[[el]][i])   # none inside holes
  x <- d[[el]][i][keep]
  g[keep][cpts(cpt.meanvar(x, method = "PELT", penalty = "MBIC", minseglen = 4))]
}

runs <- expand.grid(step = 1:5, element = els, record = c("vert", "tooth"), fish = names(recs),
                    stringsAsFactors = FALSE)
runs$width <- mapply(function(s, k) widths[[s]][k], runs$record, runs$step)
runs$spacing <- ifelse(runs$record == "vert", 10, 30) * 0.288    # um between readings
runs$spots <- round(runs$width * runs$spacing / spot, 1)
brk <- lapply(seq_len(nrow(runs)), function(j)
  breaks_at(recs[[runs$fish[j]]][[runs$record[j]]], runs$element[j], runs$width[j], gaps[[runs$record[j]]]))
runs$n_breaks <- lengths(brk)

# share of each width's breakpoints still found (within 2 spots) at the next wider width
runs$agree_next <- sapply(seq_len(nrow(runs)), function(j){
  if(runs$step[j] == 5 || !length(brk[[j]])) return(NA)
  mean(sapply(brk[[j]], function(x) any(abs(brk[[j + 1]] - x) <= 2 * spot)))
})
```

## Results

``` r
summ <- aggregate(cbind(n_breaks, agree_next) ~ record + step + width + spots, runs,
                  function(x) round(median(x), 2), na.action = na.pass)
summ <- summ[order(summ$record, summ$step), ]
rownames(summ) <- NULL
knitr::kable(summ)
```

| record | step | width | spots | n_breaks | agree_next |
|:-------|-----:|------:|------:|---------:|-----------:|
| tooth  |    1 |     1 |   0.2 |     20.5 |       0.85 |
| tooth  |    2 |     7 |   1.2 |     32.0 |       0.95 |
| tooth  |    3 |    11 |   1.9 |     38.0 |       0.91 |
| tooth  |    4 |    15 |   2.6 |     39.0 |       0.90 |
| tooth  |    5 |    21 |   3.6 |     40.0 |         NA |
| vert   |    1 |     1 |   0.1 |      6.0 |       0.75 |
| vert   |    2 |    17 |   1.0 |      7.0 |       1.00 |
| vert   |    3 |    31 |   1.8 |     10.0 |       0.92 |
| vert   |    4 |    45 |   2.6 |     11.0 |       0.92 |
| vert   |    5 |    61 |   3.5 |     12.5 |         NA |

- **width**: rolling median, in points; **spots**: the same in 50 µm
  spot widths.
- **n_breaks**, **agree_next**: medians over the six elements and three
  fish.

``` r
p <- ggplot(runs[!is.na(runs$agree_next), ], aes(spots, agree_next, colour = element)) +
  geom_line() + geom_point() +
  facet_grid(record ~ fish) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(plot.caption = element_text(hjust = 0)) +
  labs(x = "Rolling median width (spot widths, 50 um)", y = "Breakpoints still found at the next width",
       caption = paste(strwrap(paste(
         "For each fish (columns), record (rows) and element: the share of breakpoints (changepoint::cpt.meanvar,",
         "PELT) found with one rolling median width that are still found, within 100 um, with the next wider one.",
         "Values are taken once per 50 um spot width. High and flat = the smoothing no longer moves the breakpoints."),
         115), collapse = "\n"))
p
```

![](06a2_smooth_test_files/figure-gfm/figure-1.png)<!-- -->

``` r
dir.create('figures/dtw', showWarnings = FALSE)
ggsave('figures/dtw/smooth_test.png', p, width = 10, height = 7, dpi = 300)
save(runs, summ, file = 'data/output/dtw_cache/06a2.RData')
```

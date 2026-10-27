06f · DTW comparison
================
Nikki Biskis
Rendered 27 Oct 2026 22:10

- [Build-out by version](#build-out-by-version)
- [Placements](#placements)
  - [SMA001W](#sma001w)
  - [SMA001B](#sma001b)
  - [SMA002W](#sma002w)
- [Tables](#tables)
  - [SMA001W](#sma001w-1)
  - [SMA001B](#sma001b-1)
  - [SMA002W](#sma002w-1)

- **Paper:** The potential of elemental analysis in determining origin
  of sawfish rostra
- **DOI:** tbd
- **Written:** summary Oct 2026
- **Content:** where each method puts the tooth on the vertebra

Starts from 06a, 06c, 06d and 06e (`data/output/dtw_cache/`).

I’;; make it prettier eventuall! For now a lot of base R.

Each specimens alignment test:

- **A · 06c cut**: the tooth is cut at its start and end (06c). An
  element that tracks growth correctly starts later on the vertebra when
  the tooth’s start is cut, and ends earlier when its end is cut; one
  that doesn’t (e.g. always ending at the vertebra’s end) is just
  warping to fit. Elements whose placement moves the right way for both
  (Spearman ρ over all cuts above `rho_min`) are kept, and their
  placements of the whole, uncut tooth are averaged.
- **B · 06d build**: after identifying a solid anchor in 6c,
  systematically adding to the pattern
- **C · 06d free**: that same stretch of vertebra from B, but with the
  tooth placed freely by all five final elements (06d step 2), an
  independent check on B.
- **D · 06e mix**: `align_multi()` on the whole records with the largest
  mix of elements that keeps the placement within `shift_max` of Mg
  alone (06e’s shift_um), best fit at that size.
- **Mean**: the average of A–D (versions averaged within B and C first).

``` r
source('helpers/microHelp.R')
source('helpers/functions/fxn_mavs.R')
source('helpers/functions/fxn_dtw.R')
library(patchwork)
load('data/output/dtw_cache/06a.RData')  # objects from 06a
load('data/output/dtw_cache/06c.RData')  # fish (records, bands), stacks (cut runs and paths), gap_um from 06c
load('data/output/dtw_cache/06d.RData')  # anchors, anchor_pairs, up, built, final from 06d
load('data/output/dtw_cache/06e.RData')  # combo_tab from 06e
```

``` r
cmp_settings <- list(rho_min   = 0,             # A: an element is kept if both cut responses are above this
                     shift_max = 500,           # D: largest mix within this (um) of Mg alone
                     main      = "Mg_Ca_mav")   # element laid on the vertebra in the alignment panels
fish_of <- function(v) sub("[0-9]+$", "", v)
els_lab <- function(e) paste(sub("_Ca.*$", "", e), collapse = "+")
built$v_mid <- (built$v_from + built$v_to) / 2
built$t_mid <- (built$tooth_start + built$tooth_end) / 2
built$fish  <- fish_of(built$set)
built$version <- sprintf("%s (%s)", built$set, sapply(built$set, function(n) els_lab(anchor_pairs[[n]])))
```

## Build-out by version

``` r
anch <- do.call(rbind, lapply(names(up), function(n) data.frame(fish = fish_of(n), from = anchors[[n]][1],
  to = anchors[[n]][2], version = sprintf("%s (%s)", n, els_lab(anchor_pairs[[n]])))))
bands <- do.call(rbind, lapply(names(fish), function(f) data.frame(fish = f, band = fish[[f]]$bands)))
p_ver <- ggplot() +
  geom_rect(data = anch, aes(xmin = from, xmax = to, ymin = -Inf, ymax = Inf, fill = version), alpha = 0.12) +
  geom_vline(data = bands, aes(xintercept = band), colour = "brown", linewidth = 0.3) +
  geom_line(data = built, aes(v_mid, t_mid, colour = version), linewidth = 0.6) +
  geom_point(data = built, aes(v_mid, t_mid, colour = version, size = signal)) +
  scale_size_continuous(range = c(0.5, 3), name = "Signal") +
  facet_wrap(~ fish, ncol = 1, scales = "free") +
  theme_bw() +
  theme(plot.caption = element_text(hjust = 0)) +
  labs(x = "Vertebra (um)", y = "Tooth (um)", colour = "Version (anchor pair)", fill = "Version (anchor pair)",
       caption = paste(strwrap(paste(
         "Each specimen versions overlaid. Points = 06d build-out windows (x = vertebra window centre, y = centre",
         "of the matched tooth piece; size = how much the vertebra window varies), joined in order; shaded = each",
         "version's anchor; colour = version and the element pair it was anchored with. Brown lines = vertebral bands."),
         115), collapse = "\n"))
p_ver
```

![](06f_dtw_plots_files/figure-gfm/versions-1.png)<!-- -->

## Placements

``` r
# A: each element's response to cutting the tooth, and where it puts the whole (uncut) tooth
rho <- function(x, y) if(sd(x) == 0 || sd(y) == 0) NA else cor(x, y, method = "spearman")
cut_tab <- do.call(rbind, lapply(names(stacks), function(f){
  r <- stacks[[f]]$runs
  do.call(rbind, lapply(split(r, r$element), function(d){
    u <- d[d$t_start == min(d$t_start) & d$t_end == max(d$t_end), ]
    data.frame(fish = f, element = d$element[1], rho_start = rho(d$t_start, d$v_start), rho_end = rho(d$t_end, d$v_end),
               vert_from = u$v_start, vert_to = u$v_end, tooth_from = u$t_start, tooth_to = u$t_end, cut = u$cut)
  }))
}))
cut_tab$kept <- with(cut_tab, !is.na(rho_start) & !is.na(rho_end) &
                       rho_start > cmp_settings$rho_min & rho_end > cmp_settings$rho_min)

# each element's warp path for the uncut tooth, and the mean of the kept ones
el_paths <- do.call(rbind, lapply(split(cut_tab, seq_len(nrow(cut_tab))), function(x){
  p <- stacks[[x$fish]]$paths
  p <- p[p$cut == x$cut & p$element == x$element, c("tooth", "vert")]
  cbind(fish = x$fish, element = sub("_Ca.*$", "", x$element), kept = x$kept, p)
}))
mean_path <- function(f){
  k <- el_paths[el_paths$fish == f & el_paths$kept, ]
  if(nrow(k) == 0) return(NULL)
  tg <- seq(min(k$tooth), max(k$tooth), length.out = 500)
  v <- sapply(split(k, k$element), function(e) approx(e$tooth, e$vert, xout = tg, ties = mean)$y)
  data.frame(tooth = tg, vert = rowMeans(as.matrix(v), na.rm = TRUE))
}

# D: the largest mix within shift_max of Mg alone, best fit at that size; its align_multi() path
best_mix <- do.call(rbind, lapply(split(combo_tab, combo_tab$set), function(d){
  d <- d[d$shift_um <= cmp_settings$shift_max, ]
  d <- d[d$n_elements == max(d$n_elements), ]
  d[which.min(d$fit), ]
}))
mix_path <- function(f, k = 1000){
  set <- fish[[f]]; els <- paste0(strsplit(best_mix$elements[best_mix$set == f], "\\+")[[1]], "_Ca_mav")
  invisible(capture.output(al <- align_multi(set$tooth, set$vert, els, n = k, open_begin = TRUE,
                                             open_end = TRUE, gap_um = gap_um, plot = FALSE)))
  tt <- sort(set$tooth$true_dist)[round(seq(1, nrow(set$tooth), length.out = k))]
  vv <- sort(set$vert$true_dist)[round(seq(1, nrow(set$vert), length.out = k))]
  v  <- tapply(vv[al$index2], al$index1, mean)
  data.frame(tooth = tt[as.integer(names(v))], vert = as.numeric(v))
}

# every placement as a tooth -> vertebra path, one per row of the summary
paths <- list(); place <- list()
add <- function(f, method, label, p){
  p <- p[order(p$tooth), ]
  paths[[length(paths) + 1]] <<- cbind(fish = f, method = method, label = label, p)
  place[[length(place) + 1]] <<- data.frame(fish = f, method = method, label = label,
    vert_from = min(p$vert), vert_to = max(p$vert), tooth_from = min(p$tooth), tooth_to = max(p$tooth))
}
for(f in names(fish)){
  k <- cut_tab[cut_tab$fish == f & cut_tab$kept, ]
  if(nrow(k)) add(f, "A · 06c cut", sprintf("A · 06c cut (%s)", els_lab(k$element)), mean_path(f))
  for(n in names(up)[fish_of(names(up)) == f]){
    b <- built[built$set == n, ]; b <- b[order(b$v_from), ]
    add(f, "B · 06d build", sprintf("B · 06d build %s (%s)", n, els_lab(anchor_pairs[[n]])),
        data.frame(vert  = c(b$v_from[1], b$v_mid, b$v_to[nrow(b)]),
                   tooth = cummax(c(b$tooth_start[1], b$t_mid, b$tooth_end[nrow(b)]))))
  }
  for(n in names(up)[fish_of(names(up)) == f]){
    x <- final[final$set == n, ]
    add(f, "C · 06d free", sprintf("C · 06d free %s (%s)", n, els_lab(up_settings$final_elements)),
        data.frame(vert = c(x$v_from, x$v_to), tooth = c(x$tooth_start, x$tooth_end)))
  }
  add(f, "D · 06e mix", sprintf("D · 06e mix (%s)", best_mix$elements[best_mix$set == f]), mix_path(f))
}
paths <- do.call(rbind, paths); place <- do.call(rbind, place)

# mean: versions averaged within each method, then the methods
num <- c("vert_from", "vert_to", "tooth_from", "tooth_to")
by_m <- aggregate(place[num], place[c("fish", "method")], mean)
mean_row <- aggregate(by_m[num], by_m["fish"], mean)
place <- rbind(place, cbind(mean_row["fish"], method = "Mean", label = "Mean of A-D", mean_row[num]))
place$R <- round((place$tooth_to - place$tooth_from) / (place$vert_to - place$vert_from), 1)
place[num] <- lapply(place[num], round)
```

``` r
cols <- c("A · 06c cut" = "#1b9e77", "B · 06d build" = "#d95f02", "C · 06d free" = "#7570b3",
          "D · 06e mix" = "#e7298a", "Mean" = "black")
zs <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

plot_fish <- function(f){
  rec <- fish[[f]]; bd <- data.frame(band = rec$bands); el <- cmp_settings$main
  pa  <- paths[paths$fish == f, ]
  # 1: the tooth laid on the vertebra by each method (main element, z-scores)
  vert <- data.frame(vert = rec$vert$true_dist, value = zs(rec$vert[[el]]))
  tz <- data.frame(tooth = rec$tooth$true_dist, value = zs(rec$tooth[[el]]))
  lay <- do.call(rbind, lapply(split(pa, pa$label), function(p){
    t <- tz[tz$tooth >= min(p$tooth) & tz$tooth <= max(p$tooth), ]
    data.frame(method = p$method[1], label = p$label[1], value = t$value,
               vert = approx(p$tooth, p$vert, xout = t$tooth, ties = mean)$y)
  }))
  vert_m <- merge(vert, data.frame(method = unique(pa$method)))
  p1 <- ggplot() +
    geom_vline(data = bd, aes(xintercept = band), colour = "brown", linewidth = 0.3) +
    geom_line(data = vert_m, aes(vert, value), colour = "grey55", linewidth = 0.3) +
    geom_line(data = lay, aes(vert, value, colour = method, group = label), linewidth = 0.35) +
    scale_colour_manual(values = cols, guide = "none") +
    facet_wrap(~ method, nrow = 1) +
    theme_bw() + labs(x = "Vertebra (um)", y = paste(sub("_Ca.*$", "", el), "z"))
  # 2: start and end of each placement, with the tooth it covers
  pl <- place[place$fish == f, ]
  pl$label <- factor(pl$label, levels = rev(pl$label))
  p2 <- ggplot(pl, aes(y = label)) +
    geom_vline(xintercept = rec$bands, colour = "brown", linewidth = 0.3) +
    geom_segment(aes(x = vert_from, xend = vert_to, yend = label, colour = method), linewidth = 2.5) +
    geom_text(aes(x = vert_to, label = sprintf("  tooth %.1f-%.1f mm, R %.1f", tooth_from / 1000, tooth_to / 1000, R)),
              hjust = 0, size = 2.8) +
    scale_colour_manual(values = cols, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.35))) +
    theme_bw() + labs(x = "Vertebra (um)", y = NULL)
  # 3: each element's path for the uncut tooth (06c), and the mean of those kept
  ep <- el_paths[el_paths$fish == f, ]
  mp <- mean_path(f)
  p3 <- ggplot() +
    geom_vline(data = bd, aes(xintercept = band), colour = "brown", linewidth = 0.3) +
    geom_path(data = ep, aes(vert, tooth, colour = element, linetype = kept), linewidth = 0.6) +
    scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "22"), name = "Kept (A)") +
    theme_bw() + labs(x = "Vertebra (um)", y = "Tooth (um)", colour = "Element")
  if(!is.null(mp)) p3 <- p3 + geom_path(data = mp, aes(vert, tooth), colour = "black", linewidth = 1.2)
  (p1 / p2 / p3) + plot_layout(heights = c(1, 0.9, 1.4)) +
    plot_annotation(caption = paste(c(sprintf("%s.", f), strwrap(paste(
      "Top: the tooth laid on the vertebra by each method (A-D; ", sub("_Ca.*$", "", el), " as z-scores; grey =",
      "vertebra, colour = tooth; one line per version in B and C). Middle: where each puts the tooth on the vertebra,",
      "with the stretch of tooth it covers and the growth ratio R (um of tooth per um of vertebra); Mean = A-D",
      "averaged, versions averaged within B and C first. Bottom: each element's warp path for the whole, uncut tooth",
      "(06c); solid = kept for A (its placement moves the right way when the tooth is cut), dashed = not kept; black =",
      "mean of the kept paths (A). Brown lines = vertebral bands."), 115)), collapse = "\n"),
      theme = theme(plot.caption = element_text(hjust = 0)))
}
figs <- lapply(setNames(names(fish), names(fish)), plot_fish)
```

### SMA001W

``` r
figs$SMA001W
```

![](06f_dtw_plots_files/figure-gfm/sma001w-1.png)<!-- -->

### SMA001B

``` r
figs$SMA001B
```

![](06f_dtw_plots_files/figure-gfm/sma001b-1.png)<!-- -->

### SMA002W

``` r
figs$SMA002W
```

![](06f_dtw_plots_files/figure-gfm/sma002w-1.png)<!-- -->

## Tables

``` r
for(f in names(fish)){
  cat("\n\n###", f, "\n\n")
  print(knitr::kable(place[place$fish == f, c("label", num, "R")], row.names = FALSE))
  ct <- cut_tab[cut_tab$fish == f, c("element", "rho_start", "rho_end", "kept", "vert_from", "vert_to")]
  ct$element <- sub("_Ca.*$", "", ct$element); ct[c("rho_start", "rho_end")] <- round(ct[c("rho_start", "rho_end")], 2)
  ct[c("vert_from", "vert_to")] <- round(ct[c("vert_from", "vert_to")])
  cat("\n\n")
  print(knitr::kable(ct, row.names = FALSE))
}
```

### SMA001W

| label | vert_from | vert_to | tooth_from | tooth_to | R |
|:---|---:|---:|---:|---:|---:|
| A · 06c cut (Mg+Zn) | 2702 | 10457 | 864 | 54832 | 7.0 |
| B · 06d build SMA001W1 (Mg+Sr) | 5911 | 9061 | 40300 | 49700 | 3.0 |
| B · 06d build SMA001W3 (Zn+Sr) | 9711 | 10461 | 46900 | 51150 | 5.7 |
| C · 06d free SMA001W1 (Mg+Sr+Ba+Mn+Zn) | 5911 | 9061 | 36964 | 42031 | 1.6 |
| C · 06d free SMA001W3 (Mg+Sr+Ba+Mn+Zn) | 9711 | 10461 | 30464 | 32355 | 2.5 |
| D · 06e mix (Mg+Ba+Mn+Zn) | 3687 | 9064 | 864 | 54832 | 10.0 |
| Mean of A-D | 5503 | 9761 | 19761 | 49321 | 6.9 |

| element | rho_start | rho_end | kept  | vert_from | vert_to |
|:--------|----------:|--------:|:------|----------:|--------:|
| Ba      |      0.37 |      NA | FALSE |      7162 |   10647 |
| Li      |      0.33 |   -0.64 | FALSE |      3529 |    4843 |
| Mg      |      0.64 |    0.99 | TRUE  |      2185 |   10258 |
| Mn      |     -0.10 |   -0.85 | FALSE |      3827 |    8825 |
| Sr      |     -0.35 |    0.42 | FALSE |      6913 |   10657 |
| Zn      |      0.27 |    0.26 | TRUE  |      3220 |   10657 |

### SMA001B

| label | vert_from | vert_to | tooth_from | tooth_to | R |
|:---|---:|---:|---:|---:|---:|
| A · 06c cut (Mg+Mn+Zn) | 4499 | 9802 | 12002 | 52836 | 7.7 |
| B · 06d build SMA001B1 (Mg+Sr) | 5546 | 8296 | 34602 | 48400 | 5.0 |
| B · 06d build SMA001B3 (Zn+Sr) | 5746 | 9896 | 34602 | 46250 | 2.8 |
| C · 06d free SMA001B1 (Mg+Sr+Ba+Mn+Zn) | 5546 | 8296 | 34002 | 40099 | 2.2 |
| C · 06d free SMA001B3 (Mg+Sr+Ba+Mn+Zn) | 5746 | 9896 | 34302 | 39715 | 1.3 |
| D · 06e mix (Mg+Sr+Ba+Li+Mn) | 4096 | 9928 | 12002 | 52836 | 7.0 |
| Mean of A-D | 4972 | 9481 | 23189 | 48226 | 5.6 |

| element | rho_start | rho_end | kept  | vert_from | vert_to |
|:--------|----------:|--------:|:------|----------:|--------:|
| Ba      |     -0.64 |    0.05 | FALSE |      7834 |    9956 |
| Li      |      0.21 |   -0.58 | FALSE |      2702 |    4437 |
| Mg      |      0.29 |    0.94 | TRUE  |      3865 |    9799 |
| Mn      |      0.06 |    0.63 | TRUE  |      4400 |    9965 |
| Sr      |      0.72 |   -0.26 | FALSE |      1337 |    6495 |
| Zn      |      0.05 |    0.29 | TRUE  |      5231 |    9642 |

### SMA002W

| label | vert_from | vert_to | tooth_from | tooth_to | R |
|:---|---:|---:|---:|---:|---:|
| A · 06c cut (Zn) | 814 | 7286 | 156 | 42969 | 6.6 |
| B · 06d build SMA002W1 (Mg+Sr) | 1534 | 7884 | 4956 | 37650 | 5.1 |
| B · 06d build SMA002W2 (Mg+Ba) | 4134 | 6284 | 26300 | 41150 | 6.9 |
| B · 06d build SMA002W3 (Zn+Sr) | 334 | 2084 | 1556 | 16900 | 8.8 |
| C · 06d free SMA002W1 (Mg+Sr+Ba+Mn+Zn) | 1534 | 7884 | 23256 | 33471 | 1.6 |
| C · 06d free SMA002W2 (Mg+Sr+Ba+Mn+Zn) | 4134 | 6284 | 27756 | 33177 | 2.5 |
| C · 06d free SMA002W3 (Mg+Sr+Ba+Mn+Zn) | 334 | 2084 | 32356 | 37834 | 3.1 |
| D · 06e mix (Mg+Sr+Ba+Mn+Zn) | 3931 | 5769 | 156 | 42969 | 23.3 |
| Mean of A-D | 2187 | 5972 | 9759 | 38166 | 7.5 |

| element | rho_start | rho_end | kept  | vert_from | vert_to |
|:--------|----------:|--------:|:------|----------:|--------:|
| Ba      |      0.83 |   -0.95 | FALSE |      1912 |    5777 |
| Li      |      0.09 |      NA | FALSE |      5548 |    7949 |
| Mg      |     -0.40 |   -0.01 | FALSE |      3635 |    5899 |
| Mn      |      0.23 |   -0.04 | FALSE |       357 |    1882 |
| Sr      |     -0.55 |    0.23 | FALSE |      7126 |    7949 |
| Zn      |      0.89 |    0.18 | TRUE  |       814 |    7286 |

- **First table**: each placement: where on the vertebra (**vert_from,
  vert_to**) and which stretch of tooth (**tooth_from, tooth_to**), µm;
  **R** = µm of tooth per µm of vertebra.
- **Second table**, per element (06c): **rho_start** = how consistently
  the vertebra start moves later as more is cut from the tooth’s start
  (Spearman ρ over all cuts; 1 = always); **rho_end** = the same for the
  end moving earlier (NA = it never moves); **kept** = both above
  `rho_min` (used for A); **vert_from, vert_to** = where it puts the
  whole, uncut tooth.

``` r
dir.create('data/output/dtw_final', showWarnings = FALSE)
dir.create('figures/dtw', showWarnings = FALSE)
write.csv(place, 'data/output/dtw_final/placements.csv', row.names = FALSE)
write.csv(cut_tab, 'data/output/dtw_final/cut_response.csv', row.names = FALSE)
ggsave('figures/dtw/buildout_versions.png', p_ver, width = 10, height = 10, dpi = 300)
for(f in names(figs)) ggsave(paste0('figures/dtw/placements_', f, '.png'), figs[[f]], width = 10, height = 12, dpi = 300)
save(place, cut_tab, el_paths, paths, built, best_mix, cmp_settings, file = 'data/output/dtw_cache/06f.RData')
```

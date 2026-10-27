06g · DTW summary statistics
================
Nikki Biskis
Rendered 27 Oct 2026 22:02

- [1. Where the tooth sits](#1-where-the-tooth-sits)
  - [SMA001W](#sma001w)
  - [SMA001B](#sma001b)
  - [SMA002W](#sma002w)
- [2. Elements](#2-elements)
  - [SMA001W](#sma001w-1)
  - [SMA001B](#sma001b-1)
  - [SMA002W](#sma002w-1)
- [3. Growth ratio along the tooth](#3-growth-ratio-along-the-tooth)

- **Paper:** The potential of elemental analysis in determining origin
  of sawfish rostra
- **DOI:** tbd
- **Written:** 30 Aug 2026 (summary tables Oct 2026)
- **Content:** dynamic time warping VII: the numbers to report

Starts from 06c, 06d, 06e and 06f (`data/output/dtw_cache/`). The
previous version (warping and misfit along the vertebra) is in
`original/06g_dtw_stats_old.Rmd`.

Summary statistics, as other studies of paired records report them: how
consistently methods and elements place the tooth (agreement, mean ± SD;
e.g. Clark et al. 2020 compared weaning ages between Sr and Ba and
between methods), which elements carry the alignment, the growth ratio
and how it changes along the tooth, and correlation of the aligned
profiles. The per-element measures follow `original/06e_dtw_rates.R`.

1.  **Where the tooth sits**: the vertebra position each method (06f,
    A–D) gives points along the tooth, every 5 mm, and their mean ± SD.
2.  **Elements**: for each element’s alignment of the whole, uncut tooth
    (06c), where it starts and ends, its growth ratio, how far it
    departs from constant growth, how far it sits from the other
    elements, how well the aligned profiles correlate, and whether it
    was used by A, B or D.
3.  **Growth ratio along the tooth**: from 06d’s build windows, the
    growth ratio and whether it changes along the tooth.

``` r
source('helpers/microHelp.R')
load('data/output/dtw_cache/06c.RData')  # fish (records, bands) from 06c
load('data/output/dtw_cache/06d.RData')  # anchor_pairs, up from 06d
load('data/output/dtw_cache/06f.RData')  # place, cut_tab, el_paths, paths, built, best_mix from 06f
ref_um <- 5000                           # table 1: a point every 5 mm of tooth
fish_of <- function(v) sub("[0-9]+$", "", v)
short <- function(e) sub("_Ca.*$", "", e)
```

## 1. Where the tooth sits

``` r
where <- do.call(rbind, lapply(names(fish), function(f){
  pa <- paths[paths$fish == f, ]
  tr <- seq(ref_um, max(pa$tooth), by = ref_um)
  v <- do.call(rbind, lapply(split(pa, pa$label), function(p) data.frame(method = substr(p$method[1], 1, 1),
    tooth = tr, vert = approx(p$tooth, p$vert, xout = tr, ties = mean)$y)))   # NA outside what it covers
  v <- aggregate(vert ~ method + tooth, v, mean)                                # versions averaged within B, C
  w <- reshape(v, idvar = "tooth", timevar = "method", direction = "wide")
  names(w) <- sub("^vert\\.", "", names(w))
  m <- as.matrix(w[intersect(c("A", "B", "C", "D"), names(w))])
  data.frame(fish = f, tooth_mm = w$tooth / 1000, round(m), mean = round(rowMeans(m, na.rm = TRUE)),
             sd = round(apply(m, 1, sd, na.rm = TRUE)), n = rowSums(!is.na(m)))
}))
where <- where[order(where$fish, where$tooth_mm), ]
for(f in names(fish)){
  cat("\n\n###", f, "\n\n")
  print(knitr::kable(where[where$fish == f, -1], row.names = FALSE))
}
```

### SMA001W

| tooth_mm |     A |     B |    C |    D | mean |   sd |   n |
|---------:|------:|------:|-----:|-----:|-----:|-----:|----:|
|        5 |  3848 |    NA |   NA | 3916 | 3882 |   48 |   2 |
|       10 |  5019 |    NA |   NA | 4623 | 4821 |  280 |   2 |
|       15 |  5535 |    NA |   NA | 5347 | 5441 |  133 |   2 |
|       20 |  5779 |    NA |   NA | 5789 | 5784 |    7 |   2 |
|       25 |  6244 |    NA |   NA | 5789 | 6017 |  322 |   2 |
|       30 |  7006 |    NA |   NA | 5948 | 6477 |  749 |   2 |
|       35 |  7452 |    NA |   NA | 5959 | 6706 | 1056 |   2 |
|       40 |  8112 |    NA | 7799 | 6890 | 7600 |  635 |   3 |
|       45 |  9372 |  7986 |   NA | 8603 | 8654 |  694 |   3 |
|       50 | 10093 | 10258 |   NA | 9064 | 9805 |  647 |   3 |

### SMA001B

| tooth_mm |    A |    B |    C |    D | mean |   sd |   n |
|---------:|-----:|-----:|-----:|-----:|-----:|-----:|----:|
|       15 | 4711 |   NA |   NA | 4124 | 4418 |  415 |   2 |
|       20 | 5142 |   NA |   NA | 4124 | 4633 |  719 |   2 |
|       25 | 5670 |   NA |   NA | 4133 | 4901 | 1086 |   2 |
|       30 | 6366 |   NA |   NA | 5804 | 6085 |  398 |   2 |
|       35 | 6729 | 5882 | 6139 | 6138 | 6222 |  359 |   4 |
|       40 | 8252 | 7751 | 8251 | 8036 | 8073 |  237 |   4 |
|       45 | 9103 | 8472 |   NA | 8701 | 8759 |  319 |   3 |
|       50 | 9787 |   NA |   NA | 9921 | 9854 |   95 |   2 |

### SMA002W

| tooth_mm |    A |    B |    C |    D | mean |   sd |   n |
|---------:|-----:|-----:|-----:|-----:|-----:|-----:|----:|
|        5 | 1020 | 1330 |   NA | 3931 | 2094 | 1599 |   3 |
|       10 | 1618 | 1832 |   NA | 3963 | 2471 | 1297 |   3 |
|       15 | 2789 | 2376 |   NA | 4214 | 3126 |  964 |   3 |
|       20 | 3345 | 4175 |   NA | 4251 | 3924 |  502 |   3 |
|       25 | 4126 | 4294 | 2618 | 4764 | 3951 |  928 |   4 |
|       30 | 5310 | 4884 | 5375 | 4862 | 5108 |  273 |   4 |
|       35 | 6951 | 5449 | 1179 | 4870 | 4612 | 2451 |   4 |
|       40 | 7065 | 6081 |   NA | 4870 | 6005 | 1099 |   3 |

- **tooth_mm**: position along the tooth (mm); **A–D**: the vertebra
  position (µm) each method puts it at (blank = the method doesn’t cover
  that part of the tooth; B and C averaged over versions); **mean, sd,
  n**: over the methods that cover it.

``` r
ends <- do.call(rbind, lapply(split(place[place$method != "Mean", ], place$fish[place$method != "Mean"]), function(p){
  w <- p[p$method %in% c("A · 06c cut", "D · 06e mix"), ]   # the two that place the whole tooth
  data.frame(fish = p$fish[1], vert_start = sprintf("%.0f ± %.0f", mean(w$vert_from), sd(w$vert_from)),
             vert_end = sprintf("%.0f ± %.0f", mean(w$vert_to), sd(w$vert_to)),
             R = sprintf("%.1f ± %.1f", mean(p$R), sd(p$R)), n_R = nrow(p))
}))
knitr::kable(ends, row.names = FALSE)
```

| fish    | vert_start  | vert_end    | R         | n_R |
|:--------|:------------|:------------|:----------|----:|
| SMA001B | 4298 ± 285  | 9865 ± 89   | 4.3 ± 2.6 |   6 |
| SMA001W | 3194 ± 697  | 9760 ± 985  | 5.0 ± 3.2 |   6 |
| SMA002W | 2372 ± 2204 | 6528 ± 1073 | 7.2 ± 6.9 |   8 |

- **vert_start, vert_end**: where the whole tooth starts and ends on the
  vertebra (µm, mean ± SD of A and D, the two methods that place the
  whole tooth).
- **R**: µm of tooth per µm of vertebra, mean ± SD over every placement
  in 06f (A–D, each version); **n_R**: how many.

## 2. Elements

``` r
el_tab <- do.call(rbind, lapply(names(fish), function(f){
  ep <- el_paths[el_paths$fish == f, ]
  tg <- seq(max(tapply(ep$tooth, ep$element, min)), min(tapply(ep$tooth, ep$element, max)), length.out = 500)
  V  <- sapply(split(ep, ep$element), function(e) approx(e$tooth, e$vert, xout = tg, ties = mean)$y)
  rec <- fish[[f]]
  do.call(rbind, lapply(split(ep, ep$element), function(e){
    el <- e$element[1]; col <- paste0(el, "_Ca_mav"); ct <- cut_tab[cut_tab$fish == f & short(cut_tab$element) == el, ]
    vhat <- e$vert[1] + (e$vert[nrow(e)] - e$vert[1]) * (e$tooth - e$tooth[1]) / (e$tooth[nrow(e)] - e$tooth[1])
    tz <- approx(rec$tooth$true_dist, rec$tooth[[col]], xout = e$tooth, ties = mean)$y   # tooth profile ...
    vz <- approx(rec$vert$true_dist, rec$vert[[col]], xout = e$vert, ties = mean)$y      # ... and the vertebra where it's matched
    data.frame(fish = f, element = el,
               vert_start = round(min(e$vert)), vert_end = round(max(e$vert)),
               R = round(diff(range(e$tooth)) / diff(range(e$vert)), 1),
               delta_um = round(mean(abs(e$vert - vhat))),
               Dbar_um = round(mean(colMeans(abs(V[, setdiff(colnames(V), el), drop = FALSE] - V[, el]), na.rm = TRUE))),
               r = round(cor(tz, vz, use = "complete.obs"), 2),
               rho_start = round(ct$rho_start, 2), rho_end = round(ct$rho_end, 2),
               A = ct$kept,
               B = sum(sapply(names(up)[fish_of(names(up)) == f], function(n) el %in% short(anchor_pairs[[n]]))),
               D = el %in% strsplit(best_mix$elements[best_mix$set == f], "\\+")[[1]])
  }))
}))
for(f in names(fish)){
  cat("\n\n###", f, "\n\n")
  print(knitr::kable(el_tab[el_tab$fish == f, -1], row.names = FALSE))
}
```

### SMA001W

| element | vert_start | vert_end | R | delta_um | Dbar_um | r | rho_start | rho_end | A | B | D |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|:---|---:|:---|
| Ba | 7162 | 10647 | 15.5 | 806 | 2920 | 0.59 | 0.37 | NA | FALSE | 0 | TRUE |
| Li | 3529 | 4843 | 41.1 | 233 | 4034 | 0.70 | 0.33 | -0.64 | FALSE | 0 | FALSE |
| Mg | 2185 | 10258 | 6.7 | 712 | 2300 | 0.71 | 0.64 | 0.99 | TRUE | 1 | TRUE |
| Mn | 3827 | 8825 | 10.8 | 514 | 2164 | 0.60 | -0.10 | -0.85 | FALSE | 0 | TRUE |
| Sr | 6913 | 10657 | 14.4 | 1256 | 3205 | 0.77 | -0.35 | 0.42 | FALSE | 2 | FALSE |
| Zn | 3220 | 10657 | 7.3 | 575 | 2108 | 0.67 | 0.27 | 0.26 | TRUE | 1 | TRUE |

### SMA001B

| element | vert_start | vert_end | R | delta_um | Dbar_um | r | rho_start | rho_end | A | B | D |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|:---|---:|:---|
| Ba | 7834 | 9956 | 19.2 | 414 | 3401 | 0.50 | -0.64 | 0.05 | FALSE | 0 | TRUE |
| Li | 2702 | 4437 | 23.5 | 345 | 3796 | 0.86 | 0.21 | -0.58 | FALSE | 0 | TRUE |
| Mg | 3865 | 9799 | 6.9 | 693 | 1952 | 0.82 | 0.29 | 0.94 | TRUE | 1 | TRUE |
| Mn | 4400 | 9965 | 7.3 | 365 | 1911 | 0.82 | 0.06 | 0.63 | TRUE | 0 | TRUE |
| Sr | 1337 | 6495 | 7.9 | 1232 | 2393 | 0.78 | 0.72 | -0.26 | FALSE | 2 | TRUE |
| Zn | 5231 | 9642 | 9.3 | 341 | 2033 | 0.83 | 0.05 | 0.29 | TRUE | 1 | FALSE |

### SMA002W

| element | vert_start | vert_end | R | delta_um | Dbar_um | r | rho_start | rho_end | A | B | D |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|:---|---:|:---|
| Ba | 1912 | 5777 | 11.1 | 660 | 2239 | 0.87 | 0.83 | -0.95 | FALSE | 1 | TRUE |
| Li | 5548 | 7949 | 17.8 | 278 | 2961 | 0.82 | 0.09 | NA | FALSE | 0 | FALSE |
| Mg | 3635 | 5899 | 18.9 | 365 | 2265 | 0.83 | -0.40 | -0.01 | FALSE | 2 | TRUE |
| Mn | 357 | 1882 | 28.1 | 211 | 4593 | 0.70 | 0.23 | -0.04 | FALSE | 0 | TRUE |
| Sr | 7126 | 7949 | 52.0 | 283 | 3697 | 0.59 | -0.55 | 0.23 | FALSE | 2 | TRUE |
| Zn | 814 | 7286 | 6.6 | 430 | 2512 | 0.68 | 0.89 | 0.18 | TRUE | 1 | TRUE |

- **vert_start, vert_end**: where the element puts the whole tooth on
  the vertebra (µm); **R**: µm of tooth per µm of vertebra over that.
- **delta_um** (δ): mean distance (µm of vertebra) between the element’s
  path and a straight line from its start to its end: how far growth
  departs from a constant ratio.
- **Dbar_um** (D̄): mean distance (µm of vertebra) between the element’s
  path and each other element’s, over the tooth they share: how far it
  sits from the rest.
- **r**: Pearson correlation between the tooth profile and the vertebra
  profile at the matched positions (06a’s rolling medians). DTW makes
  the profiles as alike as it can, so r is high by construction; compare
  it between elements, not as a test.
- **rho_start, rho_end**: how consistently the element’s vertebra start
  moves later as the tooth’s start is cut, and its end earlier as the
  tooth’s end is cut (Spearman ρ over 06c’s cuts; NA = never moves).
- **A**: kept for 06f’s A (both ρ above `rho_min`); **B**: number of
  anchor versions using it that built (06d); **D**: in 06f’s D mix
  (06e).

## 3. Growth ratio along the tooth

``` r
R_tab <- do.call(rbind, lapply(split(built, built$set), function(b){
  if(nrow(b) < 3) return(data.frame(version = b$version[1], windows = nrow(b), tooth_mm = NA, R_median = round(median(b$ratio), 1),
                                    R_IQR = NA, rho = NA, p = NA, R_first_half = NA, R_second_half = NA))
  ct <- suppressWarnings(cor.test(b$t_mid, b$ratio, method = "spearman"))
  half <- b$t_mid < median(b$t_mid)
  data.frame(version = b$version[1], windows = nrow(b),
             tooth_mm = sprintf("%.1f-%.1f", min(b$tooth_start) / 1000, max(b$tooth_end) / 1000),
             R_median = round(median(b$ratio), 1),
             R_IQR = sprintf("%.1f-%.1f", quantile(b$ratio, 0.25), quantile(b$ratio, 0.75)),
             rho = round(ct$estimate, 2), p = signif(ct$p.value, 2),
             R_first_half = round(median(b$ratio[half]), 1), R_second_half = round(median(b$ratio[!half]), 1))
}))
knitr::kable(R_tab, row.names = FALSE)
```

| version | windows | tooth_mm | R_median | R_IQR | rho | p | R_first_half | R_second_half |
|:---|---:|:---|---:|:---|---:|---:|---:|---:|
| SMA001B1 (Mg+Sr) | 11 | 34.6-48.5 | 1.9 | 1.5-4.5 | 0.64 | 0.034 | 1.5 | 3.1 |
| SMA001B3 (Zn+Sr) | 18 | 34.6-46.2 | 1.9 | 1.0-2.9 | -0.29 | 0.250 | 2.9 | 1.5 |
| SMA001W1 (Mg+Sr) | 13 | 40.3-50.2 | 2.4 | 1.5-5.2 | 0.32 | 0.290 | 1.9 | 5.2 |
| SMA001W3 (Zn+Sr) | 1 | NA | 5.7 | NA | NA | NA | NA | NA |
| SMA002W1 (Mg+Sr) | 29 | 5.0-39.0 | 2.9 | 1.0-5.7 | 0.16 | 0.410 | 1.5 | 2.9 |
| SMA002W2 (Mg+Ba) | 8 | 26.3-41.1 | 4.0 | 2.2-6.1 | 0.27 | 0.520 | 1.9 | 5.7 |
| SMA002W3 (Zn+Sr) | 6 | 1.6-16.9 | 1.7 | 1.5-3.7 | -0.20 | 0.700 | 4.3 | 1.5 |

- **version**: 06d’s anchor version (anchor pair); **windows**: build
  windows; **tooth_mm**: the stretch of tooth they matched.
- **R_median, R_IQR**: each window’s growth ratio (µm of tooth per µm of
  vertebra), median and interquartile range.
- **rho, p**: Spearman correlation of the ratio with position along the
  tooth (positive = more tooth per µm of vertebra further along the
  tooth); a correlation needs 3+ windows.
- **R_first_half, R_second_half**: median ratio over the first and
  second half of the matched tooth.

``` r
dir.create('data/output/dtw_final', showWarnings = FALSE)
write.csv(where, 'data/output/dtw_final/tooth_on_vert.csv', row.names = FALSE)
write.csv(el_tab, 'data/output/dtw_final/elements.csv', row.names = FALSE)
write.csv(R_tab, 'data/output/dtw_final/growth_ratio.csv', row.names = FALSE)
save(where, ends, el_tab, R_tab, file = 'data/output/dtw_cache/06g.RData')
```

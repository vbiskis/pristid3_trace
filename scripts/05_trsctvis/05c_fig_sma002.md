05c · SMA002 profiles
================
Nikki Biskis
Rendered 26 Oct 2026 17:09

- [SMA002](#sma002)
  - [Bleach](#bleach)
  - [Water](#water)
  - [Vert](#vert)
  - [Figure (supp S4)](#figure-supp-s4)

- **Paper:** The potential of elemental analysis in determining origin
  of sawfish rostra
- **DOI:** tbd
- **Written:** 16 Jul 2024
- **Content:** vert and teeth pics for aquarium samples. This one’s
  manual, doesn’t use the loop, can adjust.

R markdown version of `05c_fig_sma002.R`. Rendering it also saves
`figures/supp/s4/SMA002VT.png`.

``` r
source('helpers/microHelp.R')
source('helpers/functions/fxn_mavs.R')
source('helpers/functions/fxn_profile.R')
profclean <- readRDS('data/output/prof_processed.RDS')
transsplit <- readRDS('data/output/transsplit.RDS')
names(transsplit)[54:60]
```

    #> [1] "SMA001BT-TR" "SMA001HT"    "SMA001BV"    "SMA001HV"    "SMA002BT-TR"
    #> [6] "SMA002HT"    "SMA002BV"

``` r
# want 58-60 for SMA002
#58 tooth bleach
#59 tooth water
#60 vert bleach

#issues with vert water (in general not as high performing)
```

## SMA002

### Bleach

``` r
SMA002TB <- transsplit[[58]]

# lets remove In spikes
b2 <- plot_holes(SMA002TB, trigger =  2.5, release = 1,
                 elements = c("Sr_Ca", "Ba_Ca", "Li_Ca", "Mg_Ca")) # lol welp
```

![](05c_fig_sma002_files/figure-gfm/sma002tb-holes-1.png)<!-- -->

``` r
SMA002TB <- SMA002TB %>%
  filter(# true_dist >= min(true_dist[which(In_Ca < 5)]), # this tooth is rough
         true_dist <= 49907)
```

``` r
SMA002TB_out <- profs(SMA002TB, el_all, "tooth", span = 0.01, window_size = 21, # stats
                      show_In = TRUE, In_style = "rect",
                      rect_threshold = .5)
```

![](05c_fig_sma002_files/figure-gfm/sma002tb-1.png)<!-- -->![](05c_fig_sma002_files/figure-gfm/sma002tb-2.png)<!-- -->![](05c_fig_sma002_files/figure-gfm/sma002tb-3.png)<!-- -->

``` r
SMA002TB_out$plots$Ca
```

![](05c_fig_sma002_files/figure-gfm/sma002tb-4.png)<!-- -->

``` r
sma02a2 <- SMA002TB_out$plots$Sr + theme(axis.title.x = element_blank())
sma02b2 <- SMA002TB_out$plots$Ba + theme(axis.title.x = element_blank())
sma02c2 <- SMA002TB_out$plots$Li + theme(axis.title.x = element_blank())
sma02d2 <- SMA002TB_out$plots$Mg + theme(axis.title.x = element_blank())
sma02e2 <- SMA002TB_out$plots$Mn
```

### Water

``` r
SMA002TH <- transsplit[[59]]

h2 <- plot_holes(SMA002TH, trigger =  2.5, release = 1.25,
                 elements = c("Sr_Ca", "Ba_Ca", "Li_Ca", "Mg_Ca"))
```

![](05c_fig_sma002_files/figure-gfm/sma002th-holes-1.png)<!-- -->

``` r
SMA002TH <- h2$data[!h2$data$.hole, ] # lol! k
```

#### Holes and smoothing

``` r
find_gaps(SMA002TH$true_dist) # holes > 50 um (normal spacing ~8.6 um)
```

    #>       start       end
    #> 1  1477.518  1598.484
    #> 2  2142.831  2272.437
    #> 3 18931.227 19484.214
    #> 4 20192.730 20287.776

``` r
plot_smoothing(SMA002TH, c("Sr_Ca", "Ba_Ca", "Li_Ca", "Mg_Ca"), span = 0.01)
```

![](05c_fig_sma002_files/figure-gfm/sma002th-gaps-1.png)<!-- -->

``` r
plot_smoothing(SMA002TH, c("Sr_Ca", "Ba_Ca", "Li_Ca", "Mg_Ca"), span = 0.05, xlim = c(13000, 25000))
```

![](05c_fig_sma002_files/figure-gfm/water-gaps-zoom-1.png)<!-- -->

`trim_gaps()` drops points within `buffer` um of a hole:

``` r
plot_smoothing(trim_gaps(SMA002TH, buffer = 400), c("Sr_Ca", "Ba_Ca", "Li_Ca", "Mg_Ca"),
               span = 0.05, xlim = c(13000, 25000))
```

![](05c_fig_sma002_files/figure-gfm/water-gaps-trim-1.png)<!-- -->

``` r
SMA002TH_out <- profs(SMA002TH, el_all, "tooth", 
                      span = 0.01, window_size = 21, gap_um = 200,
                      # img_path = "data/transects/imgs/SMA002/SMA002-Trsct.jpg", # still big
                      show_In = TRUE, In_style = "rect")
```

![](05c_fig_sma002_files/figure-gfm/sma002th-1.png)<!-- -->![](05c_fig_sma002_files/figure-gfm/sma002th-2.png)<!-- -->![](05c_fig_sma002_files/figure-gfm/sma002th-3.png)<!-- -->

``` r
SMA002TH_out$plots$Ca
```

![](05c_fig_sma002_files/figure-gfm/sma002th-4.png)<!-- -->

``` r
sma02a <- SMA002TH_out$plots$Sr + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
                                        legend.position = "none")
sma02b <- SMA002TH_out$plots$Ba + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
                                        legend.position = "none")
sma02c <- SMA002TH_out$plots$Li + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
                                        legend.position = "none")
sma02d <- SMA002TH_out$plots$Mg + theme(axis.title.y = element_blank(), 
                                        legend.position = "none")
sma02e <- SMA002TH_out$plots$Mn 
```

### Vert

``` r
SMA002V <- transsplit[[60]]
SMA002V <- SMA002V %>%
  filter(true_dist >= true_dist[which(In_Ca < 10)][1], # adjust to see window move (2, 5, 10)
         true_dist <= 7950) # bc it comes back down, manual here

SMA002V_out <- profs(SMA002V, el_all, "vert", span = 0.025, show_In = FALSE)
```

![](05c_fig_sma002_files/figure-gfm/sma002v-1.png)<!-- -->![](05c_fig_sma002_files/figure-gfm/sma002v-2.png)<!-- -->![](05c_fig_sma002_files/figure-gfm/sma002v-3.png)<!-- -->

``` r
SMA002V_out$plots$Ca
```

![](05c_fig_sma002_files/figure-gfm/sma002v-4.png)<!-- -->

``` r
sma02va <- SMA002V_out$plots$Sr + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
                                        legend.position = "none")
sma02vb <- SMA002V_out$plots$Ba + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
                                        legend.position = "none")
sma02vc <- SMA002V_out$plots$Li + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
                                        legend.position = "none")
sma02vd <- SMA002V_out$plots$Mg + theme(axis.title.y = element_blank(), legend.position = "none")
sma02ve <- SMA002V_out$plots$Mn
```

### Figure (supp S4)

``` r
sma002 <- (sma02a + sma02va) / (sma02b + sma02vb) / (sma02c + sma02vc) /
  (sma02d + sma02vd) & theme(plot.margin = margin(5, 5, 5, 5))

ggsave(
  "SMA002VT.png",
  plot = sma002,
  device = NULL,
  path = 'figures/supp/s4',
  scale = 1,
  width = 8,
  height = 8,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
)

sma002
```

![](05c_fig_sma002_files/figure-gfm/sma002-supp-1.png)<!-- -->

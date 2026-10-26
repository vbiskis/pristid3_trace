06b · DTW run
================
Nikki Biskis
Rendered 24 Oct 2026 19:14

- [SMA001](#sma001)
  - [Multi: bleach](#multi-bleach)
  - [Multi: water](#multi-water)
- [SMA002](#sma002)
  - [Multi](#multi)
- [Save for 06f](#save-for-06f)

- **Paper:** The potential of elemental analysis in determining origin
  of sawfish rostra

- **DOI:** tbd

- **Written:** 28 Aug 2026

- **Content:** dynamic time warping II

- application of functions developed previously

- setting the models built with parameters developed in dtw1

R markdown version of `06b_dtw_run.R`. It starts from the objects 06a
saved (`data/output/dtw_cache/06a.RData`), and saves its own objects to
`data/output/dtw_cache/06b.RData` for 06f.

``` r
source('helpers/functions/fxn_dtw.R')
source('helpers/functions/fxn_profile.R')
source('helpers/functions/fxn_mavs.R')
source('helpers/microHelp.R')
load('data/output/dtw_cache/06a.RData') # objects from 06a
```

## SMA001

``` r
align_els(SMA001T_filt %>% filter(true_dist <= 43000),
          SMA001V_filt, # Ba gets confused
          elements = names(span_vec),
          spanT = 0.08, spanV = 0.02, n = 2000)
```

    #> Ba_Ca_mav 
    #>  Tooth: 812.202 to 42994.89 
    #>  Vert:  1523.6 to 3168.166 
    #>  Normalized distance: 0.6196438

![](06b_dtw_run_files/figure-gfm/sma001-els-1.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-2.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-3.png)<!-- -->

    #> Li_Ca_mav 
    #>  Tooth: 812.202 to 42994.89 
    #>  Vert:  3896.843 to 7641.038 
    #>  Normalized distance: 0.5091649

![](06b_dtw_run_files/figure-gfm/sma001-els-4.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-5.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-6.png)<!-- -->

    #> Mg_Ca_mav 
    #>  Tooth: 812.202 to 42994.89 
    #>  Vert:  2088.109 to 6949.802 
    #>  Normalized distance: 0.2626141

![](06b_dtw_run_files/figure-gfm/sma001-els-7.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-8.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-9.png)<!-- -->

    #> Mn_Ca_mav 
    #>  Tooth: 812.202 to 42994.89 
    #>  Vert:  1834.656 to 8790.218 
    #>  Normalized distance: 0.5348143

![](06b_dtw_run_files/figure-gfm/sma001-els-10.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-11.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-12.png)<!-- -->

    #> Sr_Ca_mav 
    #>  Tooth: 812.202 to 42994.89 
    #>  Vert:  6949.802 to 10656.56 
    #>  Normalized distance: 0.2945826

![](06b_dtw_run_files/figure-gfm/sma001-els-13.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-14.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-15.png)<!-- -->

    #> Zn_Ca_mav 
    #>  Tooth: 812.202 to 42994.89 
    #>  Vert:  1880.738 to 10512.55 
    #>  Normalized distance: 0.4870439

![](06b_dtw_run_files/figure-gfm/sma001-els-16.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-17.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-18.png)<!-- -->

    #> $Ba_Ca_mav
    #> DTW alignment object
    #>  Alignment size (query x reference): 2000 x 2000
    #>  Call: dtw(x = tSmooth, y = vSmooth, step.pattern = pattern, keep.internals = TRUE, 
    #>  Call:     open.end = TRUE, open.begin = TRUE)
    #> 
    #> $Li_Ca_mav
    #> DTW alignment object
    #>  Alignment size (query x reference): 2000 x 2000
    #>  Call: dtw(x = tSmooth, y = vSmooth, step.pattern = pattern, keep.internals = TRUE, 
    #>  Call:     open.end = TRUE, open.begin = TRUE)
    #> 
    #> $Mg_Ca_mav
    #> DTW alignment object
    #>  Alignment size (query x reference): 2000 x 2000
    #>  Call: dtw(x = tSmooth, y = vSmooth, step.pattern = pattern, keep.internals = TRUE, 
    #>  Call:     open.end = TRUE, open.begin = TRUE)
    #> 
    #> $Mn_Ca_mav
    #> DTW alignment object
    #>  Alignment size (query x reference): 2000 x 2000
    #>  Call: dtw(x = tSmooth, y = vSmooth, step.pattern = pattern, keep.internals = TRUE, 
    #>  Call:     open.end = TRUE, open.begin = TRUE)
    #> 
    #> $Sr_Ca_mav
    #> DTW alignment object
    #>  Alignment size (query x reference): 2000 x 2000
    #>  Call: dtw(x = tSmooth, y = vSmooth, step.pattern = pattern, keep.internals = TRUE, 
    #>  Call:     open.end = TRUE, open.begin = TRUE)
    #> 
    #> $Zn_Ca_mav
    #> DTW alignment object
    #>  Alignment size (query x reference): 2000 x 2000
    #>  Call: dtw(x = tSmooth, y = vSmooth, step.pattern = pattern, keep.internals = TRUE, 
    #>  Call:     open.end = TRUE, open.begin = TRUE)

``` r
align_els(SMA001TB_filt %>% filter(true_dist <= 46000),
          SMA001VB_filt,
          elements = names(span_vec),
          spanT = 0.08, spanV = 0.02, n = 2000)
```

    #> Ba_Ca_mav 
    #>  Tooth: 9530.43 to 45993.12 
    #>  Vert:  861.165 to 1514.959 
    #>  Normalized distance: 0.6130417

![](06b_dtw_run_files/figure-gfm/sma001-els-19.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-20.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-21.png)<!-- -->

    #> Li_Ca_mav 
    #>  Tooth: 9530.43 to 45993.12 
    #>  Vert:  2733.263 to 4458.473 
    #>  Normalized distance: 0.3263849

![](06b_dtw_run_files/figure-gfm/sma001-els-22.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-23.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-24.png)<!-- -->

    #> Mg_Ca_mav 
    #>  Tooth: 9530.43 to 45993.12 
    #>  Vert:  1157.821 to 9066.713 
    #>  Normalized distance: 0.243817

![](06b_dtw_run_files/figure-gfm/sma001-els-25.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-26.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-27.png)<!-- -->

    #> Mn_Ca_mav 
    #>  Tooth: 9530.43 to 45993.12 
    #>  Vert:  4366.308 to 9965.319 
    #>  Normalized distance: 0.2894851

![](06b_dtw_run_files/figure-gfm/sma001-els-28.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-29.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-30.png)<!-- -->

    #> Sr_Ca_mav 
    #>  Tooth: 9530.43 to 45993.12 
    #>  Vert:  1345.03 to 9467.053 
    #>  Normalized distance: 0.2978097

![](06b_dtw_run_files/figure-gfm/sma001-els-31.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-32.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-33.png)<!-- -->

    #> Zn_Ca_mav 
    #>  Tooth: 9530.43 to 45993.12 
    #>  Vert:  2033.386 to 9786.75 
    #>  Normalized distance: 0.2121849

![](06b_dtw_run_files/figure-gfm/sma001-els-34.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-35.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-els-36.png)<!-- -->

    #> $Ba_Ca_mav
    #> DTW alignment object
    #>  Alignment size (query x reference): 2000 x 2000
    #>  Call: dtw(x = tSmooth, y = vSmooth, step.pattern = pattern, keep.internals = TRUE, 
    #>  Call:     open.end = TRUE, open.begin = TRUE)
    #> 
    #> $Li_Ca_mav
    #> DTW alignment object
    #>  Alignment size (query x reference): 2000 x 2000
    #>  Call: dtw(x = tSmooth, y = vSmooth, step.pattern = pattern, keep.internals = TRUE, 
    #>  Call:     open.end = TRUE, open.begin = TRUE)
    #> 
    #> $Mg_Ca_mav
    #> DTW alignment object
    #>  Alignment size (query x reference): 2000 x 2000
    #>  Call: dtw(x = tSmooth, y = vSmooth, step.pattern = pattern, keep.internals = TRUE, 
    #>  Call:     open.end = TRUE, open.begin = TRUE)
    #> 
    #> $Mn_Ca_mav
    #> DTW alignment object
    #>  Alignment size (query x reference): 2000 x 2000
    #>  Call: dtw(x = tSmooth, y = vSmooth, step.pattern = pattern, keep.internals = TRUE, 
    #>  Call:     open.end = TRUE, open.begin = TRUE)
    #> 
    #> $Sr_Ca_mav
    #> DTW alignment object
    #>  Alignment size (query x reference): 2000 x 2000
    #>  Call: dtw(x = tSmooth, y = vSmooth, step.pattern = pattern, keep.internals = TRUE, 
    #>  Call:     open.end = TRUE, open.begin = TRUE)
    #> 
    #> $Zn_Ca_mav
    #> DTW alignment object
    #>  Alignment size (query x reference): 2000 x 2000
    #>  Call: dtw(x = tSmooth, y = vSmooth, step.pattern = pattern, keep.internals = TRUE, 
    #>  Call:     open.end = TRUE, open.begin = TRUE)

``` r
# similar...
```

### Multi: bleach

``` r
mult1b2 <- align_multi(SMA001TB_filt %>% filter(true_dist <= 46000),
                       SMA001VB_filt,
                     elements = c('Ba_Ca_mav', 'Sr_Ca_mav'),
                     spanT = 0.08, spanV = 0.02,
                     n = 6000, open_end = TRUE, open_begin = TRUE)
```

    #> Tooth: 9530.43 to 45993.12 
    #>  Vert:  5354.199 to 9907.717 
    #>  Normalized distance: 0.6968619

![](06b_dtw_run_files/figure-gfm/sma001-multi-bleach-1.png)<!-- -->

``` r
mult1b2.1 <- align_multi(SMA001TB_filt %>% filter(true_dist <= 46000),
                       SMA001VB_filt,
                       elements = c('Mg_Ca_mav', 'Sr_Ca_mav'),
                       spanT = 0.08, spanV = 0.02,
                       n = 6000, open_end = TRUE, open_begin = TRUE)
```

    #> Tooth: 9530.43 to 45993.12 
    #>  Vert:  1163.581 to 8554.046 
    #>  Normalized distance: 0.5095142

![](06b_dtw_run_files/figure-gfm/sma001-multi-bleach-2.png)<!-- -->

``` r
mult1b3 <- align_multi(SMA001TB_filt %>% filter(true_dist <= 46000),
                         SMA001VB_filt,
                         elements = c('Mg_Ca_mav', 'Sr_Ca_mav', "Mn_Ca_mav"),
                         spanT = 0.08, spanV = 0.02,
                         n = 6000, open_end = TRUE, open_begin = TRUE)
```

    #> Tooth: 9530.43 to 45993.12 
    #>  Vert:  2597.896 to 9130.076 
    #>  Normalized distance: 0.6208782

![](06b_dtw_run_files/figure-gfm/sma001-multi-bleach-3.png)<!-- -->

``` r
mult1b4 <- align_multi(SMA001TB_filt %>% filter(true_dist <= 46000),
                       SMA001VB_filt,
                       elements = mavs1,
                       spanT = 0.08, spanV = 0.02,
                       n = 6000, open_end = TRUE, open_begin = TRUE)
```

    #> Tooth: 9530.43 to 45993.12 
    #>  Vert:  3502.263 to 6670.428 
    #>  Normalized distance: 0.7352295

![](06b_dtw_run_files/figure-gfm/sma001-multi-bleach-4.png)<!-- -->

``` r
mult1b6 <- align_multi(SMA001TB_filt %>% filter(true_dist <= 46000),
                       SMA001VB_filt,
                       elements = mavs2,
                       spanT = 0.08, spanV = 0.02,
                       n = 6000, open_end = TRUE, open_begin = TRUE)
```

    #> Tooth: 9530.43 to 45993.12 
    #>  Vert:  3505.143 to 9729.147 
    #>  Normalized distance: 0.8137965

![](06b_dtw_run_files/figure-gfm/sma001-multi-bleach-5.png)<!-- -->

``` r
# ahhh that's interesting!
```

### Multi: water

``` r
mult1h2 <- align_multi(SMA001T_filt, SMA001V_filt,
                      elements = c('Ba_Ca_mav', 'Sr_Ca_mav'),
                      spanT = 0.08, spanV = 0.02,
                      n = 6000, open_begin = TRUE, # lockV_start = 3266.09,
                      open_end = TRUE) # same result for both verts
```

    #> Tooth: 812.202 to 54901.41 
    #>  Vert:  6938.282 to 10656.56 
    #>  Normalized distance: 0.5913537

![](06b_dtw_run_files/figure-gfm/sma001-multi-water-1.png)<!-- -->

``` r
mult1h2.1 <- align_multi(SMA001T_filt, SMA001V_filt,
                         elements = c('Mg_Ca_mav', 'Sr_Ca_mav'),
                         spanT = 0.08, spanV = 0.02,
                         n = 6000, open_end = TRUE, open_begin = TRUE)
```

    #> Tooth: 812.202 to 54901.41 
    #>  Vert:  1679.128 to 10397.34 
    #>  Normalized distance: 0.6562104

![](06b_dtw_run_files/figure-gfm/sma001-multi-water-2.png)<!-- -->

``` r
mult1h3 <- align_multi(SMA001T_filt, SMA001V_filt,
                       elements = c('Mg_Ca_mav', 'Sr_Ca_mav', "Mn_Ca_mav"),
                       spanT = 0.08, spanV = 0.02,
                       n = 6000, open_end = TRUE, open_begin = TRUE)
```

    #> Tooth: 812.202 to 54901.41 
    #>  Vert:  3640.51 to 9037.911 
    #>  Normalized distance: 0.8012172

![](06b_dtw_run_files/figure-gfm/sma001-multi-water-3.png)<!-- -->

``` r
mult1h4 <- align_multi(SMA001T_filt, SMA001V_filt,
                      elements = mavs1,
                      # spanT = 0.025, spanV = 0.005,
                      spanT = 0.08, spanV = 0.02,
                      n = 6000, open_begin = TRUE,
                      open_end = TRUE) # but this one suggests longer
```

    #> Tooth: 812.202 to 54901.41 
    #>  Vert:  3628.99 to 4815.611 
    #>  Normalized distance: 0.8811592

![](06b_dtw_run_files/figure-gfm/sma001-multi-water-4.png)<!-- -->

``` r
mult1h6 <- align_multi(SMA001T_filt, SMA001V_filt,
                      elements = mavs2,
                      # spanT = 0.025, spanV = 0.005,
                      spanT = 0.08, spanV = 0.02,
                      n = 6000, open_begin = TRUE, # lockV_start = 3266.09,
                      open_end = TRUE)
```

    #> Tooth: 812.202 to 54901.41 
    #>  Vert:  3628.99 to 5100.746 
    #>  Normalized distance: 0.8263147

![](06b_dtw_run_files/figure-gfm/sma001-multi-water-5.png)<!-- -->

## SMA002

`SMA002TH_out` is made in `05c_fig_sma002.R`, not in 06a, so it’s
rebuilt here the same way 05c does it.

``` r
# rebuilt as in 05c_fig_sma002.R (water tooth, In holes removed)
SMA002TH <- subset(flag_holes(transsplit[[59]], trigger = 3, release = 1.5), !.hole)
SMA002TH_out <- profs(SMA002TH, el_all, "tooth", span = 0.01, window_size = 11,
                      show_In = TRUE, In_style = "points")
```

![](06b_dtw_run_files/figure-gfm/sma002-profile-1.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma002-profile-2.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma002-profile-3.png)<!-- -->

``` r
SMA002TH_out$plots$Ca #(never gets to the base, tooth incomplete)
```

![](06b_dtw_run_files/figure-gfm/sma002-profile-4.png)<!-- -->

``` r
align_els(SMA002T_filt, SMA002V_filt, elements = mavs2,
          spanT = 0.08, spanV = 0.02, n = 4000)
```

    #> Sr_Ca_mav 
    #>  Tooth: 112.326 to 43012.17 
    #>  Vert:  6837.476 to 7949.214 
    #>  Normalized distance: 0.5359942

![](06b_dtw_run_files/figure-gfm/sma002-els-1.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma002-els-2.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma002-els-3.png)<!-- -->

    #> Ba_Ca_mav 
    #>  Tooth: 112.326 to 43012.17 
    #>  Vert:  1699.289 to 5797.742 
    #>  Normalized distance: 0.1254943

![](06b_dtw_run_files/figure-gfm/sma002-els-4.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma002-els-5.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma002-els-6.png)<!-- -->

    #> Li_Ca_mav 
    #>  Tooth: 112.326 to 43012.17 
    #>  Vert:  5423.323 to 7473.989 
    #>  Normalized distance: 0.2563396

![](06b_dtw_run_files/figure-gfm/sma002-els-7.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma002-els-8.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma002-els-9.png)<!-- -->

    #> Mg_Ca_mav 
    #>  Tooth: 112.326 to 43012.17 
    #>  Vert:  3611.708 to 5869.746 
    #>  Normalized distance: 0.2720076

![](06b_dtw_run_files/figure-gfm/sma002-els-10.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma002-els-11.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma002-els-12.png)<!-- -->

    #> Mn_Ca_mav 
    #>  Tooth: 112.326 to 43012.17 
    #>  Vert:  334.098 to 2698.701 
    #>  Normalized distance: 0.202941

![](06b_dtw_run_files/figure-gfm/sma002-els-13.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma002-els-14.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma002-els-15.png)<!-- -->

    #> Zn_Ca_mav 
    #>  Tooth: 112.326 to 43012.17 
    #>  Vert:  639.393 to 7505.67 
    #>  Normalized distance: 0.2191346

![](06b_dtw_run_files/figure-gfm/sma002-els-16.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma002-els-17.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma002-els-18.png)<!-- -->

    #> $Sr_Ca_mav
    #> DTW alignment object
    #>  Alignment size (query x reference): 4000 x 4000
    #>  Call: dtw(x = tSmooth, y = vSmooth, step.pattern = pattern, keep.internals = TRUE, 
    #>  Call:     open.end = TRUE, open.begin = TRUE)
    #> 
    #> $Ba_Ca_mav
    #> DTW alignment object
    #>  Alignment size (query x reference): 4000 x 4000
    #>  Call: dtw(x = tSmooth, y = vSmooth, step.pattern = pattern, keep.internals = TRUE, 
    #>  Call:     open.end = TRUE, open.begin = TRUE)
    #> 
    #> $Li_Ca_mav
    #> DTW alignment object
    #>  Alignment size (query x reference): 4000 x 4000
    #>  Call: dtw(x = tSmooth, y = vSmooth, step.pattern = pattern, keep.internals = TRUE, 
    #>  Call:     open.end = TRUE, open.begin = TRUE)
    #> 
    #> $Mg_Ca_mav
    #> DTW alignment object
    #>  Alignment size (query x reference): 4000 x 4000
    #>  Call: dtw(x = tSmooth, y = vSmooth, step.pattern = pattern, keep.internals = TRUE, 
    #>  Call:     open.end = TRUE, open.begin = TRUE)
    #> 
    #> $Mn_Ca_mav
    #> DTW alignment object
    #>  Alignment size (query x reference): 4000 x 4000
    #>  Call: dtw(x = tSmooth, y = vSmooth, step.pattern = pattern, keep.internals = TRUE, 
    #>  Call:     open.end = TRUE, open.begin = TRUE)
    #> 
    #> $Zn_Ca_mav
    #> DTW alignment object
    #>  Alignment size (query x reference): 4000 x 4000
    #>  Call: dtw(x = tSmooth, y = vSmooth, step.pattern = pattern, keep.internals = TRUE, 
    #>  Call:     open.end = TRUE, open.begin = TRUE)

### Multi

``` r
mult2 <- align_multi(SMA002T_filt, SMA002V_filt,
                       elements = c('Ba_Ca_mav', 'Sr_Ca_mav'),
                       spanT = 0.08, spanV = 0.02,
                       n = 6000, open_begin = TRUE, # lockV_start = 3266.09,
                       open_end = TRUE) #hmm!
```

    #> Tooth: 112.326 to 43012.17 
    #>  Vert:  3965.967 to 4815.611 
    #>  Normalized distance: 0.7585141

![](06b_dtw_run_files/figure-gfm/sma002-multi-1.png)<!-- -->

``` r
# oh! it's calling it right before (note tooth incomplete)

mult2.1 <- align_multi(SMA002T_filt, SMA002V_filt,
                       elements = c('Mg_Ca_mav', 'Sr_Ca_mav'),
                       spanT = 0.08, spanV = 0.02,
                       n = 6000, open_begin = TRUE, # lockV_start = 3266.09,
                       open_end = TRUE)
```

    #> Tooth: 112.326 to 43012.17 
    #>  Vert:  3942.926 to 4559.277 
    #>  Normalized distance: 0.706019

![](06b_dtw_run_files/figure-gfm/sma002-multi-2.png)<!-- -->

``` r
mult3 <- align_multi(SMA002T_filt, SMA002V_filt,
                       elements = c('Mg_Ca_mav', 'Sr_Ca_mav', 'Mn_Ca_mav'),
                       spanT = 0.08, spanV = 0.02,
                       n = 6000, open_begin = TRUE, # lockV_start = 3266.09,
                       open_end = TRUE)
```

    #> Tooth: 112.326 to 43012.17 
    #>  Vert:  3614.588 to 5708.457 
    #>  Normalized distance: 0.7343895

![](06b_dtw_run_files/figure-gfm/sma002-multi-3.png)<!-- -->

``` r
mult4 <- align_multi(SMA002T_filt, SMA002V_filt,
            elements =  mavs1,
            spanT = 0.08, spanV = 0.02,
            n = 6000,  open_begin = TRUE,
            open_end = TRUE)
```

    #> Tooth: 112.326 to 43012.17 
    #>  Vert:  3940.045 to 4115.734 
    #>  Normalized distance: 0.9184684

![](06b_dtw_run_files/figure-gfm/sma002-multi-4.png)<!-- -->

``` r
mult6 <- align_multi(SMA002T_filt, SMA002V_filt,
                       elements =  mavs2,
                       spanT = 0.04, spanV = 0.01,
                       n = 6000,  open_begin = TRUE,
                       # open_begin = FALSE, lockV_start = 3800,
                       open_end = TRUE)
```

    #> Tooth: 112.326 to 43012.17 
    #>  Vert:  3942.926 to 4302.944 
    #>  Normalized distance: 0.9145874

![](06b_dtw_run_files/figure-gfm/sma002-multi-5.png)<!-- -->

``` r
# oh woah these all look like they match - not long though and they are def scrunched.
```

## Save for 06f

Only what 06f uses is saved: each `mult*` alignment carries two n × n
cost matrices (`keep = TRUE`), so saving all of them comes to ~8 GB.

``` r
save(mult1h6, SMA001T_filt, SMA001V_filt, mavs2,
     file = 'data/output/dtw_cache/06b.RData')
```

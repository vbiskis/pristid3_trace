06b · DTW run
================
Nikki Biskis
Rendered 26 Oct 2026 19:03

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
# rebuilt as in 05c_fig_sma001.R 
SMA001TH_out <- profs(SMA001T_filt, el_all, "tooth", span = 0.01, window_size = 11,
                      show_In = TRUE, In_style = "rect")
```

![](06b_dtw_run_files/figure-gfm/sma001-profw-1.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-profw-2.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-profw-3.png)<!-- -->

``` r
SMA001TH_out$plots$Ca #(never gets to the base, tooth incomplete)
```

![](06b_dtw_run_files/figure-gfm/sma001-profw-4.png)<!-- -->

``` r
SMA001TB_out <- profs(SMA001TB_filt, el_all, "tooth", span = 0.01, window_size = 11,
                      show_In = TRUE, In_style = "rect")
```

![](06b_dtw_run_files/figure-gfm/sma001-profb-1.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-profb-2.png)<!-- -->![](06b_dtw_run_files/figure-gfm/sma001-profb-3.png)<!-- -->

``` r
SMA001TB_out$plots$Ca #(never gets to the base, tooth incomplete)
```

![](06b_dtw_run_files/figure-gfm/sma001-profb-4.png)<!-- -->

``` r
# rep1
align_els(SMA001T_filt,
          SMA001V_filt, # Ba gets confused
          elements = names(span_vec),
          spanT = 0.02, spanV = 0.04, n = 2000)
```

    #> Ba_Ca_mav 
    #>  Tooth: 820.842 to 54875.49 
    #>  Vert:  7174.454 to 10656.56 
    #>  Normalized distance: 0.6549069

![](06b_dtw_run_files/figure-gfm/sma001-els-w-1.png)<!-- -->

    #> Li_Ca_mav 
    #>  Tooth: 820.842 to 54875.49 
    #>  Vert:  3548.345 to 4830.012 
    #>  Normalized distance: 0.6130704

![](06b_dtw_run_files/figure-gfm/sma001-els-w-2.png)<!-- -->

    #> Mg_Ca_mav 
    #>  Tooth: 820.842 to 54875.49 
    #>  Vert:  1975.783 to 10287.9 
    #>  Normalized distance: 0.4175162

![](06b_dtw_run_files/figure-gfm/sma001-els-w-3.png)<!-- -->

    #> Mn_Ca_mav 
    #>  Tooth: 820.842 to 54875.49 
    #>  Vert:  3821.96 to 8810.379 
    #>  Normalized distance: 0.5088213

![](06b_dtw_run_files/figure-gfm/sma001-els-w-4.png)<!-- -->

    #> Sr_Ca_mav 
    #>  Tooth: 820.842 to 54875.49 
    #>  Vert:  6990.125 to 10656.56 
    #>  Normalized distance: 0.4544548

![](06b_dtw_run_files/figure-gfm/sma001-els-w-5.png)<!-- -->

    #> Zn_Ca_mav 
    #>  Tooth: 820.842 to 54875.49 
    #>  Vert:  3585.787 to 10538.47 
    #>  Normalized distance: 0.4250766

![](06b_dtw_run_files/figure-gfm/sma001-els-w-6.png)<!-- -->

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
# rep2
align_els(SMA001TB_filt,
          SMA001VB_filt,
          elements = names(span_vec),
          spanT = 0.02, spanV = 0.04, n = 2000)
```

    #> Ba_Ca_mav 
    #>  Tooth: 12001.59 to 52879.56 
    #>  Vert:  921.648 to 2070.828 
    #>  Normalized distance: 0.6321098

![](06b_dtw_run_files/figure-gfm/sma001-els-b-1.png)<!-- -->

    #> Li_Ca_mav 
    #>  Tooth: 12001.59 to 52879.56 
    #>  Vert:  2747.664 to 4444.072 
    #>  Normalized distance: 0.4406668

![](06b_dtw_run_files/figure-gfm/sma001-els-b-2.png)<!-- -->

    #> Mg_Ca_mav 
    #>  Tooth: 12001.59 to 52879.56 
    #>  Vert:  3813.319 to 9804.031 
    #>  Normalized distance: 0.3342614

![](06b_dtw_run_files/figure-gfm/sma001-els-b-3.png)<!-- -->

    #> Mn_Ca_mav 
    #>  Tooth: 12001.59 to 52879.56 
    #>  Vert:  4441.192 to 9965.319 
    #>  Normalized distance: 0.3662405

![](06b_dtw_run_files/figure-gfm/sma001-els-b-4.png)<!-- -->

    #> Sr_Ca_mav 
    #>  Tooth: 12001.59 to 52879.56 
    #>  Vert:  1298.948 to 6889.319 
    #>  Normalized distance: 0.4536682

![](06b_dtw_run_files/figure-gfm/sma001-els-b-5.png)<!-- -->

    #> Zn_Ca_mav 
    #>  Tooth: 12001.59 to 52879.56 
    #>  Vert:  5331.158 to 9651.383 
    #>  Normalized distance: 0.3415761

![](06b_dtw_run_files/figure-gfm/sma001-els-b-6.png)<!-- -->

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
mult1b2 <- align_multi(SMA001TB_filt, SMA001VB_filt,
                      elements = c('Ba_Ca_mav', 'Sr_Ca_mav'),
                      spanT = 0.02, spanV = 0.04,
                      n = 6000, open_end = TRUE, open_begin = TRUE)
```

    #> Tooth: 12001.59 to 52879.56 
    #>  Vert:  5152.588 to 6434.256 
    #>  Normalized distance: 0.7031147

![](06b_dtw_run_files/figure-gfm/sma001-multi-b-1.png)<!-- -->

``` r
mult1b2.1 <- align_multi(SMA001TB_filt, SMA001VB_filt,
                        elements = c('Mg_Ca_mav', 'Sr_Ca_mav'),
                        spanT = 0.02, spanV = 0.04,
                        n = 6000, open_end = TRUE, open_begin = TRUE)
```

    #> Tooth: 12001.59 to 52879.56 
    #>  Vert:  4101.334 to 9083.993 
    #>  Normalized distance: 0.6313461

![](06b_dtw_run_files/figure-gfm/sma001-multi-b-2.png)<!-- -->

``` r
mult1b3 <- align_multi(SMA001TB_filt, SMA001VB_filt,
                       elements = c('Mg_Ca_mav', 'Sr_Ca_mav', "Mn_Ca_mav"),
                       spanT = 0.02, spanV = 0.04,
                       n = 6000, open_end = TRUE, open_begin = TRUE)
```

    #> Tooth: 12001.59 to 52879.56 
    #>  Vert:  4383.588 to 9092.634 
    #>  Normalized distance: 0.6614065

![](06b_dtw_run_files/figure-gfm/sma001-multi-b-3.png)<!-- -->

``` r
mult1b4 <- align_multi(SMA001TB_filt, SMA001VB_filt,
                       elements = mavs1, spanT = 0.02, spanV = 0.04,
                       n = 6000, open_end = TRUE, open_begin = TRUE)
```

    #> Tooth: 12001.59 to 52879.56 
    #>  Vert:  3482.102 to 9965.319 
    #>  Normalized distance: 0.8378342

![](06b_dtw_run_files/figure-gfm/sma001-multi-b-4.png)<!-- -->

``` r
mult1b6 <- align_multi(SMA001TB_filt, SMA001VB_filt,
                       elements = mavs2, spanT = 0.02, spanV = 0.04,
                       n = 6000, open_end = TRUE, open_begin = TRUE)
```

    #> Tooth: 12001.59 to 52879.56 
    #>  Vert:  4115.735 to 9472.814 
    #>  Normalized distance: 0.8311513

![](06b_dtw_run_files/figure-gfm/sma001-multi-b-5.png)<!-- -->

``` r
# ahhh that's interesting!
```

### Multi: water

``` r
mult1h2 <- align_multi(SMA001T_filt, SMA001V_filt,
                      elements = c('Ba_Ca_mav', 'Sr_Ca_mav'),
                      spanT = 0.02, spanV = 0.04,
                      n = 6000, open_begin = TRUE, # lockV_start = 3266.09,
                      open_end = TRUE) # same result for both verts
```

    #> Tooth: 820.842 to 54875.49 
    #>  Vert:  6946.922 to 10656.56 
    #>  Normalized distance: 0.6512319

![](06b_dtw_run_files/figure-gfm/sma001-multi-w-1.png)<!-- -->

``` r
mult1h2.1 <- align_multi(SMA001T_filt, SMA001V_filt,
                         elements = c('Mg_Ca_mav', 'Sr_Ca_mav'),
                         spanT = 0.02, spanV = 0.04,
                         n = 6000, open_end = TRUE, open_begin = TRUE)
```

    #> Tooth: 820.842 to 54875.49 
    #>  Vert:  711.397 to 1921.06 
    #>  Normalized distance: 0.744051

![](06b_dtw_run_files/figure-gfm/sma001-multi-w-2.png)<!-- -->

``` r
mult1h3 <- align_multi(SMA001T_filt, SMA001V_filt,
                       elements = c('Mg_Ca_mav', 'Sr_Ca_mav', "Li_Ca_mav"),
                       spanT = 0.02, spanV = 0.04,
                       n = 6000, open_end = TRUE, open_begin = TRUE)
```

    #> Tooth: 820.842 to 54875.49 
    #>  Vert:  3597.308 to 4847.293 
    #>  Normalized distance: 0.923183

![](06b_dtw_run_files/figure-gfm/sma001-multi-w-3.png)<!-- -->

``` r
mult1h4 <- align_multi(SMA001T_filt, SMA001V_filt,
                      elements = mavs1,
                      spanT = 0.02, spanV = 0.04,
                      n = 6000, open_begin = TRUE,
                      open_end = TRUE) # but this one suggests longer
```

    #> Tooth: 820.842 to 54875.49 
    #>  Vert:  3600.188 to 4847.293 
    #>  Normalized distance: 0.9389438

![](06b_dtw_run_files/figure-gfm/sma001-multi-w-4.png)<!-- -->

``` r
mult1h6 <- align_multi(SMA001T_filt, SMA001V_filt,
                      elements = mavs2,
                      spanT = 0.02, spanV = 0.04,
                      n = 6000, open_begin = TRUE, # lockV_start = 3266.09,
                      open_end = TRUE)
```

    #> Tooth: 820.842 to 54875.49 
    #>  Vert:  3600.188 to 5077.705 
    #>  Normalized distance: 0.8775092

![](06b_dtw_run_files/figure-gfm/sma001-multi-w-5.png)<!-- -->

## SMA002

``` r
align_els(SMA002T_filt, SMA002V_filt, elements = mavs2,
          spanT = 0.02, spanV = 0.04, n = 4000)
```

    #> Sr_Ca_mav 
    #>  Tooth: 112.326 to 43012.17 
    #>  Vert:  7122.611 to 7931.933 
    #>  Normalized distance: 0.6173225

![](06b_dtw_run_files/figure-gfm/sma002-els-1.png)<!-- -->

    #> Ba_Ca_mav 
    #>  Tooth: 112.326 to 43012.17 
    #>  Vert:  1857.697 to 5809.263 
    #>  Normalized distance: 0.2596527

![](06b_dtw_run_files/figure-gfm/sma002-els-2.png)<!-- -->

    #> Li_Ca_mav 
    #>  Tooth: 112.326 to 43012.17 
    #>  Vert:  5437.724 to 7473.989 
    #>  Normalized distance: 0.4044115

![](06b_dtw_run_files/figure-gfm/sma002-els-3.png)<!-- -->

    #> Mg_Ca_mav 
    #>  Tooth: 112.326 to 43012.17 
    #>  Vert:  3614.588 to 5884.147 
    #>  Normalized distance: 0.3622041

![](06b_dtw_run_files/figure-gfm/sma002-els-4.png)<!-- -->

    #> Mn_Ca_mav 
    #>  Tooth: 112.326 to 43012.17 
    #>  Vert:  336.978 to 1935.461 
    #>  Normalized distance: 0.3649337

![](06b_dtw_run_files/figure-gfm/sma002-els-5.png)<!-- -->

    #> Zn_Ca_mav 
    #>  Tooth: 112.326 to 43012.17 
    #>  Vert:  858.285 to 7462.469 
    #>  Normalized distance: 0.3516778

![](06b_dtw_run_files/figure-gfm/sma002-els-6.png)<!-- -->

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
                       spanT = 0.02, spanV = 0.04,
                       n = 6000, open_begin = TRUE, # lockV_start = 3266.09,
                       open_end = TRUE) #hmm!
```

    #> Tooth: 112.326 to 43012.17 
    #>  Vert:  3965.967 to 4351.907 
    #>  Normalized distance: 0.7793199

![](06b_dtw_run_files/figure-gfm/sma002-multi-1.png)<!-- -->

``` r
# oh! it's calling it right before (note tooth incomplete)

mult2.1 <- align_multi(SMA002T_filt, SMA002V_filt,
                       elements = c('Mg_Ca_mav', 'Sr_Ca_mav'),
                       spanT = 0.02, spanV = 0.04,
                       n = 6000, open_begin = TRUE, # lockV_start = 3266.09,
                       open_end = TRUE)
```

    #> Tooth: 112.326 to 43012.17 
    #>  Vert:  3954.446 to 4343.266 
    #>  Normalized distance: 0.7465818

![](06b_dtw_run_files/figure-gfm/sma002-multi-2.png)<!-- -->

``` r
mult3 <- align_multi(SMA002T_filt, SMA002V_filt,
                       elements = c('Mg_Ca_mav', 'Sr_Ca_mav', 'Mn_Ca_mav'),
                       spanT = 0.02, spanV = 0.04,
                       n = 6000, open_begin = TRUE, # lockV_start = 3266.09,
                       open_end = TRUE)
```

    #> Tooth: 112.326 to 43012.17 
    #>  Vert:  3603.067 to 4775.289 
    #>  Normalized distance: 0.7678079

![](06b_dtw_run_files/figure-gfm/sma002-multi-3.png)<!-- -->

``` r
mult4 <- align_multi(SMA002T_filt, SMA002V_filt,
                    elements =  mavs1, spanT = 0.02, spanV = 0.04,
                    n = 6000, open_begin = TRUE, open_end = TRUE)
```

    #> Tooth: 112.326 to 43012.17 
    #>  Vert:  3942.926 to 4127.255 
    #>  Normalized distance: 0.9314739

![](06b_dtw_run_files/figure-gfm/sma002-multi-4.png)<!-- -->

``` r
mult6 <- align_multi(SMA002T_filt, SMA002V_filt,
                       elements =  mavs2,
                       spanT = 0.02, spanV = 0.04,
                       n = 6000,  open_begin = TRUE,
                       open_end = TRUE)
```

    #> Tooth: 112.326 to 43012.17 
    #>  Vert:  3911.244 to 4219.42 
    #>  Normalized distance: 0.9237344

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

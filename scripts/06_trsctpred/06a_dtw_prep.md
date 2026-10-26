06a · DTW prep
================
Nikki Biskis
Rendered 24 Oct 2026 19:14

- [Check other ratios](#check-other-ratios)
- [Apply mavs](#apply-mavs)
- [Filters](#filters)
  - [SMA001 tooth (water)](#sma001-tooth-water)
  - [SMA001 tooth (bleach)](#sma001-tooth-bleach)
  - [SMA001 vert (water)](#sma001-vert-water)
  - [SMA001 vert (bleach)](#sma001-vert-bleach)
  - [SMA002 tooth (bleach)](#sma002-tooth-bleach)
  - [SMA002 tooth (water)](#sma002-tooth-water)
  - [SMA002 vert (bleach)](#sma002-vert-bleach)
- [Mav column names](#mav-column-names)
- [Structure comparison](#structure-comparison)
  - [Tooth 1](#tooth-1)
  - [Vert](#vert)
- [Save for 06b–06e](#save-for-06b06e)

- **Paper:** The potential of elemental analysis in determining origin
  of sawfish rostra
- **DOI:** tbd
- **Written:** 8 Jul 2026
- **Content:** dynamic time warping strategy

R markdown version of `06a_dtw_prep.R`. The objects made here are saved
to `data/output/dtw_cache/06a.RData` at the end, and 06b–06e load that
file.

``` r
source('helpers/functions/fxn_dtw.R')
source('helpers/functions/fxn_profile.R')
source('helpers/functions/fxn_mavs.R')
source('helpers/microHelp.R')
transsplit <- readRDS('data/output/transsplit.RDS')
names(transsplit)
```

    #>  [1] "055HTsctBase"  "055HTsctMid"   "055HTsctTip"   "3149"         
    #>  [5] "3150TrsctBase" "3150TrsctMid"  "3150TrsctTip"  "3151"         
    #>  [9] "3158TrscBase"  "3158TrscMid"   "3158TrscTip"   "4224Transect" 
    #> [13] "4228TsctBase"  "4228TsctMid"   "4228TsctTip"   "437TsctBase"  
    #> [17] "437TsctMid"    "437TsctTip"    "4462Transect"  "4468"         
    #> [21] "4469HBase"     "4469CrossT2"   "4469FullT2"    "4469HMid"     
    #> [25] "4469HTip"      "448HTsctBase"  "448HTsctMid"   "448HTsctTip"  
    #> [29] "480BTT"        "480HTT"        "480BVT"        "480HVT"       
    #> [33] "482BTT"        "482BVT"        "482HVT"        "483BTT"       
    #> [37] "483HTT"        "483BVT"        "483HVT"        "484BTT"       
    #> [41] "484HTT"        "484BVT"        "484HVT"        "485BTT"       
    #> [45] "485HTT"        "485BVT"        "485HVT"        "486BTT"       
    #> [49] "486HTT"        "486BVT"        "486HVT"        "NB028B"       
    #> [53] "NB028Transect" "SMA001BT-TR"   "SMA001HT"      "SMA001BV"     
    #> [57] "SMA001HV"      "SMA002BT-TR"   "SMA002HT"      "SMA002BV"

``` r
plot_els <- c('Sr_Ca', 'Ba_Ca', 'Mg_Ca', 'Li_Ca')
elsg <- els[!els %in% c("Cu_Ca", "Pb_Ca", 'Ni_Ca')]
SMA001TB <- subset(flag_holes(transsplit[[54]], trigger = 5, release = 2.5), !.hole)
SMA001T <- subset(flag_holes(transsplit[[55]], trigger = 3, release = 1.5), !.hole)
SMA001VB <- transsplit[[56]]
SMA001V <- transsplit[[57]]
SMA002TB <- subset(flag_holes(transsplit[[58]], trigger = 3, release = 1.5), !.hole)
SMA002T <- subset(flag_holes(transsplit[[59]], trigger = 3, release = 1.5), !.hole)
SMA002V <- transsplit[[60]]
```

## Check other ratios

``` r
plot_els_ba <- c('Sr_Ba', 'Mg_Ba', 'Li_Ba', 'Mn_Ba', 'Zn_Ba')
plot_els_mg <- c('Sr_Mg', 'Ba_Mg', 'Li_Mg', 'Mn_Mg', 'Zn_Mg')
denoms   <- c('Ba', 'Mg')
plot_els_xc <- c('Sr', 'Ba', 'Mg', 'Li', 'Mn', 'Zn') # feels inefficient w/e

dfs <- list(SMA001TB = SMA001TB, SMA001T = SMA001T, SMA001VB = SMA001VB, SMA001V = SMA001V,
            SMA002TB = SMA002TB, SMA002T = SMA002T, SMA002V = SMA002V)
dfs <- lapply(dfs, add_alt, elements = plot_els_xc, denoms = denoms)
list2env(dfs, envir = .GlobalEnv)
```

    #> <environment: R_GlobalEnv>

## Apply mavs

``` r
allmavs <- Reduce(union, list(elsg, plot_els_ba, plot_els_mg, 'In_Ca'))

SMA001TB <- add_mavs(SMA001TB, ratio_columns = allmavs, window_size = 11) # 43 um for all
SMA001T <- add_mavs(SMA001T, ratio_columns = allmavs, window_size = 11)
SMA001VB <- add_mavs(SMA001VB, ratio_columns = allmavs, window_size = 21)
SMA001V <- add_mavs(SMA001V, ratio_columns = allmavs, window_size = 21)

SMA002TB <- add_mavs(SMA002TB, ratio_columns = allmavs, window_size = 11)
SMA002T <- add_mavs(SMA002T, ratio_columns = allmavs, window_size = 11)
SMA002V <- add_mavs(SMA002V, ratio_columns = allmavs, window_size = 21)

names(SMA001TB)
```

    #>  [1] "SpecID"        "Species"       "Type"          "Loc"          
    #>  [5] "Bath"          "ID"            "ElapsedTime_s" "speed"        
    #>  [9] "distance"      "Ca43_CPS"      "IntStdWv"      "Li_Ca"        
    #> [13] "B_Ca"          "Na_Ca"         "Mg_Ca"         "Mn_Ca"        
    #> [17] "Ni_Ca"         "Cu_Ca"         "Zn_Ca"         "Sr_Ca"        
    #> [21] "In_Ca"         "Ba_Ca"         "Pb_Ca"         "Speed"        
    #> [25] "true_dist"     ".hole"         "Sr_Ba"         "Mg_Ba"        
    #> [29] "Li_Ba"         "Mn_Ba"         "Zn_Ba"         "Sr_Mg"        
    #> [33] "Ba_Mg"         "Li_Mg"         "Mn_Mg"         "Zn_Mg"        
    #> [37] "Li_Ca_mav"     "B_Ca_mav"      "Mg_Ca_mav"     "Mn_Ca_mav"    
    #> [41] "Zn_Ca_mav"     "Sr_Ca_mav"     "Ba_Ca_mav"     "Sr_Ba_mav"    
    #> [45] "Mg_Ba_mav"     "Li_Ba_mav"     "Mn_Ba_mav"     "Zn_Ba_mav"    
    #> [49] "Sr_Mg_mav"     "Ba_Mg_mav"     "Li_Mg_mav"     "Mn_Mg_mav"    
    #> [53] "Zn_Mg_mav"     "In_Ca_mav"

## Filters

### SMA001 tooth (water)

``` r
SMA001T_filt <- SMA001T %>% filter(true_dist >= 650, true_dist <= 55583) %>% na.omit()
make_plot(SMA001T_filt, 'Sr_Ca', color = "brown", # had to clip it 3000 extra due to contam
          'Sr_Ca', span = 0.05, show_In = TRUE)
```

![](06a_dtw_prep_files/figure-gfm/filt-sma001t-1.png)<!-- -->

### SMA001 tooth (bleach)

``` r
SMA001TB_filt <- SMA001TB %>% filter(true_dist >= 5000, # remove na (overcautious but to be safe)
                                     true_dist <= 55402) %>% na.omit()
make_plot(SMA001TB_filt, 'Sr_Ca', color = "brown",
          'Sr_Ca', span = 0.05, show_In = TRUE)
```

![](06a_dtw_prep_files/figure-gfm/filt-sma001tb-1.png)<!-- -->

### SMA001 vert (water)

``` r
SMA001V_filt <- SMA001V %>% filter(true_dist >= 710, true_dist <= 10659) %>% na.omit()
make_plot(SMA001V_filt, 'Sr_Ca', color = "brown",
          'Sr_Ca', span = 0.0125, show_In = TRUE, In_style = "points")
```

![](06a_dtw_prep_files/figure-gfm/filt-sma001v-1.png)<!-- -->

### SMA001 vert (bleach)

``` r
SMA001VB_filt <- SMA001VB %>% filter(true_dist >= 745, true_dist <= 9968) %>% na.omit() # from measurements
make_plot(SMA001VB_filt, 'Sr_Ca', color = "brown",
          'Sr_Ca', span = 0.0125, show_In = TRUE, In_style = "points")
```

![](06a_dtw_prep_files/figure-gfm/filt-sma001vb-1.png)<!-- -->

### SMA002 tooth (bleach)

``` r
SMA002TB_filt <- SMA002TB %>% filter(true_dist <= 49907) %>% na.omit() # yuck - toss
make_plot(SMA002TB_filt, 'Sr_Ca', color = "brown",
          'Sr_Ca', span = 0.05, show_In = TRUE, In_style = "points")
```

![](06a_dtw_prep_files/figure-gfm/filt-sma002tb-1.png)<!-- -->

### SMA002 tooth (water)

``` r
SMA002T_filt <- SMA002T %>% na.omit() # he's great tho!
make_plot(SMA002T_filt, 'Sr_Ca', color = "brown",
          'Sr_Ca', span = 0.05, show_In = TRUE, In_style = "points")
```

![](06a_dtw_prep_files/figure-gfm/filt-sma002t-1.png)<!-- -->

### SMA002 vert (bleach)

``` r
SMA002V_filt <- SMA002V %>% filter(true_dist >= 334, true_dist <= 7952) %>% na.omit() # from measurements
make_plot(SMA002V_filt, 'Sr_Ca', color = "brown",
          'Sr_Ca', span = 0.0125, show_In = TRUE, In_style = "points")
```

![](06a_dtw_prep_files/figure-gfm/filt-sma002v-1.png)<!-- -->

## Mav column names

``` r
# get the name for mavs cols!
names(SMA001TB_filt) # they arent nearly as good - lose too much detail with smoothing
```

    #>  [1] "SpecID"        "Species"       "Type"          "Loc"          
    #>  [5] "Bath"          "ID"            "ElapsedTime_s" "speed"        
    #>  [9] "distance"      "Ca43_CPS"      "IntStdWv"      "Li_Ca"        
    #> [13] "B_Ca"          "Na_Ca"         "Mg_Ca"         "Mn_Ca"        
    #> [17] "Ni_Ca"         "Cu_Ca"         "Zn_Ca"         "Sr_Ca"        
    #> [21] "In_Ca"         "Ba_Ca"         "Pb_Ca"         "Speed"        
    #> [25] "true_dist"     ".hole"         "Sr_Ba"         "Mg_Ba"        
    #> [29] "Li_Ba"         "Mn_Ba"         "Zn_Ba"         "Sr_Mg"        
    #> [33] "Ba_Mg"         "Li_Mg"         "Mn_Mg"         "Zn_Mg"        
    #> [37] "Li_Ca_mav"     "B_Ca_mav"      "Mg_Ca_mav"     "Mn_Ca_mav"    
    #> [41] "Zn_Ca_mav"     "Sr_Ca_mav"     "Ba_Ca_mav"     "Sr_Ba_mav"    
    #> [45] "Mg_Ba_mav"     "Li_Ba_mav"     "Mn_Ba_mav"     "Zn_Ba_mav"    
    #> [49] "Sr_Mg_mav"     "Ba_Mg_mav"     "Li_Mg_mav"     "Mn_Mg_mav"    
    #> [53] "Zn_Mg_mav"     "In_Ca_mav"

``` r
mavs1 <- c('Sr_Ca_mav', 'Ba_Ca_mav', 'Li_Ca_mav', 'Mg_Ca_mav')
mavs2 <- c('Sr_Ca_mav', 'Ba_Ca_mav', 'Li_Ca_mav', 'Mg_Ca_mav', 'Mn_Ca_mav', 'Zn_Ca_mav')
mavsb <- names(SMA001VB_filt)[43:47]
mavsm <- names(SMA001VB_filt)[48:51]
```

## Structure comparison

But be systematic.

``` r
spans <- c(0.1, 0.075, 0.05, 0.025, 0.01, 0.0075, 0.005) #0.1 is high but
 # and 0.0025 beyond plotting capability
ns <- c(500, 1000, 2000)
```

### Tooth 1

``` r
tooth_comp <- align_grid(SMA001TB_filt, SMA001T_filt, mavs2, spans, ns)
# quick function for summarising these
result <- select_spans(tooth_comp, n_fixed = 2000, r2_threshold = 0.9) # yes!
span_vec <- setNames(result$chosen$span, result$chosen$element) # for ca ratios

# median looks to be 0.01
```

#### Test it now

``` r
diff(range(SMA001TB_filt$true_dist, na.rm = TRUE)) # for comparison
```

    #> [1] 44601.99

``` r
result$consensus_end - result$consensus_start # better!
```

    #> [1] 31840.07

``` r
align_els(SMA001TB_filt, SMA001T_filt, elements = names(span_vec),
          spanT = 0.01, spanV = 0.01, n = 2000) # all start around 20! nice work
```

    #> Ba_Ca_mav 
    #>  Tooth: 9530.43 to 54132.42 
    #>  Vert:  4873.215 to 46425.15 
    #>  Normalized distance: 0.2507101

![](06a_dtw_prep_files/figure-gfm/tooth-test-1.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/tooth-test-2.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/tooth-test-3.png)<!-- -->

    #> Li_Ca_mav 
    #>  Tooth: 9530.43 to 54132.42 
    #>  Vert:  1157.82 to 33594.06 
    #>  Normalized distance: 0.1882462

![](06a_dtw_prep_files/figure-gfm/tooth-test-4.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/tooth-test-5.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/tooth-test-6.png)<!-- -->

    #> Mg_Ca_mav 
    #>  Tooth: 9530.43 to 54132.42 
    #>  Vert:  15423.2 to 45163.62 
    #>  Normalized distance: 0.2611469

![](06a_dtw_prep_files/figure-gfm/tooth-test-7.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/tooth-test-8.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/tooth-test-9.png)<!-- -->

    #> Mn_Ca_mav 
    #>  Tooth: 9530.43 to 54132.42 
    #>  Vert:  10247.57 to 54901.41 
    #>  Normalized distance: 0.2140595

![](06a_dtw_prep_files/figure-gfm/tooth-test-10.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/tooth-test-11.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/tooth-test-12.png)<!-- -->

    #> Sr_Ca_mav 
    #>  Tooth: 9530.43 to 54132.42 
    #>  Vert:  13384.06 to 45189.54 
    #>  Normalized distance: 0.2600666

![](06a_dtw_prep_files/figure-gfm/tooth-test-13.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/tooth-test-14.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/tooth-test-15.png)<!-- -->

    #> Zn_Ca_mav 
    #>  Tooth: 9530.43 to 54132.42 
    #>  Vert:  14585.08 to 50304.69 
    #>  Normalized distance: 0.1989351

![](06a_dtw_prep_files/figure-gfm/tooth-test-16.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/tooth-test-17.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/tooth-test-18.png)<!-- -->

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
finalt <- align_multi(SMA001T_filt, SMA001TB_filt,
                      elements = mavs1,
                      open_begin = TRUE, open_end = TRUE,
                      spanT = 0.05, spanV = 0.05,
                      n = 2000) # and this matches that number well too!
```

    #> Tooth: 812.202 to 54901.41 
    #>  Vert:  28116.03 to 46356.02 
    #>  Normalized distance: 0.5830888

![](06a_dtw_prep_files/figure-gfm/tooth-test-19.png)<!-- -->

``` r
range(SMA001T_filt$true_dist)
```

    #> [1]   812.202 54901.410

``` r
range(SMA001TB_filt$true_dist)
```

    #> [1]  9530.43 54132.42

### Vert

But most importantly: our verts.

``` r
vert_comp <- align_grid(SMA001V_filt, SMA001VB_filt, mavs2, spans, ns)
result <- select_spans(vert_comp, n_fixed = 2000, r2_threshold = 0.9) # awesome!
span_vecv <- setNames(result$chosen$span, result$chosen$element) # oh wow

# how much wigglier are we talking
span_vecv/span_vec # was span_vec2 in the .R script, which isn't defined anywhere in 06a
```

    #> Ba_Ca_mav Li_Ca_mav Mg_Ca_mav Mn_Ca_mav Sr_Ca_mav Zn_Ca_mav 
    #> 1.3333333 0.6666667 1.0000000 2.0000000 7.5000000 7.5000000

``` r
# that actually makes a lot of sense.

align_els(SMA001VB_filt, SMA001V_filt, elements = names(span_vec),
          spanT = 0.001, spanV = 0.001, n = 2000)
```

    #> Ba_Ca_mav 
    #>  Tooth: 745.959 to 9965.319 
    #>  Vert:  1538.001 to 8781.577 
    #>  Normalized distance: 0.01706575

![](06a_dtw_prep_files/figure-gfm/vert-grid-1.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/vert-grid-2.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/vert-grid-3.png)<!-- -->

    #> Li_Ca_mav 
    #>  Tooth: 745.959 to 9965.319 
    #>  Vert:  1218.304 to 10512.55 
    #>  Normalized distance: 0.1413239

![](06a_dtw_prep_files/figure-gfm/vert-grid-4.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/vert-grid-5.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/vert-grid-6.png)<!-- -->

    #> Mg_Ca_mav 
    #>  Tooth: 745.959 to 9965.319 
    #>  Vert:  1483.278 to 10616.23 
    #>  Normalized distance: 0.06997789

![](06a_dtw_prep_files/figure-gfm/vert-grid-7.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/vert-grid-8.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/vert-grid-9.png)<!-- -->

    #> Mn_Ca_mav 
    #>  Tooth: 745.959 to 9965.319 
    #>  Vert:  1831.776 to 9179.038 
    #>  Normalized distance: 0.1168543

![](06a_dtw_prep_files/figure-gfm/vert-grid-10.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/vert-grid-11.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/vert-grid-12.png)<!-- -->

    #> Sr_Ca_mav 
    #>  Tooth: 745.959 to 9965.319 
    #>  Vert:  1523.6 to 10650.8 
    #>  Normalized distance: 0.02836459

![](06a_dtw_prep_files/figure-gfm/vert-grid-13.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/vert-grid-14.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/vert-grid-15.png)<!-- -->

    #> Zn_Ca_mav 
    #>  Tooth: 745.959 to 9965.319 
    #>  Vert:  1820.255 to 10477.99 
    #>  Normalized distance: 0.2052073

![](06a_dtw_prep_files/figure-gfm/vert-grid-16.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/vert-grid-17.png)<!-- -->![](06a_dtw_prep_files/figure-gfm/vert-grid-18.png)<!-- -->

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
# plots look EXCELLENT!

finalv <- align_multi(SMA001VB_filt, SMA001V_filt,
                      elements = mavs1,
                      open_begin = TRUE, open_end = TRUE,
                      spanT = 0.01, spanV = 0.01,
                      n = 2000)
```

    #> Tooth: 745.959 to 9965.319 
    #>  Vert:  1463.117 to 10621.99 
    #>  Normalized distance: 0.1426421

![](06a_dtw_prep_files/figure-gfm/vert-grid-19.png)<!-- -->

``` r
# straight up... wow

range(SMA001VB_filt$true_dist)
```

    #> [1]  745.959 9965.319

``` r
range(SMA002T_filt$true_dist) # 4x; should match spans
```

    #> [1]   112.326 43012.170

## Save for 06b–06e

``` r
dir.create('data/output/dtw_cache', showWarnings = FALSE)
save(list = Filter(function(x) !is.function(get(x)), ls()),
     file = 'data/output/dtw_cache/06a.RData')
```

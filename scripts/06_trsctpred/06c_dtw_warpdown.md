06c · DTW warp down
================
Nikki Biskis
Rendered 27 Oct 2026 21:05

- [SMA001 (water)](#sma001-water)
- [SMA001 (bleach)](#sma001-bleach)
- [SMA002](#sma002)

- **Paper:** The potential of elemental analysis in determining origin
  of sawfish rostra
- **DOI:** tbd
- **Written:** 14 Aug 2026
- **Content:** dynamic time warping IV

R markdown version of `original/06d_dtw_warpdown.R`. It starts from the
objects 06a saved (`data/output/dtw_cache/06a.RData`).

``` r
source('helpers/microHelp.R')
source('helpers/functions/fxn_mavs.R')
source('helpers/functions/fxn_dtw.R')
source('helpers/functions/fxn_stack.R')
load('data/output/dtw_cache/06a.RData') # objects from 06a
```

Instead of building up sections, trialing cutting down the tooth

- and seeing which elements are aligned in a section
- and importantly, whether cutting it makes the tooth end to an earlier
  part of the vert

5 gives a general idea - quick 10 fits the end much better

## SMA001 (water)

``` r
s <- stack_alignments(SMA001T_filt, SMA001V_filt, # 1
                      elements = mavs2, # or mavs1 for just the main 4
                      n = 1000, n_steps = c(16, 8), # how high res the map is
                      tol_um = 600, ref_step = 1000, # how (detailed) the table prints
                      plot = "heat", # for continuous scale
                      bands = c(710, # a visual of how path relates to banding
                                1715, 2900, 5900, 
                                8500, 9175, 9600))
```

![](06c_dtw_warpdown_files/figure-gfm/sma001w-1.png)<!-- -->

    #>  tooth_ref n_agree n_elements     members n_cuts vert_median vert_lo vert_hi
    #>       1000       3          6    Ba+Li+Sr      8        6595  3528.7    7192
    #>       2000       3          6    Ba+Li+Sr      8        6844  3548.6    7381
    #>       3000       3          6    Ba+Li+Sr     16        6976  3528.7    7441
    #>       4000       2          6       Ba+Sr     24        7203  5141.4    7650
    #>       5000       2          6       Ba+Sr     32        7156  5335.0    8074
    #>       6000       2          6       Mg+Mn     40        3385  2583.0    5412
    #>       7000       2          6       Mg+Mn     40        3489  3096.3    5729
    #>       8000       2          6       Mg+Mn     48        3474  3230.0    5729
    #>       9000       2          6       Mg+Mn     56        3851  2889.6    5729
    #>      10000       0          6        none     NA          NA      NA      NA
    #>      11000       0          6        none     NA          NA      NA      NA
    #>      12000       0          6        none     NA          NA      NA      NA
    #>      13000       0          6        none     NA          NA      NA      NA
    #>      14000       0          6        none     NA          NA      NA      NA
    #>      15000       2          6       Sr+Zn     96        2653  2030.3   10517
    #>      16000       2          6       Sr+Zn    104        2798  1348.5   10517
    #>      17000       2          6       Mn+Sr    104        3499  1348.5   10517
    #>      18000       2          6       Mn+Sr    112        3901  1249.0   10517
    #>      19000       2          6       Mn+Zn    120        3489  2377.7    5789
    #>      20000       2          6       Mg+Sr    128        5122   890.6   10517
    #>      21000       2          6       Mg+Sr    128        5390   890.6   10517
    #>      22000       2          6       Mg+Sr    128        5709   890.6   10517
    #>      23000       2          6       Mg+Sr    128        5719   890.6   10517
    #>      24000       2          6       Mg+Sr    128        5729   910.5   10607
    #>      25000       2          6       Mg+Sr    128        5729   950.3   10607
    #>      26000       2          6       Mg+Sr    128        5744   950.3   10607
    #>      27000       2          6       Mg+Sr    128        5818  1000.1   10607
    #>      28000       2          6       Mg+Sr    128        6034  1020.0   10657
    #>      29000       2          6       Mg+Sr    128        6087  1059.8   10657
    #>      30000       2          6       Mg+Sr    128        6112  1069.8   10657
    #>      31000       2          6       Mg+Sr    128        6122  1169.3   10657
    #>      32000       2          6       Mg+Mn    128        5779  1268.9    8466
    #>      33000       2          6       Mg+Mn    128        5779  1288.8    8573
    #>      34000       2          6       Mg+Mn    128        5779  1288.8    9155
    #>      35000       2          6       Mg+Sr    128        6246  1497.9   10657
    #>      36000       2          6       Mg+Sr    112        6197  5808.4   10657
    #>      37000       2          6       Mg+Sr    112        6207  5808.4   10657
    #>      38000       2          6       Mg+Sr    112        6256  5923.8   10657
    #>      39000       2          6       Mg+Sr     96        6419  6007.5   10657
    #>      40000       2          6       Mg+Sr     96        6792  6117.0   10657
    #>      41000       2          6       Li+Sr     96        7638  3857.2   10657
    #>      42000       2          6       Li+Mg     80        7531  3857.2    8218
    #>      43000       2          6       Li+Mg     80        7756  3907.3    8314
    #>      44000       2          6       Mn+Zn     80        5052  4822.9   10278
    #>      45000       2          6       Ba+Zn     64       10560  5051.8   10647
    #>      46000       2          6       Sr+Zn     64       10072  5051.8   10657
    #>      47000       3          6    Ba+Sr+Zn     48       10567  5051.8   10657
    #>      48000       3          6    Ba+Sr+Zn     48       10617  5051.8   10657
    #>      49000       3          6    Ba+Sr+Zn     48       10647  5051.8   10657
    #>      50000       3          6    Ba+Sr+Zn     32       10647 10248.4   10657
    #>      51000       3          6    Ba+Sr+Zn     32       10647 10507.2   10657
    #>      52000       4          6 Ba+Mg+Sr+Zn     32       10647  9929.8   10657
    #>      53000       4          6 Ba+Mg+Sr+Zn     16       10647 10178.7   10657
    #>      54000       4          6 Ba+Mg+Sr+Zn     16       10647 10198.6   10657
    #>  uncert_um
    #>     3663.5
    #>     3832.7
    #>     3912.4
    #>     2508.7
    #>     2739.4
    #>     2829.0
    #>     2632.5
    #>     2498.7
    #>     2839.2
    #>         NA
    #>         NA
    #>         NA
    #>         NA
    #>         NA
    #>     8486.9
    #>     9168.7
    #>     9168.7
    #>     9268.2
    #>     3410.8
    #>     9626.6
    #>     9626.6
    #>     9626.6
    #>     9626.6
    #>     9696.3
    #>     9656.5
    #>     9656.5
    #>     9606.7
    #>     9636.6
    #>     9596.7
    #>     9586.8
    #>     9487.2
    #>     7197.6
    #>     7284.1
    #>     7866.3
    #>     9158.7
    #>     4848.1
    #>     4848.1
    #>     4732.8
    #>     4649.0
    #>     4539.5
    #>     6799.3
    #>     4360.3
    #>     4406.9
    #>     5455.2
    #>     5594.8
    #>     5604.7
    #>     5604.7
    #>     5604.7
    #>     5604.7
    #>      408.2
    #>      149.3
    #>      726.7
    #>      477.8
    #>      457.9

``` r
s$heat(by_element = TRUE)
```

![](06c_dtw_warpdown_files/figure-gfm/sma001w-indi-1.png)<!-- -->

``` r
s$plot(by_element = TRUE)
```

![](06c_dtw_warpdown_files/figure-gfm/sma001w-cuts-1.png)<!-- -->

a couple options here looks like: - seeing overlap between Mg, Sr, and
to lesser extent Mn and Zn from 14k - 50k - and Li, Ba at the earlier
part (4k - 35k). VERY promising

## SMA001 (bleach)

``` r
t <- stack_alignments(SMA001TB_filt, SMA001VB_filt, # removing end due to inc holes
                      elements = mavs2, # or mavs1 for just the main 4
                      n = 1000, n_steps = c(16, 8), # how high res the map is
                      tol_um = 600, ref_step = 1000, # how (detailed) the table prints
                      plot = "heat", # for continuous scale
                      bands = c(630, # birth
                                1210, 2400, 5600, 
                                8000, 8800, 9300))
```

![](06c_dtw_warpdown_files/figure-gfm/sma001b-1.png)<!-- -->

    #>  tooth_ref n_agree n_elements     members n_cuts vert_median vert_lo vert_hi
    #>      13000       2          6       Sr+Zn     16        2005  1355.0    5231
    #>      14000       3          6    Li+Sr+Zn     24        2426  1283.2    5877
    #>      15000       3          6    Li+Sr+Zn     32        2702  1066.0    5905
    #>      16000       3          6    Li+Sr+Zn     40        2712  1465.8    5951
    #>      17000       3          6    Li+Sr+Zn     48        2712  1175.0    5979
    #>      18000       3          6    Li+Sr+Zn     56        2758  1465.8    5997
    #>      19000       3          6    Li+Sr+Zn     64        2776  1050.4    6025
    #>      20000       3          6    Li+Sr+Zn     72        3007  1419.6    6052
    #>      21000       3          6    Li+Sr+Zn     80        3159  1419.6    6052
    #>      22000       2          6       Sr+Zn     88        3422  1558.1    6147
    #>      23000       2          6       Li+Zn     96        2758  2444.0    6486
    #>      24000       0          6        none     NA          NA      NA      NA
    #>      25000       2          6       Mn+Sr    112        4934  1705.2    5757
    #>      26000       0          6        none     NA          NA      NA      NA
    #>      27000       2          6       Li+Zn    128        4510  1253.5    7427
    #>      28000       3          6    Li+Sr+Zn    128        5785  1672.2    7630
    #>      29000       3          6    Li+Sr+Zn    128        5803  2757.8    7898
    #>      30000       3          6    Li+Sr+Zn    128        5812  2757.8    7907
    #>      31000       2          6       Li+Sr    128        5812  2757.8    8361
    #>      32000       2          6       Li+Sr    128        5812  2757.8    9001
    #>      33000       3          6    Li+Mn+Zn    128        6984  2757.8    9790
    #>      34000       3          6    Li+Mn+Zn    128        7086  2757.8    9845
    #>      35000       3          6    Li+Mn+Zn    128        7395  2757.8    9882
    #>      36000       2          6       Mg+Sr    128        5859  5434.1    7665
    #>      37000       2          6       Li+Sr    128        6653  2757.8    9901
    #>      38000       2          6       Li+Sr    128        7345  2757.8    9901
    #>      39000       3          6    Mn+Sr+Zn    112        8512  4668.1    9965
    #>      40000       3          6    Ba+Mn+Zn    112        9019   930.5    9965
    #>      41000       4          6 Ba+Mn+Sr+Zn     96        9351   930.5    9965
    #>      42000       3          6    Ba+Sr+Zn     96        9458   930.5    9956
    #>      43000       3          6    Ba+Sr+Zn     80        9513   930.5    9956
    #>      44000       2          6       Ba+Mn     80        9956   930.5    9965
    #>      45000       3          6    Ba+Mn+Zn     64        9956   930.5    9965
    #>      46000       3          6    Ba+Mn+Zn     64        9956  1281.2    9965
    #>      47000       3          6    Ba+Mn+Zn     48        9956  1281.2    9965
    #>      48000       2          6       Ba+Mn     48        9961  1281.2    9965
    #>      49000       4          6 Ba+Mg+Mn+Zn     32        9616  1290.4    9965
    #>      50000       4          6 Ba+Mg+Mn+Zn     32        9811  1382.7    9965
    #>      51000       4          6 Ba+Mg+Mn+Zn     16        9781  1382.7    9965
    #>      52000       4          6 Ba+Mg+Mn+Zn     16        9781  2003.4    9965
    #>  uncert_um
    #>       3876
    #>       4594
    #>       4839
    #>       4485
    #>       4804
    #>       4531
    #>       4974
    #>       4633
    #>       4633
    #>       4589
    #>       4042
    #>         NA
    #>       4052
    #>         NA
    #>       6174
    #>       5958
    #>       5140
    #>       5150
    #>       5603
    #>       6244
    #>       7032
    #>       7088
    #>       7124
    #>       2231
    #>       7143
    #>       7143
    #>       5297
    #>       9035
    #>       9035
    #>       9026
    #>       9026
    #>       9035
    #>       9035
    #>       8684
    #>       8684
    #>       8684
    #>       8675
    #>       8583
    #>       8583
    #>       7962

``` r
t$heat(by_element = TRUE)
```

![](06c_dtw_warpdown_files/figure-gfm/sma001b-indi-1.png)<!-- -->

``` r
t$plot(by_element = TRUE)
```

![](06c_dtw_warpdown_files/figure-gfm/sma001b-cuts-1.png)<!-- -->

Ooh and all of them 22k - 50k?

## SMA002

``` r
u <- stack_alignments(SMA002T_filt, SMA002V_filt, # sawfish 1 - bleach
                      elements = mavs2, # or mavs1 for just the main 4
                      n = 1000, n_steps = c(16, 8), # how high res the map is
                      tol_um = 600, ref_step = 1000, # how (detailed) the table prints
                      plot = "heat", # for continuous scale
                      bands = c(770, 1025, # birth and/or 1yr
                                1650, 2700, 3700,
                                4550, 7550))
```

![](06c_dtw_warpdown_files/figure-gfm/sma002-1.png)<!-- -->

    #>  tooth_ref n_agree n_elements  members n_cuts vert_median vert_lo vert_hi
    #>       1000       2          6    Ba+Zn      8      1012.5   879.5    1973
    #>       2000       2          6    Mn+Zn     16       993.5   349.3    6049
    #>       3000       0          6     none     NA          NA      NA      NA
    #>       4000       0          6     none     NA          NA      NA      NA
    #>       5000       0          6     none     NA          NA      NA      NA
    #>       6000       0          6     none     NA          NA      NA      NA
    #>       7000       2          6    Li+Sr     56      5727.2   852.4    7660
    #>       8000       2          6    Li+Sr     64      5845.3   852.4    7660
    #>       9000       0          6     none     NA          NA      NA      NA
    #>      10000       0          6     none     NA          NA      NA      NA
    #>      11000       2          6    Ba+Mg     88      3634.7  2521.8    4168
    #>      12000       2          6    Ba+Mg     96      3634.7  2529.4    4206
    #>      13000       2          6    Ba+Mg    104      3680.8  2224.5    4237
    #>      14000       3          6 Ba+Mg+Sr    112      3772.1   852.4    7829
    #>      15000       3          6 Ba+Mg+Sr    120      3794.8   872.1    7919
    #>      16000       3          6 Ba+Mg+Sr    128      3868.0   997.3    7919
    #>      17000       3          6 Ba+Mg+Sr    128      4076.9   997.3    7919
    #>      18000       3          6 Ba+Mg+Sr    128      4236.7   997.3    7919
    #>      19000       3          6 Ba+Mg+Sr    128      4236.9  1451.8    7919
    #>      20000       3          6 Ba+Mg+Sr    128      4252.2  2063.7    7919
    #>      21000       3          6 Ba+Mg+Sr    128      4328.4  2639.1    7919
    #>      22000       3          6 Ba+Mg+Sr    128      4477.0  3149.1    7919
    #>      23000       3          6 Mg+Sr+Zn    128      4259.8  3227.4    7919
    #>      24000       3          6 Mg+Sr+Zn    128      4305.5  3558.5    7919
    #>      25000       3          6 Mg+Sr+Zn    128      4355.6  3840.6    7949
    #>      26000       3          6 Ba+Mg+Sr    128      5060.2  4198.8    7949
    #>      27000       2          6    Ba+Sr    128      5479.4  4343.3    7949
    #>      28000       2          6    Ba+Sr    112      5479.4  4374.1    7949
    #>      29000       2          6    Ba+Sr    112      5593.8  4496.1    7949
    #>      30000       2          6    Ba+Zn    112      5479.4  5247.1    7945
    #>      31000       2          6    Ba+Zn     96      5678.5  5479.4    7747
    #>      32000       2          6    Li+Sr     96      7949.2  4694.3    7949
    #>      33000       2          6    Ba+Mg     80      5616.7  4823.9    6013
    #>      34000       2          6    Ba+Mg     80      5620.5  4823.9    6348
    #>      35000       2          6    Ba+Mg     64      5647.5  4846.8    6371
    #>      36000       2          6    Ba+Mg     64      5670.0  4846.8    6490
    #>      37000       2          6    Ba+Mg     48      5670.0  4862.0    6037
    #>      38000       2          6    Ba+Mg     48      5690.7  4862.0    6371
    #>      39000       2          6    Ba+Mg     32      5719.5  5662.4    5797
    #>      40000       3          6 Li+Sr+Zn     32      7918.7  7065.0    7949
    #>      41000       2          6    Ba+Mg     16      5776.7  5670.0    6371
    #>      42000       2          6    Ba+Mg     16      5776.7  5670.0    6510
    #>  uncert_um
    #>     1093.5
    #>     5700.1
    #>         NA
    #>         NA
    #>         NA
    #>         NA
    #>     6807.1
    #>     6807.1
    #>         NA
    #>         NA
    #>     1646.5
    #>     1677.0
    #>     2012.4
    #>     6976.4
    #>     7046.7
    #>     6921.4
    #>     6921.4
    #>     6921.4
    #>     6466.9
    #>     5855.0
    #>     5279.6
    #>     4769.7
    #>     4691.3
    #>     4360.2
    #>     4108.7
    #>     3750.4
    #>     3605.9
    #>     3575.1
    #>     3453.1
    #>     2698.4
    #>     2267.1
    #>     3254.9
    #>     1189.1
    #>     1524.5
    #>     1524.5
    #>     1642.9
    #>     1174.8
    #>     1509.3
    #>      135.0
    #>      884.2
    #>      701.3
    #>      839.7

``` r
u$heat(by_element = TRUE) 
```

![](06c_dtw_warpdown_files/figure-gfm/sma002-indi-1.png)<!-- --> Ba/Mg
10k-30k or Sr/Zn 19k-28k

``` r
u$plot(by_element = TRUE)
```

![](06c_dtw_warpdown_files/figure-gfm/sma002-cuts-1.png)<!-- -->

``` r
# for 06c2 (anchors): each fish's records, and its stack_alignments() run (paths, runs)
gap_um <- 350   # tooth holes wider than this (um), as in 06a's add_mavs()
fish <- list(
  SMA001W = list(tooth = SMA001T_filt,  vert = SMA001V_filt,  bands = c(710, 1715, 2900, 5900, 8500, 9175, 9600)),
  SMA001B = list(tooth = SMA001TB_filt, vert = SMA001VB_filt, bands = c(630, 1210, 2400, 5600, 8000, 8800, 9300)),
  SMA002W = list(tooth = SMA002T_filt,  vert = SMA002V_filt,  bands = c(770, 1025, 1650, 2700, 3700, 4550, 7550)))
stacks <- lapply(list(SMA001W = s, SMA001B = t, SMA002W = u), `[`, c("paths", "runs"))
save(fish, stacks, gap_um, file = 'data/output/dtw_cache/06c.RData')
```

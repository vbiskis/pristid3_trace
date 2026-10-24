# Details ----
#' dtw_match()
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 28 Aug 2026
#' Content: 
#' dynamic time warping III
#' testing whether a loop can calculate the best location in a sequence, rather than guess and check
#' building on span determined in dtw1 and prelim results in dtw2
#' -----------

source('helpers/functions/fxn_dtw.R')
source('helpers/functions/fxn_mavs.R')
source('helpers/functions/fxn_warps.R')
source('helpers/microHelp.R')

library(patchwork)
library(furrr)
library(progressr)

# first visual----
# forced align chose ~6k or ~3k depending on # elements
match_window(SMA001T_filt, SMA001V_filt, # option 1
             elements = c("Ba_Ca_mav", "Sr_Ca_mav"),
             v_from = 6940, win = 3710, ratio_range = c(3, 16), 
             n_ratios = 24, spanV = 0.02, plot = TRUE) # based on 6b runs

match_window(SMA001T_filt, SMA001V_filt, # option 1
             elements = c("Mg_Ca_mav", "Sr_Ca_mav"),
             v_from = 1680, win = 8720, ratio_range = c(3, 16), 
             n_ratios = 24, spanV = 0.02, plot = TRUE, 
             bands = c(710, 1715, 2900, 5900, 8500, 9175, 9600)) # based on 6b runs
# warped heavily near the band! but the line graph looks excellent

# lets go hunting!----
# start with SMA001----
## first - an ordered loop:----

range(SMA001T_filt$true_dist)
range(SMA001V_filt$true_dist)

starts  <- seq(1680, 10500, by = 100) 
# start after residency in aquarium, can always search forward once we find the anchor
# thought process is - we know its had these teeth already at least, better than getting stuck early on
results <- list()
search  <- c(-Inf, Inf)
for(v in starts){
  r <- match_window(SMA001T_filt, SMA001V_filt, elements = c("Mg_Ca_mav", "Sr_Ca_mav"),
                    v_from = v, win = 500, spanV = 0.025,
                    ratio_range = c(3, 12),
                    tooth_max = 50000, tooth_search = search)
  results[[length(results) + 1]] <- r
  if(!is.na(r$tooth_start)) search <- c(r$tooth_start, Inf)
}
results <- do.call(rbind, results)

### lets check the output----
ggplot(results, aes((v_from + v_to) / 2, (tooth_start + tooth_end) / 2, colour = disagree_pct, 
                    size = signal)) + geom_point() + theme_bw() +
  labs(x = "Vert (um)", y = "Matched tooth (um)",
       colour = "Elements\ndisagree (%)", size = "Signal") # ooh! that is pretty
# weak from 4k to 5.5k but lets...
# similar with and without Mg

### what does that warp look like there---
match_window(SMA001T_filt, SMA001V_filt, 
             elements = c("Mg_Ca_mav", "Sr_Ca_mav"),
             v_from = 1500, win = 1700, 
             ratio_range = c(3, 12), n_ratios = 24,
             spanV = 0.02, tooth_max = 45000, plot = TRUE) #checks out its just the start

## now check random order----
plan(multisession, workers = availableCores() - 1) # leave one core free bc idk thats what Dave says so it's law

handlers(handler_progress(format = ":spin :current/:total [:bar] :percent  elapsed :elapsed  ETA :eta"))
results <- list()
starts <- seq(min(SMA001V_filt$true_dist), max(SMA001V_filt$true_dist) - 1000, by = 100)
with_progress({
  p <- progressor(along = starts)
  results <- future_map_dfr(starts, function(v){
    r <- match_window(SMA001T_filt, SMA001V_filt, keep_best = 0.1,
                      elements = c("Mg_Ca_mav", "Sr_Ca_mav"),
                      v_from = v, win = 1000, spanV = 0.02, # best span from 6a
                      ratio_range = c(3, 12),
                      tooth_search = c(10000, 55000)) # keeps finding high element agreement at early tooth (even though signal is weak)
    p()        # tick the bar
    r
  }, .options = furrr_options(packages = c("dtw", "ggplot2")))
}) 

plan(sequential)  # back to normal when done
# omg its!

### sticky beak----
ggplot(results, aes((v_from + v_to) / 2, (tooth_start + tooth_end) / 2, colour = disagree_pct, 
                    size = signal)) + geom_point() + theme_bw() +
  labs(x = "Vert (um)", y = "Matched tooth (um)",
       colour = "Elements\ndisagree (%)", size = "Signal")

# that's looking like we've got a signal in the middle!
# go with an anchor (5.9-6.5k) and build **
# can redo from opposite end @ 3k and go out too and see if we get the same

# check our anchor
match_window(SMA001T_filt, SMA001V_filt, 
             elements = c("Mg_Ca_mav", "Sr_Ca_mav"),
             v_from = 1700, win = 1000, ratio_range = c(3, 20),
             spanV = 0.02,
             bands = c(710, 1715, 2900, 5900, 8500, 9175, 9600),
             plot = TRUE) 

# we've found our anchor----
anchor <- subset(results, v_from >= 1600 & v_from <= 2100) # ah should have done this actually
anchor$v_mid <- (anchor$v_from + anchor$v_to) / 2
anchor$t_mid <- (anchor$tooth_start + anchor$tooth_end) / 2
coef(lm(t_mid ~ v_mid, data = anchor))[2]   # observed slope (tooth um per vert um)
anchor$ratio                                # ratios the windows chose

# and build it out from a part I KNOW aligns----

a1 <- anchor[1, ]; aN <- anchor[nrow(anchor), ]
args <- list(datT = SMA001T_filt, datV = SMA001V_filt,
             elements = c("Mg_Ca_mav", "Sr_Ca_mav"), win = 500,
             spanV = 0.02, ratio_range = c(3, 16), tooth_max = 45000 
             ) # from the graph 

## forward: later vert must match later tooth
fwd <- list(); search <- c(aN$tooth_start, Inf)
for(v in seq(aN$v_from + 200, 10000, by = 200)){
  r <- do.call(match_window, c(args, list(v_from = v, tooth_search = search)))
  fwd[[length(fwd) + 1]] <- r
  if(!is.na(r$tooth_start)) search <- c(r$tooth_start, Inf)
}

## backward: earlier vert must match earlier tooth
bwd <- list(); search <- c(-Inf, a1$tooth_end)
for(v in rev(seq(min(SMA001V_filt$true_dist), a1$v_from - 200, by = 200))){
  r <- do.call(match_window, c(args, list(v_from = v, tooth_search = search)))
  bwd[[length(bwd) + 1]] <- r
  if(!is.na(r$tooth_end)) search <- c(-Inf, r$tooth_end)
}

built <- rbind(do.call(rbind, bwd), anchor[names(results)], do.call(rbind, fwd))

### moment of truth----
ggplot(built, aes((v_from + v_to) / 2, (tooth_start + tooth_end) / 2, colour = disagree_pct, 
                    size = signal)) + geom_point() + theme_bw() +
  labs(x = "Vert (um)", y = "Matched tooth (um)",
       colour = "Elements\ndisagree (%)", size = "Signal") 
# not always a strong signal but elements are aligned

### final plot!----
match_window(SMA001T_filt, SMA001V_filt, 
             elements = c("Sr_Ca_mav", "Li_Ca_mav", "Mg_Ca_mav", "Ba_Ca_mav"),
             v_from = 1500, win = 3500, ratio_range = c(3, 24), n_ratios = 24,
             spanV = 0.02, plot = TRUE, 
             bands = c(710, # birth
                       1715, 2900, 5900, # spread out for most; there are faint ones not included
                       8500, 9175, 9600))

# repeat for SMA002----
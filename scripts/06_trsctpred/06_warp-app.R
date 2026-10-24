# Details ----
#' 06e_register_run.R
#' Paper: Tooth vs vertebra elemental chronologies (sawfish)
#' Replaces the span/n guess-and-check in 06a-06d with:
#'   stage 1  linear registration across all elements (full cost surface)
#'   check    leave-one-element-out (fit of each element at the timeline set by the others)
#'   stage 2  banded DTW, band width chosen by held-out elements
#'   stage 3  per-element lag (Mg delay) against a timeline fitted without it
#' Needs the *_filt objects from the "filters" section of 06a (edge clipping).
#' Uses RAW ratio columns (not _mav): smoothing is done inside, in vert um.
#' -----------

source('helpers/functions/fxn_profile.R')
source('helpers/functions/fxn_app.R')
transsplit <- readRDS('data/output/transsplit.RDS')
names(transsplit)[54:60]

SMA001T <- transsplit[[55]] %>% 
  filter(true_dist >= 165, true_dist <= 55798) %>% na.omit()

SMA001V <- transsplit[[57]] %>% 
  filter(true_dist <= 10659.44) %>% na.omit()

# 0. look first ----
run_els_app(SMA001T, SMA001V,
            elements = c("Sr_Ca", "Ba_Ca", "Li_Ca", "Mg_Ca", "Mn_Ca"),
            bands = c(710, # birth
                      1715, 2900, 5900, # spread out for most; there are faint ones not included
                      8500, 9175, 9600))

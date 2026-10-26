# Details ----
#' dtw.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 14 Aug 2026
#' Content: 
#' dynamic time warping IV
#' instead of building up sections, trialling cutting down the tooth 
#' + until the same elements are aligned in a section 
#' + and then growing it
#' -----------

## try again----

source('helpers/functions/fxn_stack.R') # stack_alignments()

s <- stack_alignments(SMA001T_filt, SMA001V_filt, # sawfish 1
                      elements = mavs2, # or mavs1 for just the main 4
                      spanT = 0.05, spanV = 0.01, # adjust res again, path ~ same
                      n = 1000, n_steps = 25, # how high res the map is
                      tol_um = 800, ref_step = 1000, # how (detailed) the table prints
                      plot = "heat", # for continuous scale
                      bands = c(710, # a visual of how path relates to banding
                                1715, 2900, 5900, 
                                8500, 9175, 9600))

s$groups       # stats per tooth cut: which elements & where
s$plot()
s$heat(by_element = TRUE)

t <- stack_alignments(SMA001TB_filt, SMA001VB_filt, # sawfish 1 - bleach
                      elements = mavs2, # or mavs1 for just the main 4
                      spanT = 0.05, spanV = 0.01, # adjust res again, path ~ same
                      n = 1000, n_steps = 5, # how high res the map is
                      tol_um = 800, ref_step = 1000, # how (detailed) the table prints
                      plot = "heat", # for continuous scale
                      bands = c(630, # birth
                                1210, 2400, 5600, 
                                8000, 8800, 9300))

t$plot()
t$heat(by_element = TRUE)

u <- stack_alignments(SMA002T_filt, SMA002V_filt, # sawfish 1 - bleach
                      elements = mavs2, # or mavs1 for just the main 4
                      spanT = 0.05, spanV = 0.01, # adjust res again, path ~ same
                      n = 1000, n_steps = 5, # how high res the map is
                      tol_um = 800, ref_step = 1000, # how (detailed) the table prints
                      plot = "heat", # for continuous scale
                      bands = c(770, 1025, # birth and/or 1yr
                                1650, 2700, 3700,
                                4550, 7550))

u$plot()
u$heat(by_element = TRUE) 
#oh Ba is useful for this one - makes sense - it lived longer in the wild.

# Details ----
#' render_06.R
#' Content:
#' renders the 06a-06f .Rmd files to GitHub markdown (.md + figure folders), in order
#' (render_reports() in helpers/functions/fxn_render.R does the work)
#' + each file runs in its own fresh R session, so nothing leaks between them
#' + each file's old .md and figure folder are replaced, so no stale plots are left behind
#'   (if a render fails, the old figures are put back)
#' + 06a saves data/output/dtw_cache/06a.RData, which 06b-06f load; 06c (warp down) saves 06c.RData,
#'   06c2 (anchors) reads it and saves 06c2.RData, 06d (warp up) saves 06d.RData, 06f (plots) saves 06f.RData, and 06g (stats) uses them; 06e needs only 06a
#' + logs go to data/output/dtw_cache/logs/
#' from R (RStudio / VS Code), all of them:
#'   source('scripts/06_trsctpred/render_06.R')
#' from Terminal, any folder (caffeinate keeps the Mac awake while you're away):
#'   caffeinate -i Rscript ~/Docs/GitHub/phd_chapters/pristid3_trace/scripts/06_trsctpred/render_06.R          # all of them
#'   caffeinate -i Rscript ~/Docs/GitHub/phd_chapters/pristid3_trace/scripts/06_trsctpred/render_06.R 06d 06e  # just these (needs 06a.RData already)
#' -----------

local({

# where this script is: from Rscript (--file=) or from source() (ofile)
self <- sub('^--file=', '', grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE))
if(length(self) == 0) self <- sys.frame(1)$ofile
if(is.null(self)) stop("run the whole file with source() or Rscript, not line by line")
dir_06 <- dirname(normalizePath(self)) # the .Rmd files sit next to this script

# project root = the folder with the .Rproj (put back where you were when done)
old_wd <- setwd(dir_06)
on.exit(setwd(old_wd))
setwd(here::here())
source('helpers/functions/fxn_render.R', local = TRUE)

docs <- c('06a' = '06a_dtw_prep.Rmd',
          '06a2' = '06a2_smooth_test.Rmd',  # rolling median width: breakpoints at several widths
          '06b' = '06b_dtw_run.Rmd',
          '06c' = '06c_dtw_warpdown.Rmd',   # warp down: cutting the tooth, heat maps
          '06c2' = '06c2_dtw_anchor.Rmd',   # anchors from pair agreement (reuses 06c's runs)
          '06d' = '06d_dtw_warpup.Rmd',     # warp up: build from the anchor, random order within it
          '06e' = '06e_dtw_combos.Rmd',     # align_multi with every mix of elements containing Mg
          '06f' = '06f_dtw_plots.Rmd',      # where each method puts the tooth (needs 06c, 06d and 06e)
          '06g' = '06g_dtw_stats.Rmd')      # stats: warping and misfit
needs <- c('06a2' = '06a', '06b' = '06a', '06c' = '06a', '06c2' = '06c', '06d' = '06c2', '06e' = '06a', '06f' = '06e', '06g' = '06f')

pick <- setdiff(commandArgs(trailingOnly = TRUE), '--args')
if(length(pick) == 0) pick <- names(docs)

render_reports(dir_06, docs, needs, pick, log_dir = 'data/output/dtw_cache/logs')

})

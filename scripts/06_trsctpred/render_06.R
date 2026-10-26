# Details ----
#' render_06.R
#' Content:
#' renders the 06a-06f .Rmd files to GitHub markdown (.md + figure folders), in order
#' (render_reports() in helpers/functions/fxn_render.R does the work)
#' + each file runs in its own fresh R session, so nothing leaks between them
#' + each file's old .md and figure folder are replaced, so no stale plots are left behind
#'   (if a render fails, the old figures are put back)
#' + 06a saves data/output/dtw_cache/06a.RData, which 06b-06f load; 06c (warp down) saves 06c.RData,
#'   06d (warp up) saves 06d.RData, 06e (final table) saves 06e.RData, and 06f (figures) uses all three
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
          '06b' = '06b_dtw_run.Rmd',
          '06c' = '06c_dtw_warpdown.Rmd',   # warp down: cutting the tooth, sweeps, warp-down mix
          '06d' = '06d_dtw_warpup.Rmd',     # warp up: random-order + sequential windows, warp-up mix
          '06e' = '06e_dtw_rates.Rmd',      # final table: Mg matches and forced alignments
          '06f' = '06f_dtw_plots.Rmd')      # final figures
needs <- c('06b' = '06a', '06c' = '06a', '06d' = '06c', '06e' = '06d', '06f' = '06e')

pick <- setdiff(commandArgs(trailingOnly = TRUE), '--args')
if(length(pick) == 0) pick <- names(docs)

render_reports(dir_06, docs, needs, pick, log_dir = 'data/output/dtw_cache/logs')

})

# Details ----
#' render_05.R
#' Content:
#' renders the 05 .Rmd files to GitHub markdown (.md + figure folders)
#' (render_reports() in helpers/functions/fxn_render.R does the work)
#' + each file runs in its own fresh R session; old figure folders are replaced
#' + rendering 05c also re-saves its figures (figures/fig6, figures/supp/s4)
#' + logs go to data/output/report_logs/
#' from R (RStudio / VS Code), all of them:
#'   source('scripts/05_trsctvis/render_05.R')
#' from Terminal, any folder:
#'   Rscript ~/Docs/GitHub/phd_chapters/pristid3_trace/scripts/05_trsctvis/render_05.R           # all of them
#'   Rscript ~/Docs/GitHub/phd_chapters/pristid3_trace/scripts/05_trsctvis/render_05.R 05c_sma002 # just this one
#' -----------

local({

# where this script is: from Rscript (--file=) or from source() (ofile)
self <- sub('^--file=', '', grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE))
if(length(self) == 0) self <- sys.frame(1)$ofile
if(is.null(self)) stop("run the whole file with source() or Rscript, not line by line")
dir_05 <- dirname(normalizePath(self)) # the .Rmd files sit next to this script

# project root = the folder with the .Rproj (put back where you were when done)
old_wd <- setwd(dir_05)
on.exit(setwd(old_wd))
setwd(here::here())
source('helpers/functions/fxn_render.R', local = TRUE)

docs <- c('05c_sma001' = '05c_fig_sma001.Rmd',
          '05c_sma002' = '05c_fig_sma002.Rmd')

pick <- setdiff(commandArgs(trailingOnly = TRUE), '--args')
if(length(pick) == 0) pick <- names(docs)

render_reports(dir_05, docs, pick = pick, log_dir = 'data/output/report_logs')

})

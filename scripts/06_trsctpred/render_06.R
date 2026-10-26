# Details ----
#' render_06.R
#' Content:
#' renders the 06a-06f .Rmd files to GitHub markdown (.md + figure folders), in order
#' + each file runs in its own fresh R session, so nothing leaks between them
#' + 06a saves data/output/dtw_cache/06a.RData, which 06b-06e load; 06b saves 06b.RData for 06f
#' + logs go to data/output/dtw_cache/logs/
#' run from any folder (Terminal, so it keeps going while you're away):
#'   caffeinate -i Rscript ~/Docs/GitHub/phd_chapters/pristid3_trace/scripts/06_trsctpred/render_06.R          # all of them
#'   caffeinate -i Rscript ~/Docs/GitHub/phd_chapters/pristid3_trace/scripts/06_trsctpred/render_06.R 06d 06e  # just these (needs 06a.RData already)
#' -----------

# move to the project root (two folders up from this script), wherever it was called from
self <- sub('^--file=', '', grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE))
setwd(if(length(self) == 1) file.path(dirname(normalizePath(self)), '..', '..') else here::here())

docs <- c('06a' = '06a_dtw_prep.Rmd',
          '06b' = '06b_dtw_run.Rmd',
          '06c' = '06c_dtw_warpup.Rmd',
          '06d' = '06d_dtw_warpdown.Rmd',
          '06e' = '06e_dtw_rates.Rmd',
          '06f' = '06f_dtw_plots.Rmd')
needs <- c('06b' = '06a', '06c' = '06a', '06d' = '06a', '06e' = '06a', '06f' = '06b')

pick <- commandArgs(trailingOnly = TRUE)
if(length(pick) == 0) pick <- names(docs)

# pandoc ships inside RStudio; Rscript from Terminal doesn't find it on its own
if(Sys.getenv('RSTUDIO_PANDOC') == '')
  Sys.setenv(RSTUDIO_PANDOC = '/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools')

dir_06 <- 'scripts/06_trsctpred'
logs   <- 'data/output/dtw_cache/logs'
dir.create(logs, recursive = TRUE, showWarnings = FALSE)

status <- c()
for(k in pick){
  dep <- needs[k]
  if(!is.na(dep) && identical(status[dep], 'FAILED')){
    status[k] <- 'SKIPPED'
    message(k, ': skipped, ', dep, ' failed')
    next
  }
  log <- file.path(logs, paste0(k, '.log'))
  message(format(Sys.time(), '%H:%M'), '  ', k, ': rendering (log: ', log, ')')
  t0 <- Sys.time()
  code <- system(sprintf("Rscript -e %s > %s 2>&1", shQuote(sprintf(
    "rmarkdown::render('%s', envir = globalenv())", file.path(dir_06, docs[k]))), shQuote(log)))
  mins <- round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1)
  status[k] <- if(code == 0) 'ok' else 'FAILED'
  message(format(Sys.time(), '%H:%M'), '  ', k, ': ', status[k], ' (', mins, ' min)')
}

print(data.frame(doc = names(status), status = status), row.names = FALSE)

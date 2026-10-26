# Details ----
#' fxn_render.R
#' Content:
#' render_reports(): renders .Rmd files to GitHub markdown (.md + figure folders), in order
#' + each file runs in its own fresh R session, so nothing leaks between them
#' + each file's old .md and figure folder are replaced, so no stale plots are left behind
#'   (if a render fails, the old figures are put back)
#' + a file is skipped if one it needs (needs = c(file = 'needed file')) failed or was skipped
#' + run from the project root; logs go to log_dir
#' -----------

render_reports <- function(dir, docs, needs = character(0), pick = names(docs), log_dir){

  if(!all(pick %in% names(docs)))
    stop("unknown: ", paste(setdiff(pick, names(docs)), collapse = ", "),
         " (use any of ", paste(names(docs), collapse = " "), ")")

  # pandoc ships inside RStudio (where depends on the RStudio version); Rscript from Terminal doesn't find it on its own
  if(Sys.getenv('RSTUDIO_PANDOC') == '' && Sys.which('pandoc') == ''){
    tools <- '/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools'
    found <- file.path(tools, c(Sys.info()[['machine']], 'aarch64', 'x86_64', '.'), 'pandoc')
    found <- found[file.exists(found)]
    if(length(found) == 0) stop("can't find pandoc: install it (pandoc.org) or RStudio")
    Sys.setenv(RSTUDIO_PANDOC = dirname(found[1]))
  }

  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  status <- c()
  for(k in pick){
    dep <- needs[k]
    if(!is.na(dep) && isTRUE(status[dep] %in% c('FAILED', 'SKIPPED'))){
      status[k] <- 'SKIPPED'
      message(k, ': skipped, ', dep, ' failed')
      next
    }
    rmd  <- file.path(dir, docs[k])
    figs <- sub('\\.Rmd$', '_files', rmd)  # old figures set aside, so none go stale,
    old  <- paste0(figs, '_old')            # and put back if the render fails
    unlink(old, recursive = TRUE)
    if(dir.exists(figs)) file.rename(figs, old)
    log <- file.path(log_dir, paste0(k, '.log'))
    message(format(Sys.time(), '%H:%M'), '  ', k, ': rendering (log: ', log, ')')
    t0 <- Sys.time()
    code <- system(sprintf("Rscript -e %s > %s 2>&1", shQuote(sprintf(
      "rmarkdown::render('%s', envir = globalenv())", rmd)), shQuote(log)))
    mins <- round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1)
    status[k] <- if(code == 0) 'ok' else 'FAILED'
    if(code == 0) unlink(old, recursive = TRUE) else if(dir.exists(old)){
      unlink(figs, recursive = TRUE)
      file.rename(old, figs)
    }
    message(format(Sys.time(), '%H:%M'), '  ', k, ': ', status[k], ' (', mins, ' min)')
  }

  print(data.frame(doc = names(status), status = status), row.names = FALSE)
  invisible(status)
}

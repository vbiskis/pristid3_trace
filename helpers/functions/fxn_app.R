# Details ----
#' dtw_match()
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 10 Oct 2026
#' Content: 
#' dynamic time warping IV
#' building on dtw3 output
#' -----------

library(shiny)
library(dtw)
library(ggplot2)

# runs of In:Ca above threshold on the tooth
find_holes <- function(dat, in_col = "In_Ca", threshold = 5){
  d  <- dat[order(dat$true_dist), ]
  hi <- !is.na(d[[in_col]]) & d[[in_col]] > threshold
  if(!any(hi)) return(data.frame(start_dist = numeric(0), end_dist = numeric(0)))
  r  <- rle(hi); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  data.frame(start_dist = d$true_dist[starts[r$values]],
             end_dist   = d$true_dist[ends[r$values]])
}

in_holes <- function(x, holes, buffer){
  out <- rep(FALSE, length(x))
  for(k in seq_len(nrow(holes)))
    out <- out | (x >= holes$start_dist[k] - buffer & x <= holes$end_dist[k] + buffer)
  out
}

# loess on kept rows -> grid -> NA in holes -> z-score
smooth_grid <- function(dat, el, span, grid, grid_hole){
  d  <- dat[!dat$.hole & is.finite(dat[[el]]), ]
  sm <- predict(loess(d[[el]] ~ d$true_dist, span = span))
  y  <- approx(d$true_dist, sm, xout = grid, ties = mean, rule = 2)$y
  y[grid_hole] <- NA
  (y - mean(y, na.rm = TRUE)) / sd(y, na.rm = TRUE)
}

# standard DTW; hole rows cost 0
align_one <- function(Tz, Vz, open_begin = TRUE, open_end = TRUE){
  C   <- abs(outer(Tz, Vz, "-"))
  gap <- is.na(Tz)
  C[gap, ] <- 0
  a   <- dtw(C, step.pattern = asymmetric, open.begin = open_begin, open.end = open_end)
  kept <- !gap[a$index1]
  list(a = a, gap = gap,
       fit = mean(C[cbind(a$index1, a$index2)][kept]))   # mean |dz| over kept tooth points
}

band_layer <- function(bands, width, axis = c("x", "y")){
  axis <- match.arg(axis)
  if(is.null(bands) || length(bands) == 0) return(NULL)
  if(is.null(width) || is.na(width) || width <= 0){
    if(axis == "x") return(geom_vline(xintercept = bands, colour = "darkorange", linewidth = 0.4))
    return(geom_hline(yintercept = bands, colour = "darkorange", linewidth = 0.4))
  }
  df <- if(axis == "x")
    data.frame(xmin = bands - width / 2, xmax = bands + width / 2, ymin = -Inf, ymax = Inf) else
      data.frame(xmin = -Inf, xmax = Inf, ymin = bands - width / 2, ymax = bands + width / 2)
  geom_rect(data = df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "darkorange", alpha = 0.3, inherit.aes = FALSE)
}

# app ----
run_els_app <- function(tooth, vert, elements, bands = NULL,
                        in_col = "In_Ca", spike_fun = NULL){
  
  tooth <- tooth[order(tooth$true_dist), ]
  vert  <- vert[order(vert$true_dist), ]
  vert$.hole <- FALSE
  tr <- floor(range(tooth$true_dist)); vr <- floor(range(vert$true_dist))
  
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(width = 3,
                   checkboxGroupInput("els", "Elements", elements, elements),
                   sliderInput("tclip", "Tooth kept (um)", tr[1], tr[2] + 1, tr, step = 10),
                   sliderInput("vclip", "Vert kept (um)",  vr[1], vr[2] + 1, vr, step = 10),
                   actionButton("best", "Best vert for selected elements"),
                   fluidRow(
                     column(6, numericInput("spanT", "Tooth span", 0.04, min = 0.001, step = 0.005)),
                     column(6, numericInput("spanV", "Vert span",  0.01, min = 0.001, step = 0.005))),
                   numericInput("n", "Points per series", 1000, min = 200, step = 100),
                   checkboxInput("ob", "Open start", TRUE),
                   checkboxInput("oe", "Open end",   TRUE),
                   hr(),
                   checkboxInput("holes", "Remove tooth holes (In:Ca)", TRUE),
                   fluidRow(
                     column(6, numericInput("thr", "In:Ca above", 5, min = 0, step = 0.5)),
                     column(6, numericInput("buf", "Buffer (um)", 50, min = 0, step = 10))),
                   hr(),
                   checkboxInput("show_bands", "Show vert bands", !is.null(bands)),
                   numericInput("bw", "Band width (um; 0 = line)", 50, min = 0, step = 10),
                   sliderInput("nl", "Match lines", 10, 200, 40, step = 10)
      ),
      mainPanel(width = 9,
                plotOutput("paths", height = "380px"),
                uiOutput("two_ui"),
                tableOutput("tab"),
                tableOutput("pair"))
    ))
  
  server <- function(input, output, session){
    clips <- debounce(reactive(list(t = input$tclip, v = input$vclip)), 400)
    spans <- debounce(reactive(list(t = input$spanT, v = input$spanV, n = input$n)), 600)
    holeset <- debounce(reactive(list(on = input$holes, thr = input$thr, buf = input$buf)), 500)
    
    prep <- reactive({
      validate(need(length(input$els) > 0, "Select at least one element"))
      cl <- clips(); sp <- spans(); hs <- holeset()
      tc <- tooth[tooth$true_dist >= cl$t[1] & tooth$true_dist <= cl$t[2], ]
      vc <- vert[vert$true_dist  >= cl$v[1] & vert$true_dist  <= cl$v[2], ]
      validate(need(nrow(tc) > 30 && nrow(vc) > 30, "Too little data in the kept ranges"))
      
      holes <- data.frame(start_dist = numeric(0), end_dist = numeric(0))
      if(isTRUE(hs$on)){
        holes <- tryCatch(
          if(is.null(spike_fun)) find_holes(tc, in_col, hs$thr) else spike_fun(tc, threshold = hs$thr),
          error = function(e) holes)
        holes <- holes[, c("start_dist", "end_dist")]
      }
      tc$.hole <- in_holes(tc$true_dist, holes, if(isTRUE(hs$on)) hs$buf else 0)
      vc$.hole <- FALSE
      
      tg <- seq(min(tc$true_dist), max(tc$true_dist), length.out = sp$n)
      vg <- seq(min(vc$true_dist), max(vc$true_dist), length.out = sp$n)
      tg_hole <- in_holes(tg, holes, if(isTRUE(hs$on)) hs$buf else 0)
      
      out <- tryCatch({
        Tz <- sapply(input$els, function(el) smooth_grid(tc, el, sp$t, tg, tg_hole))
        Vz <- sapply(input$els, function(el) smooth_grid(vc, el, sp$v, vg, rep(FALSE, sp$n)))
        list(Tz = matrix(Tz, ncol = length(input$els), dimnames = list(NULL, input$els)),
             Vz = matrix(Vz, ncol = length(input$els), dimnames = list(NULL, input$els)))
      }, error = function(e) NULL)
      validate(need(!is.null(out), "Smoothing failed - span too small for the kept data"))
      
      c(out, list(tg = tg, vg = vg, holes = holes, buf = hs$buf,
                  pct_removed = 100 * mean(tc$.hole)))
    })
    
    al <- reactive({
      p <- prep()
      withProgress(message = "DTW", value = 0, {
        setNames(lapply(colnames(p$Tz), function(el){
          incProgress(1 / ncol(p$Tz))
          align_one(p$Tz[, el], p$Vz[, el], input$ob, input$oe)
        }), colnames(p$Tz))
      })
    })
    
    # matched vert position per tooth point, per element
    paths <- reactive({
      p <- prep(); A <- al()
      do.call(rbind, lapply(names(A), function(el){
        a <- A[[el]]$a
        data.frame(element = el,
                   tooth   = p$tg[a$index1],
                   vert    = p$vg[a$index2],
                   gap     = A[[el]]$gap[a$index1])
      }))
    })
    
    bands_now <- reactive(if(isTRUE(input$show_bands)) bands else NULL)
    
    observeEvent(input$best, {
      p   <- prep()                                     # current tooth cut, holes, spans
      els <- colnames(p$Tz)
      vgf <- seq(min(vert$true_dist), max(vert$true_dist), length.out = length(p$tg))
      Vzf <- sapply(els, function(el) smooth_grid(vert, el, input$spanV, vgf, rep(FALSE, length(vgf))))
      Vzf <- matrix(Vzf, ncol = length(els), dimnames = list(NULL, els))
      
      # one cost matrix: mean |dz| across selected elements; holes cost 0
      C <- Reduce(`+`, lapply(els, function(el){
        Ci <- abs(outer(p$Tz[, el], Vzf[, el], "-"))
        Ci[is.na(p$Tz[, el]), ] <- 0
        Ci
      })) / length(els)
      
      a  <- dtw(C, step.pattern = asymmetric, open.begin = TRUE, open.end = TRUE)
      vs <- vgf[min(a$index2)]; ve <- vgf[max(a$index2)]
      
      updateSliderInput(session, "vclip", value = c(floor(vs), ceiling(ve)))
      updateCheckboxInput(session, "ob", value = FALSE)   # match the chosen span exactly
      updateCheckboxInput(session, "oe", value = FALSE)
      showNotification(sprintf("Vert %.0f-%.0f um; mean |dz| = %.3f", vs, ve,
                               mean(C[cbind(a$index1, a$index2)][!is.na(p$Tz[a$index1, 1])])))
    })
    
    output$paths <- renderPlot({
      p <- prep(); d <- paths()
      d$vert_kept <- ifelse(d$gap, NA, d$vert)
      hole_df <- if(nrow(p$holes) > 0)
        data.frame(xmin = p$holes$start_dist - p$buf, xmax = p$holes$end_dist + p$buf) else NULL
      ggplot(d, aes(tooth, colour = element)) +
        { if(!is.null(hole_df)) geom_rect(data = hole_df, aes(xmin = xmin, xmax = xmax,
                                                              ymin = -Inf, ymax = Inf),
                                          fill = "grey90", inherit.aes = FALSE) } +
        band_layer(bands_now(), input$bw, "y") +
        geom_line(aes(y = vert), linetype = "dashed", linewidth = 0.3) +
        geom_line(aes(y = vert_kept), linewidth = 0.9) +
        theme_bw() +
        labs(x = "Tooth (µm)", y = "Vertebra (µm)", colour = NULL)
    })
    
    output$two_ui <- renderUI(plotOutput("two", height = paste0(200 * length(input$els), "px")))
    
    output$two <- renderPlot({
      p <- prep(); A <- al(); off <- 4
      parts <- lapply(names(A), function(el){
        a  <- A[[el]]$a; gap <- A[[el]]$gap
        v_s <- p$vg[min(a$index2)]; v_e <- p$vg[max(a$index2)]
        tx <- v_s + (p$tg - p$tg[1]) / (p$tg[length(p$tg)] - p$tg[1]) * (v_e - v_s)
        keep <- which(!gap[a$index1])
        pick <- keep[unique(round(seq(1, length(keep), length.out = min(input$nl, length(keep)))))]
        list(
          vert  = data.frame(element = el, x = p$vg, y = p$Vz[, el]),
          tooth = data.frame(element = el, x = tx,   y = p$Tz[, el] + off),
          links = data.frame(element = el,
                             x = tx[a$index1[pick]],   y = p$Tz[a$index1[pick], el] + off,
                             xend = p$vg[a$index2[pick]], yend = p$Vz[a$index2[pick], el]))
      })
      vdf <- do.call(rbind, lapply(parts, `[[`, "vert"))
      tdf <- do.call(rbind, lapply(parts, `[[`, "tooth"))
      ldf <- do.call(rbind, lapply(parts, `[[`, "links"))
      ggplot() +
        band_layer(bands_now(), input$bw, "x") +
        geom_segment(data = ldf, aes(x, y, xend = xend, yend = yend),
                     colour = "grey60", linetype = "dashed", linewidth = 0.3) +
        geom_line(data = vdf, aes(x, y), colour = "black") +
        geom_line(data = tdf, aes(x, y), colour = "brown") +
        facet_wrap(~ element, ncol = 1) +
        theme_bw() +
        labs(x = "Vertebra (µm)", y = "z (tooth offset +4)")
    })
    
    output$tab <- renderTable({
      p <- prep(); A <- al()
      data.frame(
        element       = names(A),
        vert_start    = sapply(A, function(x) p$vg[min(x$a$index2)]),
        vert_end      = sapply(A, function(x) p$vg[max(x$a$index2)]),
        mean_abs_dz   = sapply(A, `[[`, "fit"),
        norm_dist     = sapply(A, function(x) x$a$normalizedDistance),
        tooth_removed_pct = p$pct_removed)
    }, digits = 3)
    
    output$pair <- renderTable({
      d <- paths()
      validate(need(length(unique(d$element)) > 1, ""))
      w <- reshape(d[!d$gap, c("element", "tooth", "vert")], idvar = "tooth",
                   timevar = "element", direction = "wide")
      E <- sub("^vert\\.", "", names(w)[-1])
      D <- outer(E, E, Vectorize(function(i, j)
        mean(abs(w[[paste0("vert.", i)]] - w[[paste0("vert.", j)]]), na.rm = TRUE)))
      dimnames(D) <- list(E, E)
      D
    }, rownames = TRUE, digits = 0)
  }
  
  shinyApp(ui, server)
}

# Early bacterial proliferation: two groups compared with a GAM of log2 load over time
# Input: the per-fly table saved by CFU Drop Counter (Load step), one row per fly, with the
# inoculum tubes of the Time 0 plate at time 0, as in the protocol ("Early proliferation").
# The analysis itself is in prolif_core.R, which the R code tab gives with the settings chosen.
# Runs as a normal Shiny app, or in the browser with shinylive (no server).

library(shiny)
library(bslib)
library(mgcv)
library(ggplot2)
library(patchwork)
if (FALSE) { library(svglite); library(jsonlite) }   # declares them to shinylive, which bundles the packages it finds in library() calls
has_svglite <- requireNamespace("svglite", quietly = TRUE)
`%||%` <- function(a, b) if (is.null(a)) b else a
source("prolif_core.R")

APP_VERSION <- "1.0"
EXAMPLE_FILE <- "example_proliferation_CS_sexdim.csv"
EXAMPLE_TIME0 <- "example_time0_CS_sexdim.csv"
N_SIM <- 10000          # protocol, step 4
CONTROL_FRAC <- 0.1     # protocol, step 6d
PALETTE <- c("#1f5fbf", "#e8850c", "#2a9d6f", "#c0398a", "#7a5cc4", "#8c6d1f", "#3a3a3a", "#d1453b")
# Female orange and male blue in every app of the protocol, whatever the order
group_colours <- function(lv) {
  cols <- setNames(rep(NA_character_, length(lv)), lv)
  low <- tolower(lv)
  cols[low %in% c("female", "females", "f")] <- "#e8850c"
  cols[low %in% c("male", "males", "m")] <- "#1f5fbf"
  free <- setdiff(PALETTE, cols)
  cols[is.na(cols)] <- rep(free, length.out = sum(is.na(cols)))
  cols
}

# ------------------------------------------------------------------ reading
read_any <- function(path) {
  first <- readLines(path, n = 1, warn = FALSE)
  seps <- c(";" = ";", "," = ",", "\t" = "\t")
  counts <- sapply(seps, function(s) lengths(regmatches(first, gregexpr(s, first, fixed = TRUE))))
  sep <- if (max(counts) > 0) seps[[which.max(counts)]] else ""
  read.table(path, header = TRUE, sep = sep, quote = "\"'", dec = ".", stringsAsFactors = FALSE,
             check.names = FALSE, comment.char = "", strip.white = TRUE, na.strings = c("NA", ""))
}
guess <- function(cols, candidates) { hit <- candidates[candidates %in% cols]; if (length(hit)) hit[1] else "" }
num <- function(x) if (is.numeric(x)) x else suppressWarnings(as.numeric(sub(",", ".", trimws(as.character(x)), fixed = TRUE)))
chr <- function(x) { x <- trimws(as.character(x)); x[is.na(x) | x == ""] <- "NA"; x }
default_order <- function(values) {
  lv <- sort(unique(chr(values)))
  first <- grepl("^(wild|wt([^a-z]|$)|control|ctrl)", tolower(lv))
  c(lv[first], lv[!first])
}
order_levels <- function(values, text) {
  all_lv <- sort(unique(chr(values)))
  given <- trimws(strsplit(text %||% "", ",")[[1]])
  c(given[given %in% all_lv], setdiff(all_lv, given))
}
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 1e-16, "< 1e-16", ifelse(p < 0.001, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3))))
pstat <- function(p) ifelse(p < 1e-16, "p < 1e-16", paste("p =", fmt_p(p)))
PARAM_TEXT <- c(growth = "growth", lag = "lag", onset = "onset of control", control = "time to control", max = "maximum load")
f1 <- function(v, d = 2) ifelse(is.finite(v), formatC(v, format = "f", digits = d), "")

# ------------------------------------------------------------------ downloads
# Files are made in R and handed to the browser as a data URL, which works under shinylive.
dl_button <- function(id, label) actionButton(id, label, icon = icon("download"), class = "btn-outline-secondary btn-sm me-1")
DL_JS <- "
Shiny.addCustomMessageHandler('send-file', function(m) {
  var a = document.createElement('a');
  a.href = 'data:' + m.mime + ';base64,' + m.data;
  a.download = m.name;
  document.body.appendChild(a); a.click(); a.remove();
});"
send_file <- function(session, name, mime, write_fun) {
  f <- tempfile(fileext = paste0(".", tools::file_ext(name)))
  write_fun(f)
  data <- jsonlite::base64_enc(readBin(f, "raw", file.info(f)$size))
  session$sendCustomMessage("send-file", list(name = name, mime = mime, data = gsub("\n", "", data)))
}

# ------------------------------------------------------------------ reports
params_table <- function(tab, lv) {
  p <- tab$params
  out <- data.frame(Parameter = p$Parameter,
                    a = fmt_ci(p$est1, p$low1, p$high1, 2), b = fmt_ci(p$est2, p$low2, p$high2, 2),
                    d = fmt_ci(p$diff, p$diff_low, p$diff_high, 2),
                    i90 = ifelse(is.finite(p$diff_low90), paste0("[", f1(p$diff_low90), ", ", f1(p$diff_high90), "]"), ""),
                    m = ifelse(is.finite(p$margin), paste0("±", p$margin), ""),
                    dr = sprintf("%.1f", p$draws_both), v = p$verdict, check.names = FALSE)
  names(out) <- c("Parameter", paste(lv[1], "[95% interval]"), paste(lv[2], "[95% interval]"),
                  paste(lv[2], "-", lv[1], "[95% interval]"), "90% interval of the difference", "Margin",
                  "Draws with the parameter in both groups (%)", "Reading")
  out
}

# The methods paragraph, as short paragraphs: a character vector, one element per paragraph
methods_text <- function(all, tabs, s, n_left, n_t0) {
  lv <- s$lv
  first <- Filter(function(r) !inherits(r, "error"), all)[[1]]
  mg <- s$margins[is.finite(s$margins)]
  intro <- c(
    paste0("The bacterial load of each fly was estimated from its colony counts across the dilution series with CFU Drop Counter, and analysed in log2.",
           if (nzchar(s$filter_col) && length(s$filter_vals)) sprintf(" Flies with %s = %s were analysed.", s$filter_col, paste(s$filter_vals, collapse = ", ")) else "",
           if (n_left > 0) sprintf(" %d flies without a load (no colony, or no countable drop) were left out.", n_left) else ""),
    sprintf("The inoculum was measured in %d tubes injected into PBS. It entered the model once, as observations at time 0, and is the common start of both curves.", n_t0),
    sprintf("Log2 load was modelled against time post-injection with a generalized additive model (mgcv): a smooth for %s and, for %s, a difference smooth entered through an ordered factor, each with basis dimension k = %d and penalized by REML. The difference smooth was constrained to zero at time 0.",
            lv[1], lv[2], first$k_use),
    sprintf("Intervals came from %s coefficient vectors drawn from the fitted model, with the covariance including the uncertainty in the smoothing (random seed %d). Whether the groups differ was assessed by the difference smooth and by a simultaneous 95%% band of the difference curve.",
            format(N_SIM, big.mark = ","), s$seed),
    sprintf("Five parameters were read from each curve: growth, the steepest slope (doublings per hour); lag, the time at which the tangent at the steepest point crosses the starting load; onset of control, the time of the steepest slope; time to control, the first time after it at which the slope falls below %.0f%% of its maximum; and the maximum load within the window. Each difference is given with its 95%% interval%s.",
            100 * CONTROL_FRAC,
            if (length(mg)) paste0("; groups were called equivalent where the 90% interval of the difference fell within margins set before the data were examined (", paste(sprintf("%s ± %s", PARAM_TEXT[names(mg)], mg), collapse = ", "), ")") else ""),
    paste0(if (s$replicate != "") "Each replicate experiment was fitted separately. " else "",
           sprintf("Time points used: %s h. Analysis made with the Fly Proliferation app version %s.", paste(s$times_used, collapse = ", "), APP_VERSION)))
  res_par <- unlist(lapply(names(all), function(nm) {
    res <- all[[nm]]; tab <- tabs[[nm]]
    head <- if (nm != "All") paste0(nm, ". ") else ""
    if (inherits(res, "error")) return(sprintf("%sThe model could not be fitted (%s).", head, conditionMessage(res)))
    st <- res$smooth_test; p <- tab$params; kc <- res$kcheck; be <- res$band_excl
    c(paste0(head, sprintf("Results. The difference smooth (%s - %s) had edf = %.2f, F = %.1f, %s. The basis check gave k-index = %.2f (%s). ",
                           lv[2], lv[1], st$edf[2], st$F[2], pstat(st$p[2]), kc$k_index[2], pstat(kc$p[2])),
             "The simultaneous 95% band of the difference ",
             if (nrow(be)) paste0("excluded zero from ", paste(sprintf("%.1f to %.1f h (%s)", be$From, be$To, be$Sign), collapse = " and "), ".") else "included zero over the whole window."),
      paste(c("Parameters, estimate [95% interval]:",
              sprintf("- %s: %s %s; %s %s; difference %s (%s)", PARAM_TEXT[p$key], lv[1], fmt_ci(p$est1, p$low1, p$high1, 2),
                      lv[2], fmt_ci(p$est2, p$low2, p$high2, 2), fmt_ci(p$diff, p$diff_low, p$diff_high, 2), p$verdict)), collapse = "\n"),
      paste0(if (any(!tab$peak)) paste0("For ", paste(lv[!tab$peak], collapse = " and "), ", the maximum is the load at the end of the window, not a peak. ") else "",
             if (any(p$draws_both < 99.5)) paste0("Some parameters were not reached in every draw (", paste(sprintf("%s %.1f%%", PARAM_TEXT[p$key[p$draws_both < 99.5]], p$draws_both[p$draws_both < 99.5]), collapse = ", "),
                                                  " of draws in both groups); their intervals describe the draws in which they were reached.") else ""))
  }))
  out <- c(intro, res_par)
  out[nzchar(out)]
}

r_code <- function(s) {
  core <- paste(readLines("prolif_core.R", warn = FALSE), collapse = "\n")
  mg <- paste(sprintf("%s = %s", names(s$margins), ifelse(is.finite(s$margins), s$margins, "NA")), collapse = ", ")
  cols <- group_colours(s$lv)
  paste0(core, '

# ================================================================== analysis
# Runs on the table downloaded from the Flies tab of the app ("Proliferation_flies.csv")
pf <- read.csv("Proliferation_flies.csv")
d <- subset(pf, used == "yes")
d$group <- factor(d$group, levels = c(', paste0('"', s$lv, '"', collapse = ", "), '))   # first = reference, carries the inoculum
margins <- c(', mg, ')   # set before the data were examined; NA = none
cols <- c(', paste(sprintf('"%s" = "%s"', names(cols), cols), collapse = ", "), ')

for (r in unique(d$replicate)) {                  # each replicate experiment separately
  x <- d[d$replicate == r, ]
  res <- prolif_fit(x, k = ', s$k, ', n_sim = ', N_SIM, ', seed = ', s$seed, ', control_frac = ', CONTROL_FRAC, ')
  tab <- prolif_summary(res, margins)
  cat("\\n==== ", r, "\\n")
  print(res$smooth_test)        # step 2: the difference smooth tests whether the groups differ
  print(res$kcheck)             # step 3: basis check
  print(res$band_excl)          # step 5: where the simultaneous band excludes zero
  print(tab$params)             # steps 6, 8, 9
  print(tab$exist)              # draws in which each quantity exists
  ggsave(paste0("Proliferation_curves_", r, ".pdf"), plot_curves(res, tab, cols), width = 11, height = 8.5)
}
')
}

# ------------------------------------------------------------------ UI
ui <- page_sidebar(
  title = "Early proliferation: two groups compared with a GAM",
  tags$head(tags$script(HTML(DL_JS)), tags$style(HTML("select.shiny-input-select:not(.selectized) { appearance: auto; -webkit-appearance: menulist; }"))),
  sidebar = sidebar(width = 350,
    accordion(id = "steps", open = c("Data", "Comparison"),
      accordion_panel("Data", icon = icon("table"),
        fileInput("file", "1. The flies", accept = c(".csv", ".txt", ".tsv"), placeholder = "Per-fly table from CFU Drop Counter"),
        fileInput("file0", "2. Time zero", accept = c(".csv", ".txt", ".tsv"), placeholder = "Inoculum tubes, from CFU Drop Counter"),
        div(class = "small text-muted", style = "margin-top:-0.6rem; margin-bottom:0.6rem",
            "The per-fly table of the Time 0 plate: the tubes injected into PBS. It is the common start of both curves."),
        actionLink("use_example", "Try the example (Canton-S females and males, P. rettgeri, 1 to 12 h)"),
        div(class = "small text-muted mt-1", textOutput("data_info"))),
      accordion_panel("Comparison", icon = icon("scale-balanced"),
        selectInput("group", "Groups to compare (two)", NULL),
        textInput("order", "Order (comma-separated)", ""),
        div(class = "small text-muted", style = "margin-top:-0.6rem; margin-bottom:0.8rem",
            "The first is the reference: differences are the second minus the first."),
        selectInput("filter_col", "Keep only flies where (optional)", NULL, selectize = FALSE),
        selectizeInput("filter_vals", NULL, choices = NULL, multiple = TRUE, options = list(placeholder = "values to keep")),
        div(class = "small text-muted", style = "margin-top:-0.6rem; margin-bottom:0.8rem",
            "Applies to the flies; the time zero tubes are always kept."),
        selectInput("replicate", "Replicate experiment (optional)", NULL, selectize = FALSE),
        div(class = "small text-muted", style = "margin-top:-0.6rem",
            "Each replicate is fitted separately, and the differences are compared between replicates rather than pooled. If the time zero file has the same column, each replicate starts at its own dose; otherwise all share it.")),
      accordion_panel("Model", icon = icon("chart-line"),
        checkboxGroupInput("times", "Time points used (h)", NULL, inline = TRUE),
        div(class = "small text-muted", style = "margin-top:-0.6rem; margin-bottom:0.8rem",
            "Time 0 is always kept. Leaving out time points shows how much the parameters depend on the sampling density. Use only times before any mortality."),
        numericInput("k", "Basis dimension k", 6, min = 3, max = 20, step = 1),
        div(class = "small text-muted", style = "margin-top:-0.6rem; margin-bottom:0.8rem",
            "Start with 6. Increase it if the basis check in Results shows too little flexibility. It cannot exceed the number of time points minus one."),
        numericInput("seed", "Random seed", 1, min = 1, step = 1),
        p(class = "small text-muted", sprintf("Fixed, as in the protocol: %s draws from the fitted model, and time to control when the growth rate falls below %.0f%% of its maximum.",
                                             format(N_SIM, big.mark = ","), 100 * CONTROL_FRAC))),
      accordion_panel("Equivalence margins", icon = icon("arrows-left-right-to-line"),
        p(class = "small text-muted", "Optional. Set them before looking at the results, and justify them in the paper. A parameter is equivalent when the 90% interval of the difference lies within ± its margin; otherwise the comparison is inconclusive."),
        numericInput("m_growth", "Growth (doublings/h)", NA, min = 0, step = 0.1),
        numericInput("m_lag", "Lag (h)", NA, min = 0, step = 0.5),
        numericInput("m_onset", "Onset of control (h)", NA, min = 0, step = 0.5),
        numericInput("m_control", "Time to control (h)", NA, min = 0, step = 0.5),
        numericInput("m_max", "Maximum (log2 CFU)", NA, min = 0, step = 0.5)),
      accordion_panel("Columns and labels", icon = icon("table-columns"),
        p(class = "small text-muted", "Found automatically in a table from CFU Drop Counter."),
        selectInput("time", "Time post-injection (h)", NULL),
        selectInput("load_col", "Load per fly (log scale)", NULL),
        radioButtons("load_base", NULL, c("log2" = "2", "log10" = "10"), inline = TRUE),
        selectInput("status_col", "Detection status", NULL, selectize = FALSE),
        selectInput("lod_col", "Detection limit per fly (CFU, for the plot of flies)", NULL, selectize = FALSE),
        selectInput("id", "Fly identifier", NULL, selectize = FALSE),
        textInput("xlab", "x-axis label", "Time post-injection (h)"))
    )
  ),
  uiOutput("messages"),
  navset_card_underline(
    nav_panel("Curves",
      uiOutput("rep_pick_curves"),
      div(style = "max-width: 1100px", plotOutput("curves", height = "auto")),
      div(dl_button("dl_pdf", "PDF"), dl_button("dl_png", "PNG"), if (has_svglite) dl_button("dl_svg", "SVG"))),
    nav_panel("Flies over time",
      uiOutput("rep_pick_flies"),
      div(style = "max-width: 1000px", plotOutput("flies_plot", height = "auto")),
      div(dl_button("dl_pdf2", "PDF"), dl_button("dl_png2", "PNG"), if (has_svglite) dl_button("dl_svg2", "SVG"))),
    nav_panel("Results", uiOutput("results_ui"), dl_button("dl_results", "Results (CSV)")),
    nav_panel("Flies", p(class = "small text-muted", "Every fly and time zero tube, with how it was used. The R code runs on this table."),
              tableOutput("pf_table"), dl_button("dl_pf", "Table of flies (CSV)")),
    nav_panel("Methods text", div(style = "max-width: 46rem", uiOutput("methods")), dl_button("dl_methods", "Text")),
    nav_panel("R code", p(class = "small text-muted", "The functions the app runs, then the analysis with the settings chosen. Runs on the table downloaded from the Flies tab."),
              verbatimTextOutput("rcode"), dl_button("dl_code", "R script"))
  )
)

# ------------------------------------------------------------------ server
server <- function(input, output, session) {
  rv <- reactiveValues(d = NULL, name = NULL, d0 = NULL, name0 = NULL)
  observeEvent(input$file, {
    d <- tryCatch(read_any(input$file$datapath), error = function(e) NULL)
    if (is.null(d)) { showNotification("This file could not be read as a table.", type = "error"); return() }
    rv$d <- d; rv$name <- input$file$name
  })
  observeEvent(input$file0, {
    d <- tryCatch(read_any(input$file0$datapath), error = function(e) NULL)
    if (is.null(d)) { showNotification("The time zero file could not be read as a table.", type = "error"); return() }
    rv$d0 <- d; rv$name0 <- input$file0$name
  })
  observeEvent(input$use_example, {
    rv$d <- read_any(EXAMPLE_FILE); rv$name <- EXAMPLE_FILE
    rv$d0 <- read_any(EXAMPLE_TIME0); rv$name0 <- EXAMPLE_TIME0
    # the seed and the equivalence margins of the original analysis
    updateNumericInput(session, "seed", value = 99)
    updateNumericInput(session, "m_growth", value = 0.5); updateNumericInput(session, "m_lag", value = 1)
    updateNumericInput(session, "m_onset", value = 1); updateNumericInput(session, "m_control", value = 1)
    updateNumericInput(session, "m_max", value = 1)
  })

  observeEvent(rv$d, {
    cols <- names(rv$d); none <- c("(none)" = "")
    group_cols <- cols[vapply(rv$d, function(v) { n <- length(unique(v)); n >= 2 && n <= 30 }, TRUE)]
    load <- guess(cols, c("log2_CFU_per_fly", "log10_CFU_per_fly", "log2_load", "log10_load", "log2_CFU_fly", "log10_CFU_fly"))
    updateSelectInput(session, "time", choices = cols, selected = guess(cols, c("Time_post_injection", "Time.post.injection", "Time_point", "Time", "time")))
    updateSelectInput(session, "load_col", choices = cols, selected = load)
    updateRadioButtons(session, "load_base", selected = if (grepl("log10", load)) "10" else "2")
    updateSelectInput(session, "status_col", choices = c(none, cols), selected = guess(cols, c("Status", "status")))
    updateSelectInput(session, "lod_col", choices = c(none, cols), selected = guess(cols, c("LOD_per_fly", "LOD")))
    updateSelectInput(session, "id", choices = c("(none: row number)" = "", cols), selected = guess(cols, c("ID", "Fly_ID", "Sample_ID", "sampleID")))
    updateSelectInput(session, "group", choices = group_cols, selected = guess(group_cols, c("Sex", "Genotype", "Line", "Treatment")))
    updateSelectInput(session, "filter_col", choices = c(none, group_cols), selected = "")
    updateSelectInput(session, "replicate", choices = c(none, group_cols), selected = "")
  })
  observeEvent(input$filter_col, {
    vals <- if (nzchar(input$filter_col)) sort(unique(chr(rv$d[[input$filter_col]]))) else character(0)
    updateSelectizeInput(session, "filter_vals", choices = vals, selected = vals)
  }, ignoreInit = TRUE)

  times_all <- reactive({ req(rv$d, nzchar(input$time %||% "")); tt <- num(rv$d[[input$time]]); sort(unique(tt[is.finite(tt) & tt > 0])) })
  observeEvent(times_all(), updateCheckboxGroupInput(session, "times", choices = times_all(), selected = times_all(), inline = TRUE))
  fly_groups <- reactive({ req(rv$d, nzchar(input$group %||% ""), nzchar(input$time %||% "")); tt <- num(rv$d[[input$time]]); chr(rv$d[[input$group]])[is.finite(tt) & tt > 0] })
  observeEvent(list(rv$d, input$group), { req(nzchar(input$group %||% "")); updateTextInput(session, "order", value = paste(default_order(fly_groups()), collapse = ", ")) })

  settings <- reactive({
    lv_all <- order_levels(fly_groups(), input$order)
    mg <- c(growth = input$m_growth, lag = input$m_lag, onset = input$m_onset, control = input$m_control, max = input$m_max)
    mg <- vapply(mg, function(v) if (is.null(v) || is.na(v) || v <= 0) NA_real_ else as.numeric(v), 0)
    list(group = input$group, lv_all = lv_all, time = input$time, load_col = input$load_col, load_base = as.numeric(input$load_base),
         status_col = input$status_col %||% "", lod_col = input$lod_col %||% "", id = input$id %||% "",
         filter_col = input$filter_col %||% "", filter_vals = input$filter_vals %||% character(0),
         replicate = input$replicate %||% "", times = num(input$times), k = max(3, round(input$k %||% 6)), seed = round(input$seed %||% 1),
         margins = mg, xlab = input$xlab)
  })

  # the time zero tubes: same columns as a per-fly table of CFU Drop Counter, found by name
  time0 <- reactive({
    s <- settings(); d0 <- rv$d0
    validate(need(!is.null(d0), "Load the time zero file too: the per-fly table of the Time 0 plate (tubes injected into PBS), saved by CFU Drop Counter. Both curves start at the dose it measures."))
    cols <- names(d0)
    lc <- if (s$load_col %in% cols) s$load_col else guess(cols, c("log2_CFU_per_fly", "log10_CFU_per_fly", "log2_load", "log10_load", "log2_CFU_fly", "log10_CFU_fly"))
    validate(need(nzchar(lc), sprintf("The time zero file has no load column (%s).", paste(cols, collapse = ", "))))
    base <- if (lc == s$load_col) s$load_base else if (grepl("log10", lc)) 10 else 2
    y <- num(d0[[lc]]) * log2(base)
    st <- if ("Status" %in% cols) as.character(d0$Status) else ifelse(is.finite(y), "Detected", "No load")
    data.frame(ID = if ("ID" %in% cols) as.character(d0$ID) else paste0("time 0 row ", seq_len(nrow(d0))),
               role = "inoculum", group_given = "Time 0", time = 0,
               replicate = if (nzchar(s$replicate) && s$replicate %in% cols) chr(d0[[s$replicate]]) else NA_character_,
               log2_load = round(y, 4), log2_lod = if ("LOD_per_fly" %in% cols) log2(num(d0$LOD_per_fly)) else NA_real_,
               Status = st, stringsAsFactors = FALSE)
  })

  # every fly and time zero tube, with its role and whether it enters the GAM
  flies <- reactive({
    s <- settings(); d <- rv$d; req(d, nzchar(s$time), nzchar(s$load_col), nzchar(s$group))
    tt <- num(d[[s$time]])
    validate(need(all(is.finite(tt)), sprintf("%d row(s) have a time that is not a number.", sum(!is.finite(tt)))),
             need(s$group != s$time, "The groups to compare cannot be the time column."),
             need(!(s$replicate %in% c(s$group, s$time)) || s$replicate == "", "The replicate must be a different column from the groups and the time."))
    y <- num(d[[s$load_col]]) * log2(s$load_base)                     # to log2
    status <- if (nzchar(s$status_col)) as.character(d[[s$status_col]]) else ifelse(is.finite(y), "Detected", "No load")
    lod <- if (nzchar(s$lod_col)) log2(num(d[[s$lod_col]])) else NA_real_
    keep <- tt > 0                                                     # time zero comes from its own file
    if (nzchar(s$filter_col) && length(s$filter_vals)) keep <- keep & chr(d[[s$filter_col]]) %in% s$filter_vals
    pf <- data.frame(ID = if (nzchar(s$id)) as.character(d[[s$id]]) else paste0("row ", seq_len(nrow(d))),
                     role = "fly", group_given = chr(d[[s$group]]), time = tt,
                     replicate = if (nzchar(s$replicate)) chr(d[[s$replicate]]) else "All",
                     log2_load = round(y, 4), log2_lod = lod, Status = status, stringsAsFactors = FALSE)[keep, ]
    lv <- intersect(s$lv_all, unique(pf$group_given))
    validate(need(length(lv) == 2, sprintf("The analysis compares two groups; %s holds %d among the flies kept (%s). Choose two with \"Keep only flies where\", or another column.",
                                           s$group, length(lv), paste(lv, collapse = ", "))))
    # time zero: each replicate its own tubes when the file says which; otherwise shared by all
    z <- time0()
    reps <- unique(pf$replicate)
    z <- if (all(is.na(z$replicate))) do.call(rbind, lapply(reps, function(r) transform(z, replicate = r))) else z[z$replicate %in% reps, ]
    pf <- rbind(z, pf)
    pf$group <- ifelse(pf$role == "inoculum", lv[1], pf$group_given)  # the dose belongs to neither group: bookkeeping
    pf$detected <- pf$Status == "Detected" & is.finite(pf$log2_load)
    pf$used <- ifelse(!pf$detected, paste0("no (", pf$Status, ")"),
                      ifelse(pf$role == "fly" & !(pf$time %in% s$times), "no (time point not used)", "yes"))
    validate(need(sum(pf$used == "yes" & pf$role == "inoculum") > 0, "No time zero tube has a load."),
             need(length(unique(pf$time[pf$used == "yes"])) >= 4, "At least three time points after injection are needed."))
    attr(pf, "lv") <- lv
    attr(pf, "t0_in_flies") <- sum(tt == 0)
    pf
  })
  lv <- reactive(attr(flies(), "lv"))
  reps <- reactive(unique(flies()$replicate))

  # steps 2 to 6: one fit per replicate (only on data, k, seed and time points)
  results <- reactive({
    pf <- flies(); s <- settings(); l <- lv()
    used <- pf[pf$used == "yes", ]
    out <- lapply(setNames(reps(), reps()), function(r) {
      x <- used[used$replicate == r, ]
      x$group <- factor(x$group, levels = l)
      if (!any(x$time == 0)) return(simpleError("no inoculum tube at time 0 in this replicate"))
      tryCatch(prolif_fit(x[, c("time", "group", "log2_load")], k = s$k, n_sim = N_SIM, seed = s$seed, control_frac = CONTROL_FRAC),
               error = function(e) e)
    })
    out
  })
  # steps 6, 8 and 9 (margins change no fit)
  tabs <- reactive({ s <- settings(); lapply(results(), function(r) if (inherits(r, "error")) r else prolif_summary(r, s$margins)) })

  output$data_info <- renderText({
    paste(if (is.null(rv$d)) "Flies: not loaded." else sprintf("Flies: %s, %d rows.", rv$name, nrow(rv$d)),
          if (is.null(rv$d0)) "Time zero: not loaded." else sprintf("Time zero: %s, %d tubes.", rv$name0, nrow(rv$d0)))
  })

  output$messages <- renderUI({
    if (is.null(rv$d) || is.null(rv$d0)) return(div(class = "alert alert-info", "Load the two tables saved by CFU Drop Counter (log2): the flies, and the time zero tubes injected into PBS. Or try the example."))
    pf <- tryCatch(flies(), error = function(e) e)
    if (inherits(pf, "error")) return(div(class = "alert alert-danger", conditionMessage(pf)))
    w <- character(0)
    nd <- pf[pf$role == "fly" & !pf$detected & pf$time %in% settings()$times, ]
    if (nrow(nd)) {
      tb <- table(paste0(nd$group, " at ", nd$time, " h"))
      share <- tapply(!pf$detected[pf$role == "fly"], paste(pf$group, pf$time)[pf$role == "fly"], mean)
      w <- c(w, paste0("Flies without a load are left out of the GAM (", paste(sprintf("%s: %d", names(tb), tb), collapse = ", "), "), so the curve is fitted to the flies that have one",
                       if (max(share) > 0.1) "; at more than one fly in ten, this lifts the curve at those times. Report them, and consider analysing them as a binary outcome." else ". Report them."))
    }
    if (attr(pf, "t0_in_flies") > 0) w <- c(w, sprintf("%d row(s) of the flies table are at time 0 and are not used: time zero comes from its own file.", attr(pf, "t0_in_flies")))
    inoc <- pf[pf$role == "inoculum" & !duplicated(pf$ID), ]
    if (any(!inoc$detected)) w <- c(w, sprintf("%d time zero tube(s) have no load and are left out.", sum(!inoc$detected)))
    for (nm in names(results())) {
      r <- results()[[nm]]; lab <- if (nm == "All") "" else paste0(nm, ": ")
      if (inherits(r, "error")) { w <- c(w, paste0(lab, if (nzchar(lab)) "the" else "The", " model could not be fitted (", conditionMessage(r), ").")); next }
      kc <- r$kcheck
      if (any(kc$p < 0.05 & kc$edf > 0.8 * kc$k)) w <- c(w, paste0(lab, if (nzchar(lab)) "the" else "The", " basis check suggests too little flexibility (k-index below 1 with a small p and edf close to k'). Increase k and refit."))
      if (r$k_use < r$k) w <- c(w, sprintf("%sk was lowered to %d, the number of time points minus one.", lab, r$k_use))
      tb <- tabs()[[nm]]
      if (any(tb$params$draws_both < 95)) w <- c(w, paste0(lab, if (nzchar(lab)) "some" else "Some", " parameters are missing in more than 5% of the draws (see Results); their intervals describe only the draws that reached them, the ones that reached them soonest."))
    }
    if (length(reps()) > 1) w <- c(w, "Each replicate is fitted separately. Report the differences of every replicate and their agreement (Results), rather than counting how many are significant.")
    if (length(w)) div(class = "alert alert-warning", tags$ul(lapply(w, tags$li)))
  })

  rep_pick <- function(id) renderUI(if (length(reps()) > 1) selectInput(id, "Replicate shown", reps(), selected = reps()[1]))
  output$rep_pick_curves <- rep_pick("rep_curves")
  output$rep_pick_flies <- rep_pick("rep_flies")
  pick <- function(v) if (!is.null(v) && v %in% reps()) v else reps()[1]

  the_curves <- reactive({
    r <- pick(input$rep_curves); res <- results()[[r]]; tab <- tabs()[[r]]
    validate(need(!inherits(res, "error"), "The model could not be fitted for this replicate."))
    p <- plot_curves(res, tab, group_colours(lv()), xlab = settings()$xlab)
    if (r != "All") p + patchwork::plot_annotation(title = r) else p
  })
  # same shape on screen as the saved file (11 x 8.5 in), whatever the width of the window
  output$curves <- renderPlot(the_curves(), res = 96, height = function() round(session$clientData$output_curves_width * 8.5 / 11))

  the_flies_plot <- reactive({
    r <- pick(input$rep_flies); pf <- flies(); s <- settings()
    fl <- pf[pf$replicate == r & (pf$role == "inoculum" | pf$time %in% s$times) & !grepl("^no \\((Above|Estimation)", pf$used), ]
    p <- plot_flies(fl, lv(), 0, group_colours(lv()), xlab = s$xlab)
    if (r != "All") p + ggtitle(r) else p
  })
  output$flies_plot <- renderPlot(the_flies_plot(), res = 96, height = function() round(session$clientData$output_flies_plot_width * 5 / 8))

  output$results_ui <- renderUI({
    l <- lv()
    blocks <- lapply(names(results()), function(nm) {
      res <- results()[[nm]]; tab <- tabs()[[nm]]
      if (inherits(res, "error")) return(p(sprintf("%s: the model could not be fitted (%s).", nm, conditionMessage(res))))
      st <- res$smooth_test; kc <- res$kcheck; be <- res$band_excl
      tagList(if (nm != "All") h5(nm),
        h6("Does the difference smooth differ from zero? (step 2)"),
        renderTable(transform(st, edf = f1(edf), F = f1(F), p = fmt_p(p))),
        h6(sprintf("Basis check, k = %d (step 3)", res$k_use)),
        renderTable(transform(kc, edf = f1(edf), k_index = f1(k_index), p = fmt_p(p))),
        p(class = "small text-muted", "A k-index below 1 with a small p and an edf close to k' means too little flexibility: increase k."),
        h6("Where the groups differ: simultaneous 95% band of the difference (step 5)"),
        p(if (nrow(be)) paste0("The band excludes zero from ", paste(sprintf("%.1f to %.1f h (%s)", be$From, be$To, be$Sign), collapse = "; "), ".")
          else "The band includes zero over the whole window."),
        p(class = "small text-muted", "At time 0 the difference is zero by construction, so where the band first excludes zero depends on the anchoring and the smoothing; it is not an estimate of when the groups began to diverge."),
        h6(sprintf("Five parameters and their differences (steps 6, 8 and 9; seed %d, %s draws)", res$seed, format(res$n_sim, big.mark = ","))),
        renderTable(params_table(tab, l)),
        p(class = "small text-muted", "Five parameters were examined; weigh the pattern rather than one verdict. Growth, lag and onset of control are all read from the steepest point, so they are not independent. Differ: the 95% interval of the difference excludes 0. Equivalent: its 90% interval lies within the margin. A non-significant difference is not evidence of equality."),
        h6("Draws in which each quantity exists"),
        renderTable(tab$exist, digits = 1),
        if (any(!tab$peak)) p(sprintf("For %s, the maximum is the load at the last time point, not a peak: label it as such.", paste(l[!tab$peak], collapse = " and "))),
        h6("Flies and inoculum tubes in the model"), renderTable(as.data.frame(t(res$n))))
    })
    per_rep <- NULL
    ok <- Filter(function(t) !inherits(t, "error"), tabs())
    if (length(ok) > 1) {
      pr <- do.call(cbind, lapply(ok, function(t) fmt_ci(t$params$diff, t$params$diff_low, t$params$diff_high, 2)))
      per_rep <- tagList(h5(sprintf("Differences per replicate, %s - %s [95%% interval]", l[2], l[1])),
                         renderTable(data.frame(Parameter = PARAM_NAMES, pr, check.names = FALSE)))
    }
    tagList(per_rep, blocks)
  })

  results_csv <- reactive({
    do.call(rbind, lapply(names(tabs()), function(nm) {
      t <- tabs()[[nm]]; if (inherits(t, "error")) return(NULL)
      p <- t$params; p$key <- NULL
      names(p)[2:7] <- c(paste0(lv()[1], c("", "_low95", "_high95")), paste0(lv()[2], c("", "_low95", "_high95")))
      cbind(Replicate = nm, p)
    }))
  })

  pf_out <- reactive({
    pf <- flies()
    data.frame(ID = pf$ID, role = pf$role, group = pf$group, group_given = pf$group_given, time = pf$time,
               replicate = pf$replicate, log2_load = pf$log2_load, Status = pf$Status, used = pf$used)
  })
  output$pf_table <- renderTable(head(pf_out(), 60))

  methods <- reactive({
    s <- settings(); s$lv <- lv(); s$times_used <- sort(c(0, s$times)); pf <- flies()
    methods_text(results(), tabs(), s, sum(pf$role == "fly" & !pf$detected),
                 sum(pf$role == "inoculum" & pf$used == "yes" & !duplicated(pf$ID)))
  })
  output$methods <- renderUI(lapply(methods(), function(x) p(style = "white-space: pre-line", x)))
  code <- reactive({ s <- settings(); s$lv <- lv(); r_code(s) })
  output$rcode <- renderText(code())

  observeEvent(input$dl_pdf, send_file(session, "Proliferation_curves.pdf", "application/pdf", function(f) ggsave(f, the_curves(), width = 11, height = 8.5)))
  observeEvent(input$dl_png, send_file(session, "Proliferation_curves.png", "image/png", function(f) ggsave(f, the_curves(), width = 11, height = 8.5, dpi = 300)))
  observeEvent(input$dl_svg, send_file(session, "Proliferation_curves.svg", "image/svg+xml", function(f) ggsave(f, the_curves(), width = 11, height = 8.5, device = svglite::svglite)))
  observeEvent(input$dl_pdf2, send_file(session, "Proliferation_flies.pdf", "application/pdf", function(f) ggsave(f, the_flies_plot(), width = 8, height = 5)))
  observeEvent(input$dl_png2, send_file(session, "Proliferation_flies.png", "image/png", function(f) ggsave(f, the_flies_plot(), width = 8, height = 5, dpi = 300)))
  observeEvent(input$dl_svg2, send_file(session, "Proliferation_flies.svg", "image/svg+xml", function(f) ggsave(f, the_flies_plot(), width = 8, height = 5, device = svglite::svglite)))
  observeEvent(input$dl_results, send_file(session, "Proliferation_results.csv", "text/csv", function(f) write.table(results_csv(), f, sep = ";", row.names = FALSE)))
  observeEvent(input$dl_pf, send_file(session, "Proliferation_flies.csv", "text/csv", function(f) write.csv(pf_out(), f, row.names = FALSE)))
  observeEvent(input$dl_methods, send_file(session, "Proliferation_methods.txt", "text/plain", function(f) writeLines(paste(methods(), collapse = "\n\n"), f)))
  observeEvent(input$dl_code, send_file(session, "Proliferation_analysis.R", "text/plain", function(f) writeLines(code(), f)))
}

shinyApp(ui, server)

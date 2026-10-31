# Early bacterial proliferation: comparison of two groups with a GAM
# As in the protocol, "Quantification and statistical analysis", section "Early proliferation".
# The Fly Proliferation app runs these functions, and its R code tab gives this file with the
# settings chosen, so the analysis can be rerun in R without the app.
#
# Input: one row per fly (and per inoculum tube), with
#   time       hours post-injection; the inoculum, measured by injecting into PBS, at time 0
#   group      factor with two levels; the first is the reference. The inoculum rows carry the
#              reference level: the dose was measured in PBS and belongs to neither group, the
#              level is bookkeeping (protocol, step 2)
#   log2_load  log2 bacteria per fly (CFU Drop Counter, Load step), flies with a load only

library(mgcv)
library(ggplot2)
library(patchwork)

PARAM_NAMES <- c(growth = "Growth (doublings/h)", lag = "Lag (h)", onset = "Onset of control (h)",
                 control = "Time to control (h)", max = "Maximum (log2 CFU)")

# ------------------------------------------------------------------ steps 2 to 6
# Fit the GAM, draw coefficient vectors, and compute the curves, the simultaneous band of the
# difference and the five parameters of each curve, for the fit and for every draw.
prolif_fit <- function(d, k = 6, n_sim = 10000, seed = 1, control_frac = 0.1, n_grid = 200) {
  lv <- levels(d$group)
  stopifnot("two groups are needed" = length(lv) == 2)
  t0 <- min(d$time)
  if (any(d$time == t0 & d$group != lv[1]))
    stop("rows at time ", t0, " must carry the reference level (", lv[1], ")")
  # k.check and the draws are random: set the seed, and give the caller back its own random state
  old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old_seed, envir = globalenv()))
  set.seed(seed)
  d$groupO <- ordered(d$group, levels = lv)
  k_use <- min(k, length(unique(d$time)) - 1)

  # Step 2: one smooth for the reference, a difference smooth for the other group through an
  # ordered factor, each penalized by REML; the difference smooth is constrained to 0 at t0.
  # Its p-value (summary) is approximate: in simulations of the protocol's design with no
  # difference, it was below 0.05 in about 11% of experiments. The simultaneous band (step 5)
  # held its level (about 6%), so it is the one to read for whether the groups differ.
  fit <- gam(log2_load ~ s(time, k = k_use) + s(time, by = groupO, k = k_use, pc = t0),
             data = d, method = "REML")
  st <- summary(fit)$s.table
  smooth_test <- data.frame(Term = c(paste("Smooth of", lv[1]), paste("Difference smooth,", lv[2], "-", lv[1])),
                            edf = st[, "edf"], F = st[, "F"], p = st[, "p-value"], row.names = NULL)
  # Step 3: basis check
  kc <- k.check(fit)
  kcheck <- data.frame(Term = smooth_test$Term, k = kc[, "k'"], edf = kc[, "edf"],
                       k_index = kc[, "k-index"], p = kc[, "p-value"], row.names = NULL)

  # Step 4: draws, with the covariance that includes the uncertainty in the smoothing
  V <- vcov(fit, unconditional = TRUE)
  b <- coef(fit)
  draws <- rmvn(n_sim, b, V)

  tgrid <- seq(t0, max(d$time), length.out = n_grid)
  newd <- expand.grid(time = tgrid, group = factor(lv, levels = lv))
  newd$groupO <- ordered(newd$group, levels = lv)
  h <- 0.01                                                           # slope by finite difference
  X <- predict(fit, newd, type = "lpmatrix")
  D <- (predict(fit, transform(newd, time = time + h), type = "lpmatrix") - X) / h
  r1 <- newd$group == lv[1]; r2 <- newd$group == lv[2]

  band <- function(M) {                                               # pointwise 95%
    est <- drop(M %*% b); se <- sqrt(pmax(rowSums((M %*% V) * M), 0))
    data.frame(est = est, low = est - 1.96 * se, high = est + 1.96 * se)
  }
  # Step 5: simultaneous 95% band of the difference, from the draws
  Mdiff <- X[r2, , drop = FALSE] - X[r1, , drop = FALSE]
  dif <- band(Mdiff)
  se_d <- (dif$high - dif$est) / 1.96
  ok <- se_d > 1e-8                                                   # se = 0 at t0, by construction
  dev <- abs((draws - rep(b, each = n_sim)) %*% t(Mdiff[ok, , drop = FALSE])) / rep(se_d[ok], each = n_sim)
  crit <- unname(quantile(dev[cbind(seq_len(n_sim), max.col(dev, ties.method = "first"))], 0.95))
  dif$low <- dif$est - crit * se_d; dif$high <- dif$est + crit * se_d
  excl <- dif$low > 0 | dif$high < 0
  rr <- rle(excl); ends <- cumsum(rr$lengths); starts <- ends - rr$lengths + 1
  band_excl <- data.frame(From = tgrid[starts][rr$values], To = tgrid[ends][rr$values],
                          Sign = ifelse(dif$est[starts][rr$values] > 0, paste(lv[2], "above"), paste(lv[2], "below")))

  curves <- rbind(
    data.frame(panel = 1, group = rep(lv, each = n_grid), time = tgrid, rbind(band(X[r1, ]), band(X[r2, ]))),
    data.frame(panel = 2, group = rep(lv, each = n_grid), time = tgrid, rbind(band(D[r1, ]), band(D[r2, ]))),
    data.frame(panel = 3, group = "Difference", time = tgrid, dif))

  # Step 6: the five parameters, for the fitted curves and for every draw
  par_of <- function(Fm, Dm) {                    # Fm, Dm: one row per curve, one column per grid time
    n <- nrow(Fm); idx <- seq_len(n)
    i <- max.col(Dm, ties.method = "first")       # t*: steepest point
    growth <- Dm[cbind(idx, i)]
    fi <- Fm[cbind(idx, i)]; f0 <- Fm[, 1]
    below <- (col(Dm) > i) & (Dm <= control_frac * growth)
    has <- rowSums(below) > 0
    j <- max.col(below * 1, ties.method = "first")
    im <- max.col(Fm, ties.method = "first")
    data.frame(growth = growth,
               lag = ifelse(growth > 0, tgrid[i] - (fi - f0) / growth, NA),
               onset = tgrid[i],
               control = ifelse(has, tgrid[j], NA),          # NA: not reached in the window
               max = Fm[cbind(idx, im)],
               t_max = tgrid[im], y_growth = fi, y_start = f0,
               y_control = ifelse(has, Dm[cbind(idx, j)], NA))
  }
  est <- sim <- list()
  for (g in lv) {
    Xg <- X[newd$group == g, ]; Dg <- D[newd$group == g, ]
    est[[g]] <- par_of(rbind(drop(Xg %*% b)), rbind(drop(Dg %*% b)))
    sim[[g]] <- par_of(draws %*% t(Xg), draws %*% t(Dg))
  }
  list(fit = fit, lv = lv, t0 = t0, k = k, k_use = k_use, seed = seed, n_sim = n_sim,
       control_frac = control_frac, tgrid = tgrid, curves = curves, band_crit = crit,
       band_excl = band_excl, smooth_test = smooth_test, kcheck = kcheck, est = est, sim = sim,
       data = d[, c("time", "group", "log2_load")],
       n = c(inoculum = sum(d$time == t0), table(droplevels(d$group[d$time > t0]))))
}

# ------------------------------------------------------------------ steps 6, 8 and 9
# Each parameter with its 95% interval, the difference with its 95% interval, and, where a
# margin was set before the data were examined, equivalence from the 90% interval.
prolif_summary <- function(res, margins = c(growth = NA, lag = NA, onset = NA, control = NA, max = NA)) {
  lv <- res$lv; e1 <- res$est[[lv[1]]]; e2 <- res$est[[lv[2]]]; s1 <- res$sim[[lv[1]]]; s2 <- res$sim[[lv[2]]]
  end <- max(res$tgrid)
  qnt <- function(x, p) if (sum(is.finite(x)) > 1) unname(quantile(x, p, na.rm = TRUE)) else c(NA, NA)
  rows <- lapply(names(PARAM_NAMES), function(p) {
    x <- s2[[p]] - s1[[p]]; both <- is.finite(x)
    d_est <- e2[[p]] - e1[[p]]
    ci <- qnt(x, c(0.025, 0.975)); ci90 <- qnt(x, c(0.05, 0.95))
    m <- if (p %in% names(margins)) as.numeric(margins[[p]]) else NA_real_
    differ <- all(is.finite(ci)) && (ci[1] > 0 || ci[2] < 0)
    within <- is.finite(m) && all(is.finite(ci90)) && ci90[1] > -m && ci90[2] < m
    verdict <- if (!is.finite(d_est)) "not reached" else
      if (differ && within) "differ, within margin" else      # real, but smaller than the margin
      if (differ) "differ" else
      if (!is.finite(m)) "no margin set" else
      if (within) "equivalent" else "inconclusive"
    data.frame(Parameter = PARAM_NAMES[[p]], key = p,
               est1 = e1[[p]], low1 = qnt(s1[[p]], 0.025)[1], high1 = qnt(s1[[p]], 0.975)[1],
               est2 = e2[[p]], low2 = qnt(s2[[p]], 0.025)[1], high2 = qnt(s2[[p]], 0.975)[1],
               diff = d_est, diff_low = ci[1], diff_high = ci[2], diff_low90 = ci90[1], diff_high90 = ci90[2],
               margin = m, draws_both = 100 * mean(both), verdict = verdict, row.names = NULL)
  })
  params <- do.call(rbind, rows)
  # the maximum is a peak only when it is reached before the end of the window
  peak <- c(e1$t_max < end, e2$t_max < end)
  exist <- data.frame(Quantity = c("Lag defined (steepest slope > 0)", "Time to control reached", "Maximum is a peak, not the end of the window"),
                      g1 = 100 * c(mean(is.finite(s1$lag)), mean(is.finite(s1$control)), mean(s1$t_max < end)),
                      g2 = 100 * c(mean(is.finite(s2$lag)), mean(is.finite(s2$control)), mean(s2$t_max < end)))
  names(exist)[2:3] <- paste0(lv, " (% of draws)")
  list(params = params, exist = exist, peak = setNames(peak, lv))
}

# ------------------------------------------------------------------ figure
# The protocol's figure: fitted log2 load, growth rate, difference, and the five parameters.
fmt_ci <- function(est, low, high, digits = 1) {
  f <- function(v) ifelse(is.finite(v), formatC(v, format = "f", digits = digits), "n.r.")
  ifelse(is.finite(est), paste0(f(est), " [", f(low), ", ", f(high), "]"), "not reached")
}

plot_curves <- function(res, tab, cols, xlab = "Time post-injection (h)") {
  lv <- res$lv; e <- res$est; p <- tab$params
  cv <- res$curves
  pts <- res$data
  pts$group <- ifelse(pts$time == res$t0, "Start", as.character(pts$group))
  ci <- function(g, k) quantile(res$sim[[g]][[k]], c(0.025, 0.975), na.rm = TRUE)
  nudge <- setNames(c(-0.35, 0.35), lv)               # the two lag points share the start level
  marks <- do.call(rbind, lapply(lv, function(g) {
    x <- e[[g]]
    data.frame(group = g, panel = c(1, 1, 2, 2),
               x = c(x$lag, x$t_max, x$onset, x$control),
               y = c(x$y_start + nudge[[g]], x$max, x$growth, x$y_control),
               xlo = c(ci(g, "lag")[1], NA, ci(g, "onset")[1], ci(g, "control")[1]),
               xhi = c(ci(g, "lag")[2], NA, ci(g, "onset")[2], ci(g, "control")[2]),
               ylo = c(NA, ci(g, "max")[1], ci(g, "growth")[1], NA),
               yhi = c(NA, ci(g, "max")[2], ci(g, "growth")[2], NA))
  }))
  marks <- marks[is.finite(marks$x) & is.finite(marks$y), ]
  tangent <- do.call(rbind, lapply(lv, function(g) {
    x <- e[[g]]; data.frame(group = g, x = x$lag, xend = x$onset, y = x$y_start + nudge[[g]], yend = x$y_growth)
  }))
  tangent <- tangent[is.finite(tangent$x), ]
  div <- data.frame(x = res$band_excl$From); div <- div[div$x > res$t0, , drop = FALSE]

  allcols <- c(cols[lv], Start = "grey35", Difference = "black")
  sc <- function(keys, legend) list(
    scale_colour_manual(values = allcols, breaks = keys, limits = names(allcols), name = NULL),
    scale_fill_manual(values = allcols, breaks = keys, limits = names(allcols), name = NULL),
    scale_x_continuous(breaks = sort(unique(res$data$time))),
    theme_classic(base_size = 11),
    theme(plot.title = element_text(face = "bold", size = 11)),
    if (!legend) guides(colour = "none", fill = "none"))
  mk <- function(k) list(
    geom_segment(data = marks[marks$panel == k & !is.na(marks$xlo), ], inherit.aes = FALSE,
                 aes(x = xlo, xend = xhi, y = y, yend = y, colour = group), linewidth = 0.7),
    geom_segment(data = marks[marks$panel == k & !is.na(marks$ylo), ], inherit.aes = FALSE,
                 aes(x = x, xend = x, y = ylo, yend = yhi, colour = group), linewidth = 0.7),
    geom_point(data = marks[marks$panel == k, ], inherit.aes = FALSE, aes(x = x, y = y, colour = group),
               shape = 21, fill = "white", size = 2.6, stroke = 1.1))
  curve <- function(k) ggplot(cv[cv$panel == k, ], aes(time, est, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, colour = NA) + geom_line(linewidth = 0.9)

  p_load <- curve(1) +
    geom_point(data = pts, aes(time, log2_load, colour = group), inherit.aes = FALSE, size = 1, alpha = 0.45) +
    geom_segment(data = tangent, inherit.aes = FALSE, aes(x = x, xend = xend, y = y, yend = yend, colour = group),
                 linetype = "dotted", linewidth = 0.6) +
    mk(1) + sc(c(lv, "Start"), TRUE) + labs(title = "Fitted load", x = xlab, y = "log2 CFU per fly")
  p_rate <- curve(2) + geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
    mk(2) + sc(lv, FALSE) + labs(title = "Growth rate", x = xlab, y = "Doublings per hour")
  p_diff <- curve(3) + geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
    geom_vline(data = div, aes(xintercept = x), linetype = "dotted", colour = "grey30") +
    geom_text(data = div, inherit.aes = FALSE, aes(x = x, y = Inf, label = paste0("differ from ", round(x, 1), " h")),
              hjust = 1.05, vjust = 1.5, size = 3.3) +
    sc("Difference", FALSE) + labs(title = paste0("Difference, ", lv[2], " - ", lv[1]), x = xlab, y = "log2 difference")

  # the five parameters: the difference, its 95% interval and the reading
  short <- c(growth = "Growth (doublings/h)", lag = "Lag (h)", onset = "Onset of control (h)", control = "Time to control (h)", max = "Maximum (log2 CFU)")
  tb <- data.frame(y = rev(seq_len(nrow(p))), par = short[p$key], d = fmt_ci(p$diff, p$diff_low, p$diff_high), v = p$verdict)
  hdr <- nrow(p) + 1.2
  p_tab <- ggplot(tb) +
    geom_text(aes(0, y, label = par), hjust = 0, size = 3.4) +
    geom_text(aes(0.45, y, label = d), hjust = 0, size = 3.4) +
    geom_text(aes(0.82, y, label = v, fontface = ifelse(startsWith(v, "differ"), "bold", "plain")), hjust = 0, size = 3.2) +
    annotate("text", x = c(0, 0.45, 0.82), y = hdr, label = c("Parameter", "Difference\n[95% interval]", "Reading"), vjust = 0.3, lineheight = 0.9,
             hjust = 0, size = 3.4, fontface = "bold") +
    annotate("segment", x = 0, xend = 1.02, y = hdr - 0.65, yend = hdr - 0.65, linewidth = 0.3) +
    scale_x_continuous(limits = c(0, 1.05), expand = c(0, 0)) + scale_y_continuous(limits = c(0.4, hdr + 0.7)) +
    labs(title = paste0("Five parameters examined, ", lv[2], " - ", lv[1])) + theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11), plot.margin = margin(5, 5, 5, 15))

  (p_load + p_rate) / (p_diff + p_tab) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(caption = paste0(
      "GAM fits with 95% bands (simultaneous for the difference), including smoothing uncertainty; both groups start at the dose measured at time 0.\n",
      "Fitted load: end of the lag (dotted tangent at the steepest point) and maximum. Growth rate: fastest growth = onset of control, and time to control\n",
      "(growth below ", 100 * res$control_frac, "% of its maximum). Bars: 95% intervals from ", format(res$n_sim, big.mark = ","), " draws. ",
      "Difference: dotted line from where the simultaneous band excludes 0.\n",
      "Readings: differ = 95% interval excludes 0; equivalent = 90% interval within the margin set beforehand; otherwise inconclusive. n.r.: not reached by ",
      max(res$tgrid), " h."),
      theme = theme(plot.caption = element_text(hjust = 0, size = 8.5))) &
    theme(legend.position = "top")
}

# Individual flies over time, with medians from the common start and the detection limit.
# 'fl': time, group, log2_load (NA when not detected), log2_lod, detected (TRUE/FALSE)
plot_flies <- function(fl, lv, t0, cols, xlab = "Time post-injection (h)") {
  fl$grp <- ifelse(fl$time == t0, "Start", as.character(fl$group))
  off <- c(setNames(c(-0.17, 0.17), lv), Start = 0)
  fl$y <- ifelse(fl$detected, fl$log2_load, fl$log2_lod)
  fl <- fl[is.finite(fl$y), ]
  fl$x <- fl$time + off[fl$grp]
  fl$Status <- ifelse(fl$detected, "Detected", "Not detected (at its detection limit)")
  med <- aggregate(y ~ grp + time, data = fl, FUN = median)
  med$x <- med$time + off[med$grp]
  m0 <- med[med$grp == "Start", ]
  lines <- rbind(med[med$grp != "Start", ], transform(m0, grp = lv[1]), transform(m0, grp = lv[2]))
  lod <- if (any(is.finite(fl$log2_lod))) aggregate(log2_lod ~ time, data = fl, FUN = median) else data.frame(time = numeric(0), log2_lod = numeric(0))
  allcols <- c(cols[lv], Start = "grey35")
  ggplot(fl, aes(x, y, colour = grp)) +
    geom_segment(data = lod, inherit.aes = FALSE, aes(x = time - 0.4, xend = time + 0.4, y = log2_lod, yend = log2_lod),
                 linetype = "dashed", colour = "grey55") +
    geom_point(aes(shape = Status), size = 1.8, alpha = 0.6, position = position_jitter(width = 0.08, height = 0, seed = 1)) +
    geom_line(data = lines, aes(x, y, group = grp), linewidth = 0.9) +
    geom_point(data = med, aes(x, y), size = 3) +
    scale_colour_manual(values = allcols, breaks = names(allcols), name = NULL) +
    scale_shape_manual(values = c("Detected" = 16, "Not detected (at its detection limit)" = 1), name = NULL) +
    scale_x_continuous(breaks = sort(unique(fl$time))) +
    labs(x = xlab, y = expression("Bacterial load (" * log[2] ~ "CFU per fly)"),
         caption = "Points: individual flies, and the injected dose at time 0. Lines: medians from the common start. Dashed: detection limit.") +
    theme_classic(base_size = 13) + theme(legend.position = "top", plot.caption = element_text(hjust = 0, size = 8))
}

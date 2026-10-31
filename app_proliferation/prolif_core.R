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
  set.seed(seed)                     # k.check and the draws are random
  d$groupO <- ordered(d$group, levels = lv)
  k_use <- min(k, length(unique(d$time)) - 1)

  # Step 2: one smooth for the reference, a difference smooth for the other group through an
  # ordered factor, each penalized by REML; the difference smooth is constrained to 0 at t0
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
  q <- function(x, p) if (sum(is.finite(x)) > 1) unname(quantile(x, p, na.rm = TRUE)) else c(NA, NA)
  rows <- lapply(names(PARAM_NAMES), function(p) {
    x <- s2[[p]] - s1[[p]]; both <- is.finite(x)
    d_est <- e2[[p]] - e1[[p]]
    ci <- q(x, c(0.025, 0.975)); ci90 <- q(x, c(0.05, 0.95))
    m <- margins[[p]]
    verdict <- if (!is.finite(d_est)) "not reached" else
      if (all(is.finite(ci)) && (ci[1] > 0 || ci[2] < 0)) "differ" else
      if (!is.finite(m)) "no margin set" else
      if (all(is.finite(ci90)) && ci90[1] > -m && ci90[2] < m) "equivalent" else "inconclusive"
    data.frame(Parameter = PARAM_NAMES[[p]], key = p,
               est1 = e1[[p]], low1 = q(s1[[p]], 0.025)[1], high1 = q(s1[[p]], 0.975)[1],
               est2 = e2[[p]], low2 = q(s2[[p]], 0.025)[1], high2 = q(s2[[p]], 0.975)[1],
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
# The three panels of the protocol's figure: fitted log2 load, growth rate, difference.
fmt_ci <- function(est, low, high, digits = 1) {
  f <- function(v) ifelse(is.finite(v), formatC(v, format = "f", digits = digits), "n.r.")
  ifelse(is.finite(est), paste0(f(est), " [", f(low), ", ", f(high), "]"), "not reached")
}

plot_curves <- function(res, tab, cols, xlab = "Time post-injection (h)") {
  lv <- res$lv; e <- res$est; p <- tab$params
  panels <- c("Fitted log2 load", "Growth rate (doublings per hour)",
              paste0("Difference: ", lv[2], " - ", lv[1], " (log2)"))
  cv <- res$curves; cv$panel <- factor(panels[cv$panel], levels = panels)
  pts <- res$data
  pts$group <- ifelse(pts$time == res$t0, "Start", as.character(pts$group))
  pts$panel <- factor(panels[1], levels = panels)
  ci <- function(g, k) quantile(res$sim[[g]][[k]], c(0.025, 0.975), na.rm = TRUE)
  nudge <- setNames(c(-0.35, 0.35), lv)               # the two lag points share the start level
  marks <- do.call(rbind, lapply(lv, function(g) {
    x <- e[[g]]
    data.frame(group = g, panel = factor(panels[c(1, 1, 2, 2)], levels = panels),
               x = c(x$lag, x$t_max, x$onset, x$control),
               y = c(x$y_start + nudge[[g]], x$max, x$growth, x$y_control),
               xlo = c(ci(g, "lag")[1], NA, ci(g, "onset")[1], ci(g, "control")[1]),
               xhi = c(ci(g, "lag")[2], NA, ci(g, "onset")[2], ci(g, "control")[2]),
               ylo = c(NA, ci(g, "max")[1], ci(g, "growth")[1], NA),
               yhi = c(NA, ci(g, "max")[2], ci(g, "growth")[2], NA))
  }))
  marks <- marks[is.finite(marks$x) & is.finite(marks$y), ]
  tangent <- do.call(rbind, lapply(lv, function(g) {
    x <- e[[g]]
    data.frame(group = g, panel = factor(panels[1], levels = panels),
               x = x$lag, xend = x$onset, y = x$y_start + nudge[[g]], yend = x$y_growth)
  }))
  tangent <- tangent[is.finite(tangent$x), ]
  dd <- cv[cv$group == "Difference", ]
  div <- data.frame(group = "Difference", panel = factor(panels[3], levels = panels), x = res$band_excl$From)
  div <- div[div$x > res$t0, , drop = FALSE]
  zero <- data.frame(panel = factor(panels[2:3], levels = panels), y = 0)

  # above the panels: the difference for each parameter, its 95% interval and the reading
  short <- c(growth = "Growth (dbl/h)", lag = "Lag (h)", onset = "Ctrl onset (h)", control = "Control (h)", max = "Max (log2)")
  w <- 19
  ab <- if (substr(lv[1], 1, 1) != substr(lv[2], 1, 1)) substr(lv, 1, 1) else abbreviate(lv, 3)
  row_txt <- function(name, vals) paste0(formatC(name, width = -10), paste(formatC(vals, width = -w), collapse = ""))
  subtitle <- paste(row_txt("", short[p$key]),
                    row_txt(paste0(ab[2], "-", ab[1]), fmt_ci(p$diff, p$diff_low, p$diff_high)),
                    row_txt("", p$verdict), sep = "\n")
  allcols <- c(cols[lv], Start = "grey35", Difference = "black")
  ggplot(cv, aes(time, est, colour = group, fill = group)) +
    geom_hline(data = zero, aes(yintercept = y), colour = "grey60", linetype = "dashed") +
    geom_point(data = pts, aes(time, log2_load, colour = group), inherit.aes = FALSE, size = 1, alpha = 0.45) +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 0.9) +
    geom_segment(data = tangent, inherit.aes = FALSE, aes(x = x, xend = xend, y = y, yend = yend, colour = group),
                 linetype = "dotted", linewidth = 0.6) +
    geom_segment(data = marks[!is.na(marks$xlo), ], inherit.aes = FALSE,
                 aes(x = xlo, xend = xhi, y = y, yend = y, colour = group), linewidth = 0.7) +
    geom_segment(data = marks[!is.na(marks$ylo), ], inherit.aes = FALSE,
                 aes(x = x, xend = x, y = ylo, yend = yhi, colour = group), linewidth = 0.7) +
    geom_point(data = marks, inherit.aes = FALSE, aes(x = x, y = y, colour = group),
               shape = 21, fill = "white", size = 2.8, stroke = 1.1) +
    geom_vline(data = div, aes(xintercept = x), linetype = "dotted", colour = "grey30") +
    geom_text(data = div, inherit.aes = FALSE, aes(x = x, y = Inf, label = paste0("differ from ", round(x, 1), " h")),
              hjust = -0.08, vjust = 1.5, size = 3.3) +
    facet_wrap(~panel, ncol = 1, scales = "free_y") +
    scale_colour_manual(values = allcols, breaks = names(allcols), limits = names(allcols), name = NULL) +
    scale_fill_manual(values = allcols, breaks = names(allcols), limits = names(allcols), name = NULL) +
    scale_x_continuous(breaks = sort(unique(res$data$time))) +
    labs(x = xlab, y = NULL, subtitle = subtitle,
         caption = paste0(
           "GAM fits with 95% bands (simultaneous for the difference), including smoothing uncertainty; both groups start\n",
           "at the dose measured at time 0. Top: end of the lag (dotted tangent at the steepest point) and maximum load.\n",
           "Middle: fastest growth = onset of control, and time to control (growth below ", 100 * res$control_frac,
           "% of its maximum). Bottom: difference;\n",
           "dotted line: from when the simultaneous band excludes 0. Bars: 95% intervals from ", format(res$n_sim, big.mark = ","),
           " draws. n.r.: not reached by ", max(res$tgrid), " h.\n",
           "Above: ", lv[2], " - ", lv[1], " [95% interval] for the five parameters examined; differ = 95% interval excludes 0;\n",
           "equivalent = 90% interval within the margin set beforehand; otherwise inconclusive.")) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top", strip.background = element_blank(),
          strip.text = element_text(hjust = 0, face = "bold"),
          plot.subtitle = element_text(family = "mono", size = 7.5),
          plot.caption = element_text(hjust = 0, size = 8))
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

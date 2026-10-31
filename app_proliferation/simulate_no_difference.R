# Checks how often the Fly Proliferation analysis finds a difference between two groups when
# there is none. Both groups are drawn from the same curve: the curve fitted to the females of the
# example, with the within-group SD observed at each time, on the protocol's design (8 time zero
# tubes, 8 flies per group per hour from 1 to 12 h). Each test should call about 5% of these
# experiments different. Run from the repository folder: Rscript tools/simulate_no_difference.R
source("app_proliferation/prolif_core.R")
a <- read.table("app_proliferation/example_proliferation_CS_sexdim.csv", header = TRUE, sep = ";")
z <- read.table("app_proliferation/example_time0_CS_sexdim.csv", header = TRUE, sep = ";")
a <- a[a$Status == "Detected" & a$Sex == "Female", ]
d <- rbind(data.frame(time = 0, log2_load = z$log2_CFU_per_fly), data.frame(time = a$Time_post_injection, log2_load = a$log2_CFU_per_fly))
truth <- gam(log2_load ~ s(time, k = 6), data = d, method = "REML")
sd_t <- tapply(residuals(truth), d$time, sd)
mu <- function(t) predict(truth, data.frame(time = t))

R <- 200                      # experiments simulated
set.seed(2026); out <- NULL
for (r in seq_len(R)) {
  t <- rep(1:12, each = 8)
  x <- rbind(data.frame(time = 0, group = "A", log2_load = rnorm(8, mu(0), sd_t[["0"]])),
             data.frame(time = t, group = "A", log2_load = rnorm(length(t), mu(t), sd_t[as.character(t)])),
             data.frame(time = t, group = "B", log2_load = rnorm(length(t), mu(t), sd_t[as.character(t)])))
  x$group <- factor(x$group, levels = c("A", "B"))
  res <- prolif_fit(x, k = 6, n_sim = 2000, seed = r)
  p <- prolif_summary(res)$params
  out <- rbind(out, data.frame(difference_smooth = res$smooth_test$p[2] < 0.05, band = nrow(res$band_excl) > 0,
                               t(setNames(p$diff_low > 0 | p$diff_high < 0, p$key))))
}
cat("Share of", R, "experiments with no difference that were called different (%):\n")
print(round(100 * colMeans(out, na.rm = TRUE), 1))

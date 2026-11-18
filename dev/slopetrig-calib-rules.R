# Lane slopetrig: candidate rules for when the default engages, scored
# on the base arm of both calibration grids. Each rule is "spread < 0.05
# OR <condition>"; the ratio is over every row (ratio_all, recomputed
# from the design: equal to the nonzero-row ratio except on the dose,
# two-value and ulp designs, whose zero rows bring it to about 1.4).
# Prints per rule: engaged, short under autoscale = FALSE left
# unengaged (and the largest), engaged where autoscale = FALSE was not
# short, and short under autoscale = TRUE where engaged.
# Usage: Rscript slopetrig-calib-rules.R <tsv> [<tsv> ...]
fs <- commandArgs(TRUE)
d <- do.call(rbind, lapply(fs, utils::read.delim))
best <- pmax(d$ll_ref, d$ll_F, d$ll_T, d$ll_D, na.rm = TRUE)
d$sh_F <- best - d$ll_F
d$sh_T <- best - d$ll_T
zero_rows <- d$design %in% c("dose", "twonz", "ulpind")
# x is 0 on half the rows there: rms_all^2 = rms_nz^2 / 2, and
# var_all = rms_nz^2 / 2 - (mean_nz / 2)^2 with mean_nz^2 =
# rms_nz^2 - spread^2 (n large; the sample sd's n - 1 is ignored)
m_nz <- sqrt(pmax(d$rms^2 - d$spread^2, 0))
d$ratio_all <- ifelse(zero_rows,
                      sqrt(d$rms^2 / 2) / sqrt(d$rms^2 / 2 - (m_nz / 2)^2),
                      d$ratio)
sp <- d$spread < 0.05
rules <- list(
  "none (spread rule only)" = rep(FALSE, nrow(d)),
  "ratio_all > 2" = d$ratio_all > 2,
  "ratio_all > 2 & rms > 2" = d$ratio_all > 2 & d$rms > 2,
  "ratio_all > 3" = d$ratio_all > 3,
  "ratio_all > 5" = d$ratio_all > 5,
  "ratio_all > 10" = d$ratio_all > 10,
  "ratio_all > 20" = d$ratio_all > 20,
  "ratio_nz > 2" = d$ratio > 2)
k <- 1e-3
cat("cases", nrow(d), "; spread rule alone engages", sum(sp), "\n")
for (nm in names(rules)) {
  eng <- sp | rules[[nm]]
  miss <- d$sh_F > k & !eng
  cat(sprintf(paste0("%-26s engaged %4d (by rule only %4d) ; short F left",
                     " %3d (max %6.3g, silent %3d) ; engaged not needed",
                     " %4d ; short T where engaged %d\n"),
              nm, sum(eng), sum(eng & !sp), sum(miss),
              if (any(miss)) max(d$sh_F[miss]) else 0,
              sum(miss & d$w_F == 0), sum(eng & d$sh_F <= k),
              sum(eng & d$sh_T > k)))
}
cat("dose-type designs (dose, twonz, ulpind): ratio_all range",
    signif(range(d$ratio_all[zero_rows]), 3), "\n")

# Reviewer of lane surface, claim 4: frm_multiple() pooling against
# mice::pool(lm) (estimates, SE, Barnard-Rubin df, CI) and against
# brm_multiple()'s combined posterior (dev/surface-rev-brms.R, same
# imputations: nhanes, m = 5, mice seed 1).
#   Rscript dev/surface-rev-multiple.R > dev/surface-rev-out/multiple.txt
source("dev/surface-rev-env.R")
rev_env("lane")
suppressPackageStartupMessages(library(frmtmb))
data("nhanes", package = "mice")
imp <- mice::mice(nhanes, m = 5, print = FALSE, seed = 1)
fm <- frm_multiple(bmi ~ age * chl, data = imp)
fx <- fixef(fm)
print(fx)
mp <- summary(mice::pool(with(imp, lm(bmi ~ age * chl))), conf.int = TRUE)
print(mp)
pl <- fm$pooled
cat("\nfm$pooled columns:", names(pl), "\n")
print(pl)
cat("\nmax |Estimate diff| vs mice:",
    format(max(abs(fx[, "Estimate"] - mp$estimate))), "\n")
cat("SE ratio frmtmb / mice:", format(fx[, "Est.Error"] / mp$std.error,
                                      digits = 6), "\n")
if ("df" %in% names(pl)) {
  cat("df frmtmb:", format(pl$df[1:4], digits = 6), "\n")
  cat("df mice  :", format(mp$df, digits = 6), "\n")
}
cat("CI lower frmtmb:", format(fx[, 3], digits = 6), "\n")
cat("CI lower mice  :", format(mp$`2.5 %`, digits = 6), "\n")
cat("df.residual(fit 1):", df.residual(fm$fits[[1]]),
    "| lm:", df.residual(lm(bmi ~ age * chl, mice::complete(imp, 1))), "\n")
s <- summary(fm)
print(s)

br <- "dev/surface-rev-out/brms-multiple.rds"
if (file.exists(br)) {
  b <- readRDS(br)
  cat("\nbrm_multiple combined posterior (10 chains, 20000 draws):\n")
  print(b$fixef)
  cat("Estimate: frmtmb - brms posterior mean, in brms posterior SDs:",
      format((fx[, 1] - b$fixef[, 1]) / b$fixef[, 2], digits = 3), "\n")
  cat("SE ratio frmtmb pooled / brms posterior SD:",
      format(fx[, 2] / b$fixef[, 2], digits = 4), "\n")
  # brms's combined draws are the mixture of the m posteriors, whose
  # variance is the mean within + (m - 1) / m * between, not Rubin's
  # (1 + 1 / m) * between
  dr <- b$draws
  cat("draws columns:", head(names(dr), 6), "\n")
  ce_f <- conditional_effects(fm, "chl")[[1]]
  ce_b <- b$ce
  cat("CE grid: frmtmb", nrow(ce_f), "rows, chl range",
      format(range(ce_f$chl)), "| brms", nrow(ce_b), "rows, chl range",
      format(range(ce_b$chl)), "\n")
  cat("CE conditions: age frmtmb", unique(ce_f$age), "| brms",
      unique(ce_b$age), "\n")
  i <- match(signif(ce_b$chl, 8), signif(ce_f$chl, 8))
  if (!anyNA(i)) {
    cat("CE estimate: max |frmtmb - brms| / brms se:",
        format(max(abs(ce_f$estimate__[i] - ce_b$estimate__) / ce_b$se__),
               digits = 3), "\n")
    cat("CE band width ratio frmtmb / brms: median",
        format(median((ce_f$upper__[i] - ce_f$lower__[i]) /
                        (ce_b$upper__ - ce_b$lower__)), digits = 4), "\n")
  } else cat("CE grids do not match point for point\n")
}
# imputation ranges of chl, the grid question
cat("\nchl range per imputation:\n")
for (k in 1:5) cat(" ", k, format(range(mice::complete(imp, k)$chl)), "\n")

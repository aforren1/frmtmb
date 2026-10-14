# Reviewer, lane dblbar: fixed effects of (1 + f || g), frmtmb against
# lme4 dummies and glmmTMB diag(). Seed 20261013 (dev/dblbar-data.R).
.libPaths(c("C:/Users/adf44/source/r/wt-dblbar-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages({
  library(frmtmb); library(lme4); library(glmmTMB)
})
dd <- dblbar_data()
dd$fb <- as.numeric(dd$f == "b")
dd$fc <- as.numeric(dd$f == "c")
ff <- frm(bf(y ~ f + (1 + f || g)), family = gaussian(), data = dd)
print(fixef(ff))
lm4 <- lmer(y ~ f + (1 | g) + (0 + fb | g) + (0 + fc | g), data = dd,
            REML = FALSE)
gt <- glmmTMB(y ~ f + diag(1 + f | g), data = dd, REML = FALSE)
b <- fixef(ff)[, 1]
cat(sprintf("%.8f ", b), "\n")
cat(sprintf("%.8f ", fixef(lm4)), "\n")
cat(sprintf("%.8f ", fixef(gt)$cond), "\n")
cat("se frm", sprintf("%.6f ", fixef(ff)[, 2]), "\n")
cat("se lme4", sprintf("%.6f ", sqrt(diag(as.matrix(vcov(lm4))))), "\n")
cat("se glmmTMB", sprintf("%.6f ", sqrt(diag(vcov(gt)$cond))), "\n")

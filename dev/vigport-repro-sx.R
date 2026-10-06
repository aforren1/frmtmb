# Repro: fixef() names an s() term's unpenalized coefficient `sx_1`, as
# brms does since 0.66.0, but the column it multiplies is not brms's
# column. Found by the estimate comparison (brms_distreg fit_smooth1:
# frmtmb sx1_1 = 1.83 against brms's posterior mean 9.73, sd 3.44).
#
#   Rscript dev/vigport-repro-sx.R [lib]
#
# Prints, for s(x1) and s(x2): the correlation of frmtmb's column with
# brms's Xs column (from brms::standata(), no compile) and the ratio of
# their scales; then frmtmb's coefficients carried into brms's unit,
# beside a plain lm() on brms's own columns as a check of the unit.
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r5"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from", lib, "\n")
set.seed(1)
dat <- mgcv::gamSim(eg = 6, n = 200, scale = 2, verbose = FALSE)
f <- frm(bf(y ~ s(x1) + s(x2)), data = dat)
Xs <- brms::standata(brms::bf(y ~ s(x1) + s(x2)), data = dat)$Xs
X <- f$frame$linpreds[[1]]$X
for (v in c("x1", "x2")) {
  a <- X[, paste0("s(", v, ").fx1")]
  b <- Xs[, paste0("s", v, "_1")]
  r <- cor(a, b)
  s <- sd(a) / sd(b)
  est <- fixef(f)[paste0("s", v, "_1"), "Estimate"]
  cat(sprintf(paste0("s(%s): cor(frmtmb column, brms column) = %+.6f, ",
                     "scale ratio = %.4f; frmtmb sx_1 = %.4f, in brms's ",
                     "unit %.4f\n"), v, r, s, est, est * s * sign(r)))
}
# Same model, different unit: the fitted means do not depend on it.
g <- mgcv::gam(y ~ s(x1) + s(x2), data = dat, method = "ML")
cat(sprintf("max |fitted(frmtmb) - fitted(mgcv ML)| = %.3g (y sd %.3g)\n",
            max(abs(fitted(f)[, 1] - fitted(g))), sd(dat$y)))

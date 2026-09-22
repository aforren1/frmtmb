# An OUTSIDE judge for the ordinal threshold rows MAJOR 3 added:
# MASS::polr fits the same model with its own machinery and reports
# its own standard errors, so nothing on its path comes through
# frmtmb. polr parameterizes logit P(Y <= k) = zeta_k - x'b, brms and
# frmtmb write the same model as eta_k = tau_k + x'b with the
# cumulative link on P(Y > k), so the POINT estimates differ by a sign
# convention and the STANDARD ERRORS must agree.
#
#   Rscript dev/shapes-p1-ordjudge.R

.libPaths(c("C:/Users/adf44/source/r/shapes-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))

set.seed(20260917)
n <- 120
dd <- data.frame(x = rnorm(n), z = rnorm(n))
dd$ord <- factor(cut(1 + 0.8 * dd$x + rnorm(n), c(-Inf, -0.3, 0.8, Inf),
                     labels = 1:3), ordered = TRUE)

fit <- frm(bf(ord ~ x) + cumulative(), data = dd)
fe <- fixef(fit)
pm <- MASS::polr(ord ~ x, data = dd, Hess = TRUE, method = "logistic")
sp <- summary(pm)$coefficients

cat("frmtmb fixef()\n")
print(round(fe, 6))
cat("\nMASS::polr\n")
print(round(sp[, 1:2], 6))

# fixef() row order is brms's: the thresholds first, then the slope
pr <- c("1|2", "2|3", "x")
cmp <- data.frame(
  row = rownames(fe),
  polr_row = pr,
  frm_est = unname(fe[, "Estimate"]),
  polr_est = unname(sp[pr, 1]),
  frm_se = unname(fe[, "Est.Error"]),
  polr_se = unname(sp[pr, 2])
)
cmp$est_diff <- cmp$frm_est - cmp$polr_est
cmp$se_ratio <- cmp$frm_se / cmp$polr_se
cat("\nfrmtmb against polr\n")
print(cmp, row.names = FALSE, digits = 7)

ok <- max(abs(cmp$se_ratio - 1)) < 1e-3 && max(abs(cmp$est_diff)) < 1e-3
cat("\nmax |se ratio - 1| =", signif(max(abs(cmp$se_ratio - 1)), 3),
    "  max |est diff| =", signif(max(abs(cmp$est_diff)), 3), "\n")
cat(if (ok) "AGREES with polr\n" else "DISAGREES with polr\n")
if (!ok) quit(status = 1L)

# Lane catmulti, fix round 2: what the reference twin's frame costs in a
# refcat = NA fit on 1e5 rows. Times (elapsed, minimum of 3 rounds,
# interleaved) the whole NA fit, the twin's count alone
# (refcat_first_outer() with the fit's own arguments), and a CONTROL:
# the same NA fit's frame by frm(dry_run = "frame"), which the twin's
# count should cost about the same as. seed 101.
# Output: dev/catmulti-out/twin-time.txt
.libPaths(c("C:/Users/adf44/source/r/wt-catmulti-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(101)
n <- 1e5
d <- data.frame(x = rnorm(n), x2 = exp(rnorm(n)))
E <- cbind(0, 0.3 + 0.8 * d$x, -0.4 - 0.5 * log(d$x2))
P <- exp(E) / rowSums(exp(E))
u <- runif(n)
d$y <- factor(c("a", "b", "c")[1L + (u > P[, 1]) + (u > P[, 1] + P[, 2])])
form <- bf(y ~ poly(x, 2) + log(x2))
fam <- categorical(refcat = NA)
el <- function(expr) system.time(expr)[["elapsed"]]
t_fit <- t_twin <- t_ctrl <- numeric(3)
for (r in 1:3) {
  t_fit[r] <- el(fit <- suppressWarnings(frm(form, family = fam, data = d)))
  bform <- ns$resolve_deferred_families(ns$as_bform(form, fam), d)
  spec <- ns$parse_spec(bform)
  t_twin[r] <- el(rec <- ns$refcat_first_outer(spec, bform, d,
                                                stats::na.omit, list(),
                                                TRUE))
  t_ctrl[r] <- el(frm(form, family = fam, data = d, dry_run = "frame"))
}
cat(sprintf("n = %d rows; elapsed seconds, minimum of 3 rounds\n", n))
cat(sprintf("whole refcat = NA fit: %.2f (rounds %s)\n", min(t_fit),
            paste(sprintf("%.2f", t_fit), collapse = ", ")))
cat(sprintf("twin count alone:      %.2f (rounds %s)\n", min(t_twin),
            paste(sprintf("%.2f", t_twin), collapse = ", ")))
cat(sprintf("control, the NA frame: %.2f (rounds %s)\n", min(t_ctrl),
            paste(sprintf("%.2f", t_ctrl), collapse = ", ")))
cat(sprintf("twin share of the fit: %.1f%%\n",
            100 * min(t_twin) / min(t_fit)))
cat("df", attr(logLik(fit), "df"), " AIC", AIC(fit), "\n")

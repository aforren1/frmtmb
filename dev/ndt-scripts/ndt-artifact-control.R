# ndt lane: design B with BOTH parameterizations.
#
# Design B is the falsification design for the per-group bound: `ndt`
# is CONSTANT across subjects and the boundary separation varies, so
# the per-subject floors vary and the truth does not. The per-group fit
# reports a spread that is not there and a population `ndt` biased low.
# The question this answers is whether the GLOBAL bound does better on
# the same data, because a limitation only counts against the change if
# the thing it replaced did not have it.
#
# 15 subjects x 250 trials, seeds 771 + 0:7. `ndt` 0.25 for every
# subject, log boundary spread 0.25.
# Run: Rscript --vanilla dev/ndt-scripts/ndt-artifact-control.R
.libPaths(c("C:/Users/adf44/source/r/ndt-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
})

ns <- 15L
nt <- 250L

one <- function(seed) {
  set.seed(seed)
  u_nd <- rnorm(ns, 0, 0)
  u_bs <- rnorm(ns, 0, 0.25)
  s <- rep(seq_len(ns), each = nt)
  t0 <- 0.25 * exp(u_nd)
  d <- ddm_simulate(ns * nt, mu = 1.1, bs = 1.4 * exp(u_bs),
                    ndt = t0[s], bias = 0.5)
  d$s <- factor(s)
  key <- d[match(levels(d$s), as.character(d$s)), , drop = FALSE]
  fl <- as.numeric(tapply(d$rt, d$s, min))
  out <- numeric(0)
  for (grp in c(TRUE, FALSE)) {
    form <- if (grp) {
      bf(rt | dec(upper) + ndt_group(s) ~ 1, bs ~ 1 + (1 | s),
         ndt ~ 1 + (1 | s), bias = 0.5)
    } else {
      bf(rt | dec(upper) ~ 1, bs ~ 1 + (1 | s), ndt ~ 1 + (1 | s),
         bias = 0.5)
    }
    fit <- frm(form, family = wiener(), data = d)
    hat <- as.numeric(ndt_time(fit, newdata = key))
    tag <- if (grp) "pg" else "gl"
    out <- c(out, stats::setNames(
      c(fit$opt$convergence, 1000 * mean(hat - t0), 1000 * sd(hat),
        1000 * sqrt(mean((hat - t0)^2)), as.numeric(logLik(fit))),
      paste0(tag, c("_conv", "_bias_ms", "_sd_ms", "_rmse_ms", "_ll"))))
  }
  c(out, floor_sd_ms = 1000 * sd(fl), floor_min_ms = 1000 * min(fl))
}

m <- t(vapply(0:7, function(k) one(771L + k), numeric(12)))
cat("== design B, ndt CONSTANT at 0.25, log bs spread 0.25 ==\n")
cat("   ", ns, "subjects x", nt, "trials, 8 seeds\n")
for (nm in colnames(m)) {
  cat(sprintf("  %-14s mean %10.4f  sd %9.4f  range %10.4f %10.4f\n",
              nm, mean(m[, nm]), stats::sd(m[, nm]), min(m[, nm]),
              max(m[, nm])))
}
cat("\nper-seed log-likelihood difference, per-group less global:\n")
print(round(m[, "pg_ll"] - m[, "gl_ll"], 3))

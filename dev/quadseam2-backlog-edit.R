# Lane quadseam2's edit of dev/test-backlog.md: CLOSED lines under the
# items it closes, and its filings in "Filed in the round of 2026-11-13".
# Run once from the worktree root:  Rscript dev/quadseam2-backlog-edit.R
f <- "dev/test-backlog.md"
x <- readLines(f)
pre <- "  CLOSED in the round of 2026-11-13 (lane quadseam2): "
item_end <- function(start) {
  i <- start + 1L
  while (i <= length(x) && !grepl("^(- |#)", x[i]) && nzchar(x[i])) {
    i <- i + 1L
  }
  i - 1L
}
close <- list(
  "561" = c(
    "not reproduced. 0 breakdowns in 240 one-block fits with a group sd",
    "  of 0 or 0.05 over six families, 145 of them ending below an sd of",
    "  1e-3 (dev/quadseam2-break.R). The rule's anchor-dependent error is",
    "  4846, set by TMBad's integrator tolerances (upstream TMB-7)."),
  "4836" = c(
    "misdiagnosed: the stall stop had the group sd at exp(-20.4), and the",
    "  rule was right there to 4e-9 (dev/quadseam2-snref.R). The sd escape of",
    "  4861 restarts it; it ends 3.8e-7 below the default start's fit."),
  "4861" = c(
    "`quad_sd_up_escape()` probes rules taped at each larger sd and",
    "  restarts on a rule taped at the start, settled by `quad_settle()`."),
  "4867" = c(
    "the limit rule runs on the one-block route through `se_probe_fn()`,",
    "  its curvature by differences of a rule taped at the moved point: SE",
    "  NaN and the Laplace warning (test-quadrature-defects.R)."),
  "5222" = c(
    "`escape_stationary()` asks again from the kept optimum while it",
    "  gains, at most once per predictor on its point; hmm and mixture",
    "  seeds 1 to 3 reach the default start's fit (dev/quadseam2-hmm.R)."),
  "5408" = c(
    "the rule's escape runs only when the Laplace-side one did not",
    "  (dev/quadseam2-escdup.R); `stationary_escape_keep()` is removed."),
  "5620" = c(
    "`latres_quadratic()` drops the warning where every row's log density",
    "  is quadratic in eta: lognormal 0 of 10 fits (1 on 0.82.0), bernoulli",
    "  3 of 3 kept (dev/quadseam2-latres.R).")
)
starts <- as.integer(names(close))
stopifnot(all(grepl("^- ", x[starts])))
for (s in sort(starts, decreasing = TRUE)) {
  e <- item_end(s)
  txt <- close[[as.character(s)]]
  txt[1] <- paste0(pre, txt[1])
  x <- append(x, txt, after = e)
}
new <- c(
  "## Filed in the round of 2026-11-13",
  "",
  "### From lane quadseam2",
  "",
  "- **MEDIUM. The one-block rule's anchor error is the integrator's",
  "  tolerance, not its calibration** (item 4846). A tape whose integrand",
  "  has its exact Laplace mode and scale, so that TMB's `dx = 1`",
  "  differences and `ytol` cannot matter, has a worst log-likelihood error",
  "  of 1.19e-3 on the 120 fits of 4846, against 1.06e-3 for TMB's own",
  "  calibration; better by more than 1e-6 on 29, worse on 21",
  "  (dev/quadseam2-white.R). TMBad's `integrate()` stops at `reltol =",
  "  abstol = 1e-4` per integral, and neither those nor `dx` and `ytol`",
  "  reach it from an RTMB spec (upstream TMB-7), so route (b) is closed.",
  "  What remains is route (a): a Gauss-Hermite rule on the tape with",
  "  enough nodes, as the nested route writes it. Size L.",
  ""
)
at <- grep("^## Filed at the 0.82.0 release", x)
stopifnot(length(at) == 1L)
x <- append(x, new, after = at - 1L)
writeLines(x, f)
cat("edited", f, "\n")

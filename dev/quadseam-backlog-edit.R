# One-off: the lane's CLOSED lines and filed items in dev/test-backlog.md
# (lane-env r082, "Backlog"). Checks each anchor line before inserting,
# from the bottom up so earlier line numbers hold.
f <- "C:/Users/adf44/source/r/frmtmb-wt-quadseam/dev/test-backlog.md"
x <- readLines(f)
ins <- function(x, after, expect, lines) {
  stopifnot(identical(x[after], expect))
  append(x, lines, after = after)
}
filed <- c(
  "## Filed in the round of 2026-11-11",
  "",
  "### From lane quadseam",
  "",
  "Lane record: `dev/quadseam-findings.md`.",
  "",
  "- **A skew normal under `quadrature = TRUE` started at the stall can",
  "  end 5.5e-2 short.** wt-skewinit's acceptance seed 6 with `(1 | g)`",
  "  and the stall start: logLik -244.155229735 at alpha 5.008, against",
  "  -244.100465993 at 5.267 from the default start; the Laplace fit",
  "  reaches its maximum from both starts. The same on 0.81.0 and the",
  "  lane build; no escape fires, since the warm start does not stop on",
  "  alpha = 0. The one-block route has no Newton polish",
  "  (`newton_polish()` runs only without `integrate`), and the rule taped",
  "  at the stop reads a gradient below `grad_tol` there",
  "  (`dev/quadseam-escorder.R`, `dev/quadseam-seed6.R`).",
  "- **The Gauss-Kronrod rule's error depends on its anchor, up to",
  "  9e-7.** On 10 Bernoulli fits with groups of 4 and sd 1.2 the fit's",
  "  log-likelihood against `stats::integrate()` per group at the same",
  "  estimate: 1.7e-9 to 9.0e-7 when the rule is taped at the estimate,",
  "  7.5e-9 to 3.9e-7 when taped at the Laplace optimum",
  "  (`dev/quadseam-bernref.R`). TMB's calibration takes its curvature",
  "  by differences with a step of 1 and stops its mode search at a gain",
  "  of 1e-2 (`ytol`), so the rule is not equally good at every anchor.",
  "- **The nested quadrature route still gives no boundary verdict.** An",
  "  sd at zero under `(1 | ga/gb)` keeps the \"flat\" warning: the upward",
  "  check needs the rule taped at each probed sd, and the nested rule",
  "  needs the Laplace objective for its anchor, which the fit does not",
  "  keep (`se_edge_probed()`).",
  "- **`sd_up_escape()` does not run under quadrature.** A quadrature fit",
  "  left with an sd near zero short of an interior maximum (the start at",
  "  exp(-20) of test-se-tier.R) keeps the \"flat\" warning or the",
  "  convergence warning where the Laplace fit restarts and reaches the",
  "  maximum. The upward probe now works on the one-block route",
  "  (`se_probe_fn()`); the restart would then need `quad_settle()`.",
  "")
i <- grep("^## Filed at the 0.81.0 release \\(2026-11-10\\)$", x)
stopifnot(length(i) == 1L)
x <- append(x, filed, after = i - 1L)
cl <- function(s) paste0("  CLOSED in the round of 2026-11-11 (lane quadseam): ", s)
x <- ins(x, 3618, "  unbounded quadrature fits are not detected.", c(
  cl("the one-block rule is taped"),
  "  again at the estimate until it is stationary (`quad_settle()`); a",
  "  runaway warns `frmtmb_quad_unsettled`."))
x <- ins(x, 3481, "  boundary rule do not run on a quadrature objective.", c(
  cl("on the one-block route the"),
  "  boundary rule probes rules taped at each sd (`se_probe_fn()`): the",
  "  same message and NaN SE as under Laplace. The nested route is",
  "  filed under \"Filed in the round of 2026-11-11\"."))
x <- ins(x, 3473, "  spurious maximum.", c(
  cl("`quad_settle()` re-tapes the"),
  "  rule at the estimate until it is stationary; the two-response",
  "  runaway that stopped at log sigma -4.0 now warns",
  "  `frmtmb_quad_unsettled` (`dev/quadseam-runaway.R`)."))
x <- ins(x, 836, "  are changes to `quad_fit()` rather than to the escape.", c(
  cl("`quad_fit()` escapes on the"),
  "  Laplace objective before it tapes the rule, and re-tapes after any",
  "  later move; the case lands 1.0e-8 from the unstalled fit",
  "  (`dev/quadseam-punch1c.R`)."))
writeLines(x, f)
cat("done\n")

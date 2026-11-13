# One-off for the fix round: reword the 3466 and 3614 CLOSED lines, give
# the anchor-error item the measured range, and file the review's MINOR
# 7 (dev/test-backlog.md). Each anchor is checked before it is replaced.
f <- "C:/Users/adf44/source/r/frmtmb-wt-quadseam/dev/test-backlog.md"
x <- readLines(f)
swap <- function(x, old, new) {
  i <- which(x == old[1])
  i <- i[vapply(i, function(j) identical(x[j:(j + length(old) - 1L)], old),
                NA)]
  stopifnot(length(i) == 1L)
  c(x[seq_len(i - 1L)], new, x[-seq_len(i + length(old) - 1L)])
}
x <- swap(x, c(
  "  CLOSED in the round of 2026-11-11 (lane quadseam): `quad_settle()`",
  "  re-tapes the rule at the estimate until it is stationary; the two-response",
  "  runaway that stopped at log sigma -4.0 now warns",
  "  `frmtmb_quad_unsettled` (`dev/quadseam-runaway.R`)."), c(
  "  CLOSED in the round of 2026-11-11 (lane quadseam): `quad_settle()`",
  "  re-tapes the rule at the estimate until it is stationary, and warns",
  "  `frmtmb_quad_unsettled` on a runaway the weights check does not",
  "  reach: a gaussian with `mi()` stopped at log sigma -4.06 with code 0",
  "  (`dev/quadseam-runaway3.R`). The weighted `(1 | obs)` cases of this",
  "  item are refused by the weights check since lane wcheck (3447)."))
x <- swap(x, c(
  "  CLOSED in the round of 2026-11-11 (lane quadseam): the one-block rule is taped",
  "  again at the estimate until it is stationary (`quad_settle()`); a",
  "  runaway warns `frmtmb_quad_unsettled`."), c(
  "  CLOSED in the round of 2026-11-11 (lane quadseam): the one-block rule",
  "  is taped again at the estimate until it is stationary",
  "  (`quad_settle()`), and an unbounded fit the weights check does not",
  "  reach (`mi()`, `se()`, a family outside its table) warns",
  "  `frmtmb_quad_unsettled`; the gaussian `(1 | obs)` case is refused by",
  "  the weights check."))
x <- swap(x, c(
  "- **The Gauss-Kronrod rule's error depends on its anchor, up to",
  "  9e-7.** On 10 Bernoulli fits with groups of 4 and sd 1.2 the fit's",
  "  log-likelihood against `stats::integrate()` per group at the same",
  "  estimate: 1.7e-9 to 9.0e-7 when the rule is taped at the estimate,",
  "  7.5e-9 to 3.9e-7 when taped at the Laplace optimum",
  "  (`dev/quadseam-bernref.R`). TMB's calibration takes its curvature",
  "  by differences with a step of 1 and stops its mode search at a gain",
  "  of 1e-2 (`ytol`), so the rule is not equally good at every anchor."), c(
  "- **The Gauss-Kronrod rule's error depends on its anchor, up to",
  "  1e-3.** On the review's 120 fits at the default `grad_tol` (Bernoulli",
  "  groups of 4 and of 2, poisson sd 1 and sd 2;",
  "  `dev/reviews/quadseam-review-acc.R`) the log-likelihood against",
  "  `stats::integrate()` per group at the fit's own estimate has a worst",
  "  error of 1.06e-3 with the rule taped at the estimate (4.48e-2 on",
  "  0.81.0); against 0.81.0 it is worse by more than 1e-6 on 19 fits and",
  "  better on 17. TMB's calibration takes its curvature by differences",
  "  with a step of 1 and stops its mode search at a gain of 1e-2",
  "  (`ytol`), so the rule is not equally good at every anchor."))
add <- c(
  "- **Cross-lane (boundfam): a coefficient that runs off under one-block",
  "  quadrature gets a finite SE and no warning.** Poisson with a factor",
  "  level whose counts are all zero, plus `(1 | g)`: under Laplace the",
  "  coefficient (-22.53) has SE NaN and a warning (\"toward -Inf\" on",
  "  boundfam's build, \"flat\" on this lane's); under `quadrature = TRUE`",
  "  it has SE 7471 and no word, on both builds",
  "  (`dev/reviews/quadseam-review-limit.R`). boundfam's",
  "  `se_edge_limit()` keeps its own `!isTRUE(fit$quadrature)` gate. On",
  "  the one-block route its direction test can read `se_probe_fn()`;",
  "  its curvature test (`se_limit_curv()`) takes `obj$he()` at a moved",
  "  point, which on the frozen rule is a read away from the anchor and",
  "  would need the Hessian of a rule taped at that point. For the",
  "  consolidation to decide; not done in either lane.")
i <- which(x == "  (`se_probe_fn()`); the restart would then need `quad_settle()`.")
stopifnot(length(i) == 1L)
x <- append(x, add, after = i)
writeLines(x, f)
cat("done\n")

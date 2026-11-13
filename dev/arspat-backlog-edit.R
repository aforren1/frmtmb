# Lane arspat: add the CLOSED lines under the items this lane closed and
# the lane's filed items, to dev/test-backlog.md as of e8abed0a. Run once;
# it refuses to run on a file it has already changed.
f <- "dev/test-backlog.md"
x <- readLines(f, warn = FALSE)
if (any(grepl("(lane arspat)", x, fixed = TRUE))) stop("already edited")
tag <- "CLOSED in the round of 2026-11-11 (lane arspat): "
# item first line -> its last line, at e8abed0a, and the one line
closed <- list(
  list(first = 3897, last = 3900, txt = paste0(
    "`sar()` takes a W with complex eigenvalues; each pair adds ",
    "log|1 - rho lambda|^2, rho's interval from the real parts ",
    "(dev/arspat-findings.md).")),
  list(first = 4104, last = 4107, txt = paste0(
    "a message of class `frmtmb_autocor_gap`; the levels reading is kept, ",
    "a question for the user in dev/arspat-findings.md.")),
  list(first = 4108, last = 4109, txt = paste0(
    "`se(sigma = TRUE)` and `mi()` beside the covariance form fit as in ",
    "brms; L1 to 1.1e-13 and 1.7e-13 (dev/arspat-findings.md).")),
  list(first = 4110, last = 4114, txt = paste0(
    "a slot giving O K is accepted by its Gram matrix and read through ",
    "K'K e (`jacobian_rotated`).")),
  list(first = 4115, last = 4119, txt = paste0(
    "`options(frmtmb.autocor_dense_max = )`; frmtmb.spatial's test sets it ",
    "to 60.")),
  list(first = 4120, last = 4121, txt = paste0(
    "fcor() reads diag(M^-1) from the block: peak 34.5 Mb to 0.2 Mb and ",
    "4.0x faster per draw at n = 1500 (dev/arspat-fcor-cost.R).")),
  list(first = 4172, last = 4185, txt = paste0(
    "fits, `D R D + diag(se^2)` per group; objective, draws and pointwise ",
    "log-density (dev/arspat-findings.md).")),
  list(first = 4299, last = 4306, txt = paste0(
    "`vcov_cluster()` and `vcov = ~ g` work when every series lies in one ",
    "cluster; a split series is refused.")),
  list(first = 4320, last = 4324, txt = paste0(
    "`get_predict()` predicts marginaleffects' padding rows apart under ",
    "cov = FALSE; the HAC refusal names the within-group order and `~ g`."))
)
stopifnot(startsWith(x[3897], "- `sar()` refuses a matrix with complex"),
          startsWith(x[4172], "- **PARITY GAP: `se()` beside"),
          startsWith(x[4299], "- **PARITY GAP: a cluster sandwich"),
          startsWith(x[4320], "- **`predictions(newdata =)` on an"))
filed <- c(
  "## Filed in the round of 2026-11-11",
  "",
  "### From lane arspat",
  "",
  "- `sar()` beside `ar(cov = FALSE)` is still refused, where brms 2.23.0",
  "  fits it (`mu[n] += Err[n, 1:Kar] * ar`, then the SAR density of `Y`;",
  "  `dev/arspat-out/brms-probe.txt`). A response holds one residual term",
  "  (`frame$autocor[[resp]]`), and 33 readers in 11 files branch on",
  "  `autocor_is_cond()` of that one slot, so carrying a `cov = FALSE`",
  "  term beside a registered one is a change of the frame's shape, not",
  "  of `sar()`. The `mi()` half of backlog 3903 is closed; `newdata2`",
  "  stays open too.",
  "- `confint_varcorr()` still gives `sar()`'s rho a raw-scale",
  "  delta-method interval (backlog 3901). An interval on the logit chart",
  "  of the term's bounds needs `lo`/`hi` carried through the five",
  "  places that rbind or back-transform these rows (R/confint.R,",
  "  R/methods-fit.R twice, R/multiple.R, R/autocor-latent.R).",
  "- On gaussian data a `student()` fit runs nu to about 1e9, where brms's",
  "  program and frmtmb's density part by 1e-5 nats (`se_ar_student` of",
  "  `dev/arspat-out/l1-all.txt`); at df 4 they agree to 1.1e-13. Likely",
  "  Stan's `multi_student_t` head at a large nu (frmtmb uses",
  "  `lgamma_shift_diff()`); not checked against an exact reference.",
  "")
anchor <- which(x == "## Filed at the 0.81.0 release (2026-11-10)")
stopifnot(length(anchor) == 1L)
# from the bottom up, so the line numbers above stay valid
ord <- order(vapply(closed, `[[`, 0, "last"), decreasing = TRUE)
for (i in ord) {
  c1 <- closed[[i]]
  x <- append(x, paste0("  ", tag, c1$txt), after = c1$last)
}
anchor <- which(x == "## Filed at the 0.81.0 release (2026-11-10)")
x <- append(x, filed, after = anchor - 1L)
writeLines(x, f)
cat("added", length(closed), "CLOSED lines and", length(filed),
    "filed lines\n")

# Lane silent2: add the CLOSED lines and the lane's filed items to
# dev/test-backlog.md. Each CLOSED line goes directly under the item whose
# first line is named, after the item's last line. Run once.
p <- "dev/test-backlog.md"
x <- readLines(p, warn = FALSE)
item_end <- function(first_pat) {
  i <- grep(first_pat, x, fixed = TRUE)
  stopifnot(length(i) == 1L)
  j <- i + 1L
  # an item runs until the next item, heading or blank line
  while (j <= length(x) && nzchar(x[j]) && !startsWith(x[j], "- ") &&
         !startsWith(x[j], "#")) {
    j <- j + 1L
  }
  j - 1L
}
closed <- list(
  list("- **LOW. The degenerate check's \"collapsed\" test is exact-double",
       c("  CLOSED in the round of 2026-11-13 (lane silent2): two thresholds",
         "  within 8 ulps are collapsed, and a threshold with a standard",
         "  error above 100 on the latent scale warns (the six flat moves);",
         "  `test-ordinal-mixture.R`, `dev/silent2-findings.md`.")),
  list("- **MEDIUM. A flat ordinal-mixture threshold with no warning at all**",
       c("  CLOSED in the round of 2026-11-13 (lane silent2): the degenerate",
         "  check reports a threshold whose standard error is above 100 on",
         "  the latent scale (1.75e7 here); `test-ordinal-mixture.R`,",
         "  `dev/silent2-findings.md`.")),
  list("- **MINOR. An addition term whose value has the length of the new rows",
       c("  CLOSED in the round of 2026-11-13 (lane silent2): a value of more",
         "  than one row that reads no column of the new rows is refused by",
         "  `aterms_for_newdata()` and `ce_aterms()`;",
         "  `tests/testthat/test-newdata-aterms.R`.")),
  list("- **MINOR. `anova()`'s same-data check compares a count of rows, not",
       c("  CLOSED in the round of 2026-11-13 (lane silent2): `anova()` compares",
         "  the values each likelihood reads, in row-name order",
         "  (`loglik_rows_key()`); `tests/testthat/test-subset-rate.R`.")),
  list("- **MEDIUM. A prediction rebuilds a group-level term's design under the",
       c("  CLOSED in the round of 2026-11-13 (lane silent2): each group-level",
         "  term stores its contrasts (`re_term_contrasts()`), read by",
         "  `pred_design()`, `mm()` and `re_formula`;",
         "  `tests/testthat/test-predict-newdata.R`.")),
  list("- **MEDIUM. Laplace draws give `conditional_effects(method =",
       c("  CLOSED in the round of 2026-11-13 (lane silent2): the watch reads",
         "  every predicted dpar on the grid, as `posterior_predict()`'s does;",
         "  frmtmb.sample `test-laplace-draws.R`.")),
  list("- **MINOR. Floors outside a support that a reader still sees.** With",
       c("  CLOSED in the round of 2026-11-13 (lane silent2): the ordinal gap",
         "  floor is exact for a reader; the rest is refiled, narrowed, under",
         "  \"From lane silent2\" (`dev/silent2-findings.md`, item 5524)."))
)
# bottom first, so earlier line numbers stay valid
ends <- vapply(closed, function(cl) item_end(cl[[1L]]), 0L)
for (k in order(ends, decreasing = TRUE)) {
  e <- ends[k]
  x <- append(x, closed[[k]][[2L]], after = e)
}
filed <- c(
  "## Filed in the round of 2026-11-13",
  "",
  "### From lane silent2",
  "",
  "- **LOW. `print()` of a draws object lists the parameters in the stored",
  "  order**: `print.frmtmb_draws()` read `x$draws`, whose template order",
  "  pairs a `gp(by = )` term's `sdgp_` and `lscale_` columns by level",
  "  (sdgp_gpxfa lscale_gpxfa sdgp_gpxfb lscale_gpxfb), where",
  "  `variables()` lists brms's order (triage-082 part A, item 1915).",
  "  CLOSED in the round of 2026-11-13 (lane silent2): `print()` reads",
  "  `variables()` order; frmtmb.sample `test-brms-draws-order.R`.",
  "- **MINOR. The extensions' floors outside a support, for an exact",
  "  reader** (what remains of the 0.81.0 item \"Floors outside a support",
  "  that a reader still sees\"). Kind S, where the exact value is `-Inf`",
  "  and a switch on `exact_support()` would give it: frmtmb.eam's",
  "  `lba_race_lpdf()` at a decision time <= 0 (lba and rdm),",
  "  `gng_lpdf()`'s go rows at y <= ndt, the wiener `allow_unreachable`",
  "  branch and the contaminant outside `crange`; frmtmb.survival's",
  "  `sv_floor_pos`. Each needs `exact_support()` in the sampling API and",
  "  a reader test of its own. Kind U, where the exact value is positive",
  "  and below what the formula resolves, cannot be fixed by a switch:",
  "  the lba and rdm densities and rdm's lcdf floored at 1e-300, gddm's",
  "  grid floor, survival's caps and `sv_squeeze()`; nor can",
  "  `log1m_exp_ad()`, whose floor meets an empty window and an",
  "  unresolved one alike (`dev/silent2-findings.md`, item 5524).",
  "- **LOW. The ordinal-mixture threshold standard-error test needs the",
  "  fit-end covariance.** `mixture_ord_thres_se()` reads",
  "  `se_fit_analysis()`, which a fit with random effects builds at the",
  "  fit only when `se_check_at_fit()` says so; otherwise the degenerate",
  "  check runs without it and a flat threshold passes until the first",
  "  `summary()` (whose standard-error warning then names it). Not",
  "  measured how often an ordinal mixture with random effects defers.",
  "- **LOW. A flat direction on an ordinal mixture's `cs()` coefficient is",
  "  not named by the degenerate check.** ordstart's cs_sr_acat seed 4",
  "  (coefficient 325.5) passes it; the standard-error check reports \"not",
  "  available for 4 of 15 parameters\" there, so the fit is not silent,",
  "  but the message does not say the component is degenerate",
  "  (`dev/silent2-out/ord-sweep-base.tsv`).",
  "")
h <- grep("^## Filed at the 0.82.0 release", x)
stopifnot(length(h) == 1L)
x <- append(x, filed, after = h - 1L)
writeLines(x, p)
cat("CLOSED lines after", paste(sort(ends), collapse = " "), "; filed above",
    h, "\n")

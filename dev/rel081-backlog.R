# dev/test-backlog.md at the 0.81.0 release: the three lanes' sections
# "Filed in the round of 2026-11-09" (merged by --union) become one
# section "Filed at the 0.81.0 release", items kept as the lanes wrote
# them, then the reviews' and the consolidation's items left for later.
#
#   Rscript dev/rel081-backlog.R
REL <- "C:/Users/adf44/source/r/frmtmb-wt-rel081"
p <- file.path(REL, "dev/test-backlog.md")
x <- sub("\r$", "", readLines(p, warn = FALSE))
h <- which(x == "## Filed in the round of 2026-11-09")
stopifnot(length(h) == 1L)
nxt <- which(x == "## Filed at the 0.80.0 release (2026-11-06)")
stopifnot(length(nxt) == 1L, nxt > h)
body <- x[(h + 1L):(nxt - 1L)]
# a blank line before each lane's sub-heading (the union merge joined them)
out <- character()
for (l in body) {
  if (grepl("^### ", l) && length(out) && nzchar(out[length(out)])) out <- c(out, "")
  out <- c(out, l)
}
while (length(out) && !nzchar(out[length(out)])) out <- out[-length(out)]
while (length(out) && !nzchar(out[1L])) out <- out[-1L]
stopifnot(sum(grepl("^### From lane ", out)) == 3L)
intro <- c("## Filed at the 0.81.0 release (2026-11-10)", "",
  "Gathered from the lanes' own sections of the round of 2026-11-09",
  "(compat, code8, hollow; translator and support filed none here), then",
  "the reviews' and the consolidation's items left for later. Lane",
  "records: `dev/<lane>-findings.md`, reviews",
  "`dev/reviews/2026-11-10-<lane>.md` (translator:",
  "`dev/reviews/2026-11-09-translator.md`).", "")
more <- c("",
"### From the reviews and the lanes' records, left for later", "",
"- **MINOR. The live-row instrument does not count a conditional row's",
"  refused half.** `dev/equiv-flow.R` scores \"refusal correct\" when both",
"  packages refuse and the live row is \"conditional\". A conditional row",
"  whose spec was meant to meet the working half, where both packages",
"  refuse, then passes without a flag; the reason text names it, the",
"  summary does not count it. A generated count line would make such",
"  pairs readable (review of lane compat, round 2, m2 cost note).", "",
"- **MINOR. Floors outside a support that a reader still sees.** With",
"  `row_lpdf(exact = TRUE)` in place, the eam floors",
"  (`wiener-family.R:986`, `:1059`, `:1070`, `:1108`, the contaminant,",
"  `wiener-gng.R:642`, `lba.R`, `rdm.R`, `gddm.R:1527`), frmtmb.survival's",
"  `sv_floor_pos` and caps, the ordinal mixture's gap floor",
"  (`.gap_floor`) and `log1m_exp_ad()`'s window floor can be made exact",
"  for a reader through `exact_support()` without touching their fits",
"  (`dev/support-findings.md`, the survey table, rows \"not changed\").", "",
"- **Translator gaps the survey still has** (`dev/translator-findings.md`,",
"  section 9): an mv mixture's `theta_<resp>` (p005), `mo()` outside mu",
"  (`simo_<dpar>`, p030), an unpredicted `nu` in a multivariate student",
"  model (p068), `tmp_xi` with a predicted xi, an ordinal mixture with",
"  shared thresholds (`order = \"mu\"`), and `gr(by = )`'s sd matrix (the",
"  audit's by-group branch has no pair to read). The first three agree at",
"  L4 by the optimizing route.", "",
"- **The code-8 message on ubuntu needs CI.** `test-se-reparam.R:554`",
"  warned \"false convergence (8)\" at the maximum on every ubuntu job of",
"  0.80.0; on this machine the fit stops at code 0 under reference BLAS",
"  and OpenBLAS 0.3.26 alike (`dev/code8-aexp.R`), so the new message",
"  there is not seen here (`dev/code8-findings.md`, OpenBLAS runs).", "",
"- **MINOR. `loo()` and `psis()` of draws stop on a `-Inf` cell.** With",
"  the exact reader, a draw that puts a `gen_extreme_value()` row outside",
"  its support gives `-Inf` in `log_lik()`. `loo()` and `psis()` then",
"  stop with the loo package's \"All input values must be finite or",
"  -Inf\": the importance log ratio is `-log_lik()`, `+Inf` in that cell,",
"  and loo 2.10.1's `validate_ll()` refuses `+Inf` (it accepts `-Inf`).",
"  `waic()` does not stop: its `elpd_waic` is `NaN`. brms 2.23.0's",
"  `dgen_extreme_value()` gives `NaN` (\"NaNs produced\") in the same",
"  cell, and `loo()` and `psis()` of that matrix stop with \"NAs not",
"  allowed in input\" (`waic()` too, by the same check). Seen at the",
"  consolidation on a constructed draws object (`dev/rel081-loo-exact.R`,",
"  `dev/rel081-loo-inf.R`). A sampled scalar `xi` does not reach it: the",
"  sampler reads it through brms's `tmp_xi` chart, which keeps every",
"  fitted row inside the support. Where that chart is blocked (a",
"  predicted `xi`, a `cens()` or `trunc()` term, a second GEV response,",
"  a mixture component) the sampler reads the fitted objective, whose",
"  barrier is finite just outside a row's support, so a draw can reach",
"  it there. brms rejects such a draw in its sampler.")
y <- c(x[seq_len(h - 1L)], intro, out, more, "", x[nxt:length(x)])
con <- file(p, open = "wb")
writeLines(y, con, sep = "\r\n")
close(con)
cat("backlog: section of", sum(grepl("^- ", c(out, more))), "items\n")

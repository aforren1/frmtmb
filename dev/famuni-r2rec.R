f <- "C:/Users/adf44/source/r/frmtmb-wt-famuni/dev/famuni-findings.md"
x <- paste(readLines(f), collapse = "\n")
ed <- function(old, new) {
  if (!grepl(old, x, fixed = TRUE)) stop("not found: ", substr(old, 1, 50))
  x <<- sub(old, new, x, fixed = TRUE)
}
ed("- Not covered, and warned by `frm_sample()`: a GEV response with
  `cens()` or `trunc()` (the chart carries no CDF), and a model with
  two GEV responses (one REPORTed name). These keep the log1p scale.",
"- Not covered, and warned by `frm_sample()`: a GEV response with
  `cens()` or `trunc()` (the chart carries no CDF), a model with two
  GEV responses (one REPORTed name), and the GEV components of a
  `mixture()` (the chart maps a plain response's `xi` only; the
  mixture's warning was added in fix round 2, section 13). These keep
  the log1p scale.")
ed("   usage for `dpar_has_predictor()`. A `@usage` line fixed it, a change
   to the Rd only. Core was reinstalled from it (the last line of
   `final-install-sample2.txt`, which lists the installs as they stand at
   the end: eam 19:48, sample 20:07, core from the Rd fix). No test reads
   a package Rd file (grep),
   so the suite runs above stand for it; four core files were rerun on
   the reinstall as a check (`final-core-*.txt`):",
"   usage for `dpar_has_predictor()`. A `@usage` line fixed it. Core was
   reinstalled from it (the last line of `final-install-sample2.txt`,
   which lists the installs as they stand at the end: eam 19:48, sample
   20:07, core from that fix). That the 20:39 source change was the
   `@usage` line alone is NOT verifiable by diff: no pre-edit copy
   survives, and roxygen rewrote every Rd file and NAMESPACE at 20:39:46.
   It is covered by reruns on the 20:46 build instead: the final core
   R CMD check (tests OK, 1 NOTE), the review's three extension test
   reruns (round 2), and these four core files (`final-core-*.txt`):")
writeLines(strsplit(x, "\n")[[1]], f)

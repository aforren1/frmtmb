# Fill PASTE-TESTS of dev/famuni-findings.md: prose with the RESULT lines
# and suite totals pasted from the logs.
setwd("C:/Users/adf44/source/r/frmtmb-wt-famuni")
o <- "dev/famuni-out/"
res <- function(f) grep("^RESULT", readLines(file.path(o, f)), value = TRUE)
head_of <- function(d) readLines(file.path(d, "SUMMARY.txt"))[1:5]
perpkg <- function(d) {
  x <- readLines(file.path(d, "SUMMARY.txt"))
  x <- x[(which(x == "per file (tests pass fail error skip warn):") + 1):
           length(x)]
  p <- sub("^  ([^ ]+)--.*$", "\\1", x)
  m <- do.call(rbind, lapply(strsplit(trimws(x), " +"), function(v) {
    as.numeric(v[-1])
  }))
  agg <- rowsum(m, p)
  sprintf("%-16s files %3d tests %4d pass %5d fail %d error %d skip %3d warn %d",
          rownames(agg), as.integer(table(p)[rownames(agg)]), agg[, 1],
          agg[, 2], agg[, 3], agg[, 4], agg[, 5], agg[, 6])
}
ck <- function(p) {
  f <- file.path("dev/famuni-check", p, "check.log")
  if (!file.exists(f)) return(paste(p, ": not run"))
  con <- file(f, encoding = "UTF-16LE")
  x <- readLines(con, warn = FALSE)
  close(con)
  x <- gsub("[^ -~]", "", x)
  paste0(p, ": ", paste(grep("^Status|NOTE$|WARNING$|ERROR$", x,
                             value = TRUE), collapse = " | "))
}
txt <- c(
"`dev/famuni-test1.R` runs one file per R process against the lane",
"library and prints the `lib:` line of the frmtmb it loaded;",
"`dev/famuni-suite.sh`, `dev/famuni-gated.sh` and",
"`dev/famuni-sample-suite.sh` drive it, and `dev/famuni-summary.R`",
"writes each `SUMMARY.txt`.",
"",
"The new test file, seen to fail on rellib-r8 first (0.70.0; the",
"behavioral form is the first test, `frm(family = \"frechet\")` refused",
"as \"frechet is not a supported family\", and the wiener message), then",
"on the lane build, and under OpenBLAS 0.3.26 with LAPACK",
"(`dev/ciharden-openblas.sh 0.3.26 lapack`):",
"",
"```text",
paste("rellib-r8:", res("t-families-extreme-base.txt")),
paste("lane:     ", res("t-families-extreme.txt")),
paste("OpenBLAS: ", res("ob-test-families-extreme.R.txt")),
paste("OpenBLAS: ", res("ob-test-simulate-density.R.txt")),
"```",
"",
"Within that file, two tests caught defects of the lane's own first",
"versions, recorded rather than lost: the gradient at `xi = 0` was",
"off by up to 8.2e-08 (the penalty written as `zc - z` summed -1e8 and",
"+1e8 into the gradient beside its O(1) terms;",
"`dev/famuni-gradprobe3.R`), and the GEV CDF was `NaN` at 16 of 60",
"points beyond the upper end (`dev/famuni-cdfprobe.R`).",
"",
"The full suite, core and the seven extensions, ungated, one file per",
"process (`dev/famuni-suite-full/SUMMARY.txt`):",
"",
"```text", head_of("dev/famuni-suite-full"), perpkg("dev/famuni-suite-full"),
"```",
"",
"The one failure is the frmtmb.sample test this lane added for the",
"scalar-`xi` note, run before the fix it pins was installed (the note",
"still fired with the user's own `xi` prior): seen to fail, then fixed.",
"Against the fixed build, every frmtmb.sample file, gated",
"(`dev/famuni-suite-sample/SUMMARY.txt`):",
"",
"```text", head_of("dev/famuni-suite-sample"), "```",
"",
"The gated tier, the 29 files gated on `FRMTMB_BRMS_FIT_TESTS`",
"(`dev/famuni-suite-gated/SUMMARY.txt`), whose only skips were",
"`test-drmtmb-agreement.R`'s 13, gated on `FRMTMB_DRMTMB_FIT_TESTS`;",
"that file run with it:",
"",
"```text", head_of("dev/famuni-suite-gated"),
paste("drmTMB gate:", res("t-drmtmb-gated.txt")), "```",
"",
"A trap worth recording: `FRMTMB_DRMTMB_FIT_TESTS=true timeout ...",
"Rscript` lost the variable (13 skips again); exported in the shell and",
"run without the Rtools `timeout`, it reached R.",
"",
"**The ported brms suite.** No row flips. The gated tier ran all 15",
"ported files with no failure, so no row went to a stale \"now HOLDS\";",
"`test-brms-suite-families.R` passes 84 assertions in core and in",
"frmtmb.sample, as at 0.70.0 (76 of 84 rows pass, 8 cannot transfer).",
"brms's own suite exercises these three families only in",
"tests.stancode.R, tests.log_lik.R, tests.posterior_epred.R,",
"tests.posterior_predict.R and tests.distributions.R, none of which is",
"ported.",
"",
"**R CMD check --as-cran** (built with vignettes, in",
"`dev/famuni-check/<pkg>/`):",
"",
"```text", ck("frmtmb"), ck("frmtmb.eam"), ck("frmtmb.sample"), "```",
"",
"The NOTE is the HTML manual's math-rendering check (\"package 'V8'",
"unavailable\"), environmental on this machine.",
"",
"The first core run (`dev/famuni-check/frmtmb/check-run1.log`) ended",
"`Status: 2 WARNINGs, 1 NOTE`: `link_sigma` and `link_shape` were not",
"documented on `?frmtmb-extreme-families` (`@inheritParams` did not",
"reach them), and the new test called `evd::dgev`, an undeclared",
"dependency. Both fixed (explicit `@param`; the test reads evd by name",
"after `skip_if_not_installed()`), and the run above is the second;",
"its tests passed in the first run as well."
)
f <- "dev/famuni-findings.md"
md <- readLines(f)
b0 <- which(md == "<!-- BEGIN generated: TESTS -->")
b1 <- which(md == "<!-- END generated: TESTS -->")
new <- c("<!-- BEGIN generated: TESTS -->", txt, "<!-- END generated: TESTS -->")
i <- which(md == "PASTE-TESTS")
if (length(b0) == 1L) {
  md <- c(md[seq_len(b0 - 1L)], new, md[-seq_len(b1)])
} else {
  md <- c(md[seq_len(i - 1L)], new, md[-seq_len(i)])
}
writeLines(md, f)
cat("ok\n")

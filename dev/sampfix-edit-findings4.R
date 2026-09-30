# Lane sampfix nits round: final counts and the check-found helper.
f <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix/dev/sampfix-findings.md"
s <- paste(gsub("\r", "", readLines(f)), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(s, gregexpr(old, s, fixed = TRUE)))
  if (n != 1L) stop("found ", n, ":\n", old)
  s <<- sub(old, new, s, fixed = TRUE)
}
rep1("round, 40 files, 332 tests, 2323 pass, 0 fail, 0 error, 0 warn, 1 skip",
     "round, 40 files, 333 tests, 2324 pass, 0 fail, 0 error, 0 warn, 1 skip")
rep1("| `test-laplace-draws.R` | 14 pass, 18 fail, 9 error | 89 pass |",
     "| `test-laplace-draws.R` | 14 pass, 18 fail, 10 error | 90 pass |")
rep1("`R CMD check --as-cran` (before the nits round, which changed
frmtmb.sample's code and docs; not rerun), built with vignettes inside",
"`R CMD check --as-cran` (frmtmb before the nits round, which did not
touch core; frmtmb.sample rerun after it), built with vignettes inside")
rep1("- **S4.** The corrections above",
"- **Found by R CMD check in this round:** the rewrite of the probe block
  dropped `is_na_re_form()`, which `conditional_effects()` calls at
  `re_formula = NULL` with a group-level term. The probe then failed at
  both fills, took that for an error unrelated to the fill, and let the
  call through, and a smooth's curves came back `NA` without a word.
  Restored, and a test (`test-laplace-draws.R`, a curve at
  `re_formula = NULL` that reads a smooth) was seen to fail on the build
  without it (`dev/sampfix-log/nits-laplace-noisna.txt`). The residual
  risk stays: a probe whose own code errors at both fills passes the
  call on, and only the per-draw check stands behind it, which cannot
  see a result that is non-finite at every draw.
- **S4.** The corrections above")
writeLines(s, f)
cat("ok\n")

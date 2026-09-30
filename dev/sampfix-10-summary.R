# Lane sampfix, script 10: the numbers dev/sampfix-findings.md quotes,
# read from the logs rather than typed.
#   Rscript dev/sampfix-10-summary.R > dev/sampfix-log/10-summary.txt
LOG <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix/dev/sampfix-log"
rd <- function(f) readLines(file.path(LOG, f), warn = FALSE)
pick <- function(f, pat) {
  x <- rd(f)
  x[grepl(pat, x)]
}
say <- function(...) cat(..., "\n", sep = "")

say("== 01 repro, reference build (rellib-r3)")
for (l in pick("01-ref.txt", "non-finite|ERROR: invalid|no more scalars|variables\\(ds\\)|nchains\\]|VarCorr\\]|summary\\]")) say(l)
say("== 01 repro, lane build")
for (l in pick("01-lane.txt", "non-finite|laplace = TRUE with|y ~ x, laplace|y ~ 0 \\+ x|variables\\(ds\\)|nchains\\]|VarCorr\\]|summary\\]")) say(substr(l, 1, 120))
say("== 02 one-parameter model")
for (f in c("02-ref.txt", "02-lane.txt")) {
  for (l in pick(f, "^\\[|ARM")) say(l)
}
say("== 03 sweep: calls whose result changed, ref -> lane")
a <- rd("03-ref.txt")
b <- rd("03-lane.txt")
strip <- function(x) sub("  \\[[0-9]+ warn.*$", "", x)
arm <- ""
n_changed <- 0L
for (i in seq_along(a)) {
  if (grepl("^====", a[i])) arm <- sub("^==== ([a-z]+):.*", "\\1", a[i])
  if (!grepl("^  [a-zA-Z]", a[i])) next
  if (!identical(strip(a[i]), strip(b[i]))) n_changed <- n_changed + 1L
}
say("calls: ", sum(grepl("^  [a-zA-Z]", a)), ", changed: ", n_changed)
say("ref NONFINITE or bare-error lines on the laplace object:")
lap_a <- a[seq_len(which(grepl("^==== nul", a)) - 1L)]
for (l in lap_a[grepl("NONFINITE|missing values|NAs not allowed|STATS is longer", lap_a)]) {
  say(substr(l, 1, 110))
}
nul_a <- a[which(grepl("^==== nul", a)):length(a)]
say("ref calls on the stanfit = NULL object that died on `@`: ",
    sum(grepl("for `@` applied to an object of class", nul_a, fixed = TRUE)))
nul_b <- b[which(grepl("^==== nul", b)):length(b)]
say("lane calls on the stanfit = NULL object that died on `@`: ",
    sum(grepl("for `@` applied to an object of class", nul_b, fixed = TRUE)))
say("== 04 laplace validation")
for (l in rd("04-lane.txt")) if (nzchar(l)) say(l)
say("== 05 pp_check against brms")
for (l in pick("05-lane.txt", "log_lik max|nchains|VarCorr sd|^  loo_|^-- draws")) say(l)
say("== 06 probes")
for (f in c("06-ref.txt", "06-lane.txt")) {
  for (l in pick(f, "ARM|equals|columns: b_Intercept b_x sigma$|ERROR: frm_sample|non-finite")) {
    say(substr(l, 1, 120))
  }
}
say("== 07 ordinal names, lane against reference")
for (l in pick("07-lane.txt", "names TRUE|names FALSE")) say(l)
say("== 08 conditional_effects at an observed level (defect, not fixed)")
for (l in rd("08-lane.txt")) if (grepl("full|laplace|fit|ranef", l)) say(l)
say("== new test files, reference build (seen to fail)")
for (f in c("t-laplace-ref.txt", "t-ppcheck-ref.txt", "t-nosf-ref.txt",
            "t-ordnames-ref.txt")) {
  for (l in pick(f, "^RESULT")) say(l)
}

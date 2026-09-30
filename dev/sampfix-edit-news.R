# Lane sampfix nits round: NEWS wording (S2, M5).
f <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix/extensions/frmtmb.sample/NEWS.md"
s <- paste(gsub("\r", "", readLines(f)), collapse = "\n")
rep1 <- function(old, new) {
  stopifnot(lengths(regmatches(s, gregexpr(old, s, fixed = TRUE))) == 1L)
  s <<- sub(old, new, s, fixed = TRUE)
}
rep1("  smooth's coefficients, `mi()` values) is refused by name, as
  `log_lik()` already refused; one that does not, such as",
"  smooth's coefficients, `mi()` values) is refused, naming the function
  you called, as `log_lik()` already refused; one that does not, such as")
rep1("  `tau_raw_1` from the draws must select `b_Intercept[1]` now. The
  renaming needs the frmtmb that declares each ordinal family's inverse
  threshold map; with an older frmtmb the draws keep the internal names,
  as before.",
"  `tau_raw_1` from the draws must select `b_Intercept[1]` now. The
  names come from each ordinal family's inverse threshold map, which the
  frmtmb this release requires declares.")
writeLines(s, f)
cat("ok\n")

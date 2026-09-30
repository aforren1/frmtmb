# Lane sampfix nits round: suite counts after the round.
f <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix/dev/sampfix-findings.md"
s <- paste(gsub("\r", "", readLines(f)), collapse = "\n")
old <- "`test-fuzz.R`, gated on `FRMTMB_FUZZ`); frmtmb.sample, 40 files, 2275
pass, 0 fail, 0 error, 0 warn, 1 skip (`test-scale.R`, gated);"
new <- "`test-fuzz.R`, gated on `FRMTMB_FUZZ`); frmtmb.sample after the nits
round, 40 files, 332 tests, 2323 pass, 0 fail, 0 error, 0 warn, 1 skip
(`test-scale.R`, gated);"
stopifnot(grepl(old, s, fixed = TRUE))
s <- sub(old, new, s, fixed = TRUE)
old2 <- "`R CMD check --as-cran`, built with vignettes inside"
new2 <- "`R CMD check --as-cran` (before the nits round, which changed
frmtmb.sample's code and docs; not rerun), built with vignettes inside"
stopifnot(grepl(old2, s, fixed = TRUE))
s <- sub(old2, new2, s, fixed = TRUE)
writeLines(s, f)
cat("ok\n")

f <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix/dev/sampfix-findings.md"
s <- paste(gsub("\r", "", readLines(f)), collapse = "\n")
old <- "(`test-scale.R`, gated);"
new <- "(`test-scale.R`, gated). The last change of the round (the probe's
re-raise going through `frm_stop()`, which `test-conditions-census.R`
requires) was followed by a rerun of the 8 files it can reach, not of
all 40;"
stopifnot(lengths(regmatches(s, gregexpr(old, s, fixed = TRUE))) == 1L)
writeLines(sub(old, new, s, fixed = TRUE), f)
cat("ok\n")

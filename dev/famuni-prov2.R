f <- "C:/Users/adf44/source/r/frmtmb-wt-famuni/dev/famuni-findings.md"
x <- paste(readLines(f), collapse = "\n")
old <- "4. R CMD check --as-cran of core and sample ran on the final source;"
new <- "4. The core check of step 4's first pass (`check-run3.log`, 2 NOTEs)
   found `?frmtmb-gev-xi-chart` documenting `spec` and `lp` without a
   usage for `dpar_has_predictor()`. A `@usage` line fixed it, a change
   to the Rd only. Core was reinstalled from it (the last line of
   `final-install-sample2.txt`); no test reads an Rd file, so the suite
   runs above stand for it.
5. R CMD check --as-cran of core and sample ran on the final source;"
stopifnot(grepl(old, x, fixed = TRUE))
x <- sub(old, new, x, fixed = TRUE)
writeLines(strsplit(x, "\n")[[1]], f)

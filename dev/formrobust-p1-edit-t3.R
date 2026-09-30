# Punch round 1: comment of the lme4#682 test. Record of the edit.
p <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/tests/testthat/test-open-issues.R"
x <- paste(readLines(p), collapse = "\n")
old <- "  # brms codes any two values 0 and 1, and frmtmb does too since lane
  # formrobust, but not fractions, which are a proportion response"
new <- "  # brms codes any two values 0 and 1, and frmtmb does too since lane
  # formrobust; two values inside (0, 1) are then fitted, not silently:
  # they warn that they look like a proportion"
stopifnot(grepl(old, x, fixed = TRUE))
x <- sub(old, new, x, fixed = TRUE)
x <- sub("non-integer responses are rejected, not silently fit",
         "non-integer responses are not silently fit", x, fixed = TRUE)
con <- file(p, "wb"); writeLines(strsplit(x, "\n")[[1]], con, sep = "\r\n")
close(con)

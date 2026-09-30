# Round 3: the function-call guard case uses log(wt) + 1 against the
# precomputed column. Record of the edit.
p <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/tests/testthat/test-aterm-expr.R"
x <- readLines(p)
i <- grep("weights(exp(log(wt)))", x, fixed = TRUE)
stopifnot(length(i) == 1L)
x[i] <- "  expect_identical(logLik(frm(bf(y | weights(log(wt) + 1) ~ x),"
x[i + 2L] <- "                   logLik(frm(bf(y | weights(lw) ~ x), data = d)))"
writeLines(x, p)

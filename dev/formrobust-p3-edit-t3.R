# Round 3: a full-length vector outside the data is refused too.
p <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/tests/testthat/test-aterm-expr.R"
x <- readLines(p)
i <- which(x == "  # the guard-absent cases: the same constant as a data column")
stopifnot(length(i) == 1L)
x <- append(x, c(
"  # a full-length vector outside the data too, which 0.66.0 read",
"  w_out <- ae_data$wt",
"  expect_error(frm(bf(y | weights(w_out) ~ x), data = ae_data),",
"               \"reads `w_out`, which is not a column of `data`\",",
"               fixed = TRUE)"), i - 1L)
writeLines(x, p)

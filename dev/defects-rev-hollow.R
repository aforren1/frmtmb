# Reviewer of lane defects: does brms_hollow_empty() fire on every
# vacuous all()/any() form, and never on an empty-but-correct one?
# Each case is run through brms_port_run(), the harness's own entry.
.libPaths(c("C:/Users/adf44/source/r/wt-defects-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(testthat))
source("tests/testthat/helper-brms-suite.R")
e <- new.env()
e$sdata <- list(X = matrix(1, 3, 1, dimnames = list(NULL, "Intercept")))
e$empty_ok <- character(0)   # e.g. "no warnings were raised"
e$x <- c(1, 2)
cases <- list(
  # vacuous: a NULL operand through several spellings
  a1 = quote(expect_true(all(sdata$idxl_y_x_1 %in% 9:5))),
  a2 = quote(expect_true(all(integer(0) == x))),
  a3 = quote(expect_true(isTRUE(all(sdata$idxl %in% 1)))),
  a4 = quote(expect_true(all(sdata$idxl %in% 1, na.rm = TRUE))),
  a5 = quote(expect_true(!any(sdata$idxl > 3))),
  a6 = quote(expect_false(!all(sdata$idxl > 3))),
  a7 = quote(expect_equal(all(sdata$idxl > 3), TRUE)),
  a8 = quote(expect_identical(any(sdata$idxl > 3), FALSE)),
  a9 = quote(testthat::expect_true(all(sdata$idxl > 3))),
  a10 = quote(expect_true(all(sdata$idxl > 3) && TRUE)),
  a11 = quote(expect_true(!anyNA(sdata$idxl))),
  a12 = quote(expect_true(all(dim(sdata$idxl) == c(3, 1)))),
  # empty but correct: "no element matches" is the claim itself
  b1 = quote(expect_false(any(grepl("error", empty_ok)))),
  b2 = quote(expect_true(all(nchar(empty_ok) > 0))),
  # controls: non-empty operand, must not be flagged
  c1 = quote(expect_true(all(x > 0))),
  c2 = quote(expect_false(any(x > 5)))
)
for (k in names(cases)) {
  r <- brms_port_run(cases[[k]], e)
  cat(sprintf("%-4s held=%-5s vacuous=%-5s raw=%-5s %s | %s\n", k, r$held,
              r$vacuous, r$raw_held, deparse1(cases[[k]]), r$msg))
}

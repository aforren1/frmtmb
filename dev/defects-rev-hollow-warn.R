# Reviewer of lane defects: does brms_hollow_empty()'s second evaluation
# of an operand let a warning escape, and is the operand evaluated again?
.libPaths(c("C:/Users/adf44/source/r/wt-defects-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(testthat))
source("tests/testthat/helper-brms-suite.R")
e <- new.env()
e$n <- 0
e$f <- function() { n <<- n + 1; warning("operand warned"); c(1, 2) }
environment(e$f) <- e
w <- 0
r <- withCallingHandlers(
  brms_port_run(quote(expect_true(all(f() > 0))), e),
  warning = function(cnd) { w <<- w + 1; invokeRestart("muffleWarning") })
cat("held =", r$held, " warnings escaping brms_port_run =", w,
    " operand evaluations =", e$n, "\n")

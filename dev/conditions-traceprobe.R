# Lane wt-conditions: does a trace() on testthat's expect_error() see
# the condition a passing expectation caught? The sweep depends on it.
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(testthat)
seen <- list()
for (w in list(asNamespace("testthat"),
              as.environment("package:testthat"))) {
  suppressMessages(trace("expect_error", where = w,
  exit = quote(seen[[length(seen) + 1L]] <<- returnValue()),
  print = FALSE))
}
f <- function() stop("boom", call. = FALSE)
expect_error(f(), "boom")
expect_error(stop("x"))
expect_error(1, NA)
local(expect_error(stop("inner"), "inner"))
test_that("inside a test", expect_error(f(), "boom"))
cat("recorded:", length(seen), "\n")
for (s in seen) cat(paste(class(s), collapse = "/"), "\n")
cat("the attached binding is traced:",
    inherits(get("expect_error", "package:testthat"), "functionWithTrace"),
    "\n")

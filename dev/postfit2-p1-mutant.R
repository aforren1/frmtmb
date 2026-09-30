# Punch round 1, B7: the reuse key's `kept` field is pinned by a test.
# Mutant `keydrop` removes `kept = kept` from ce_boot_key()'s key; the
# test file must then fail, and must pass without the mutant.
#   Rscript dev/postfit2-p1-mutant.R <keydrop|none>
mutant <- commandArgs(TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/wt-postfit2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages({
  library(testthat)
  library(frmtmb)
})
drop_arg <- function(e, argname, count) {
  if (is.call(e)) {
    if (identical(e[[1]], as.name("list")) && argname %in% names(e)) {
      count$n <- count$n + 1L
      e <- e[names(e) != argname]
    }
    for (i in seq_along(e)) {
      if (!is.null(e[[i]]) && !identical(e[[i]], quote(expr = ))) {
        e[[i]] <- drop_arg(e[[i]], argname, count)
      }
    }
  }
  e
}
if (identical(mutant, "keydrop")) {
  ns <- asNamespace("frmtmb")
  fn <- get("ce_boot_key", envir = ns)
  count <- new.env()
  count$n <- 0L
  body(fn) <- drop_arg(body(fn), "kept", count)
  # a mutant that changed nothing would pass for the wrong reason
  stopifnot(count$n == 1L)
  environment(fn) <- ns
  utils::assignInNamespace("ce_boot_key", fn, ns = "frmtmb")
}
cat("mutant:", mutant, "\n")
r <- test_file("C:/Users/adf44/source/r/frmtmb-wt-postfit2/tests/testthat/test-ce-levels.R",
               package = "frmtmb", env = test_env("frmtmb"),
               reporter = ListReporter$new(), stop_on_failure = FALSE)
df <- as.data.frame(r)
cat(sprintf("RESULT keydrop-run %s: failed=%d error=%d passed=%d\n", mutant,
            sum(df$failed), sum(df$error), sum(df$passed)))
print(df[df$failed > 0 | df$error, c("test", "failed", "error")])

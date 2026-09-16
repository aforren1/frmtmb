## FALSIFICATION: the structural guard is not total. Its detector asks
## whether the literal three dots appear in deparse(body(fn)), so a
## method whose only "..." is inside a STRING is scored as a dots user.
## print.frmtmb_par_template prints "... N more" and never reads its
## dots, so it still swallows every argument it is given, and both
## test-arg-refusal.R's last block and dev/argspell-report.R's
## "METHODS STILL SWALLOWING THEIR DOTS: 0" miss it.
.libPaths(c("C:/Users/adf44/source/r/asrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
set.seed(3)
dd <- data.frame(x = rnorm(60))
dd$y <- 1 + 0.6 * dd$x + rnorm(60)
fit <- frm(bf(y ~ x) + gaussian(), data = dd)
pt <- par_template(fit)
cat("class:", paste(class(pt), collapse = "/"), "\n")
cat("formals:",
    paste(names(formals(getFromNamespace("print.frmtmb_par_template",
                                         "frmtmb"))), collapse = ", "), "\n")
o1 <- utils::capture.output(print(pt))
r <- tryCatch({
  o2 <- utils::capture.output(print(pt, nosucharg = 1, re.form = NA))
  list(ok = TRUE, same = identical(o1, o2))
}, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
if (r$ok) {
  cat("SWALLOWED: print(pt, nosucharg = 1, re.form = NA) returned ",
      "with no error; output identical to print(pt): ", r$same, "\n",
      sep = "")
} else {
  cat("refused:", r$msg, "\n")
}
# the line that fools the detector
src <- deparse(body(getFromNamespace("print.frmtmb_par_template", "frmtmb")))
cat("\nlines of the body matching the literal three dots:\n")
cat(paste0("  ", grep("...", src, fixed = TRUE, value = TRUE)), sep = "\n")
cat("\nany `...` SYMBOL in the parse tree: ",
    "..." %in% all.names(body(getFromNamespace(
      "print.frmtmb_par_template", "frmtmb"))), "\n")

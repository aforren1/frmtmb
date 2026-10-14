# Lane optstart: run one test file on the lane build with one of the
# lane's changes taken out, to tell which change moves it.
#   Rscript dev/optstart-toggle.R <lib> <pkg> <file> <noretake|nostart|norestart>
# noretake: fd_hessian() without the retake (tau_band = Inf)
# nostart: start_ls_plan() returns nothing (0.70.0's start)
# norestart: restart_improves() as 0.70.0 (objective no higher)
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(testthat); library(a[2], character.only = TRUE)})
ns <- asNamespace("frmtmb")
patch <- function(name, f) {
  environment(f) <- ns
  utils::assignInNamespace(name, f, "frmtmb")
}
fd0 <- get("fd_hessian", ns)
switch(a[4],
  noretake = patch("fd_hessian", function(par, gr, tau_band = Inf,
                                          tau_step = 1e-4) {
    fd0(par, gr, tau_band = Inf, tau_step = tau_step)
  }),
  nostart = patch("start_ls_plan", function(frame, start = NULL,
                                            prior_entries = NULL) list()),
  norestart = patch("restart_improves", function(opt, opt2, optCtrl = NULL) {
    isTRUE(opt2$objective <= opt$objective)
  }),
  stop("unknown toggle"))
cat("lib:", dirname(find.package("frmtmb")), "toggle", a[4], "\n")
res <- testthat::test_file(a[3], package = a[2],
                           env = testthat::test_env(a[2]),
                           reporter = "silent")
r <- as.data.frame(res)
for (t in res) for (x in t$results) {
  if (inherits(x, c("expectation_failure", "expectation_error"))) {
    cat("  FAIL [", t$test, "]", substr(conditionMessage(x), 1, 160), "\n")
  }
}
cat(sprintf("RESULT pass=%d fail=%d error=%d skip=%d warn=%d\n",
            sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped),
            sum(r$warning)))

# Lane wt-conditions: can user input reach the stopifnot() calls in
# R/priors.R? Run on the lane library.
.libPaths(c("C:/Users/adf44/source/r/conditions-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
for (p in c("normal(0)", "student_t(3, 0)", "exponential(-1)",
            "gamma(1)", "lkj()")) {
  e <- tryCatch({ set_prior(p); d <- data.frame(y = rnorm(20), x = rnorm(20))
    frm(y ~ x, data = d, prior = set_prior(p, class = "b")); NULL },
    error = identity)
  cat(sprintf("%-18s %s | %s\n", p,
              if (is.null(e)) "no error" else paste(class(e), collapse = "/"),
              if (is.null(e)) "" else substr(conditionMessage(e), 1, 70)))
}
e <- tryCatch(frmtmb.latent::hmm_starts(1), error = identity)
cat("hmm_starts(1):", class(e), "|", conditionMessage(e), "\n")

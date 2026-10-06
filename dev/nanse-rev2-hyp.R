# Reviewer, punch round 1 (B2): hypothesis() and emmeans() on the
# spread ridge y ~ a + b, a ~ 0 + f (k = 10, 5 rows each). a_f1 + b is
# identified with ML SE sigma_ML / sqrt(5); a_f1 alone is not.
#   Rscript dev/nanse-rev2-hyp.R lane|merge
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-nanse-lib",
           "C:/Users/adf44/source/r/rellib-r5"),
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
set.seed(1)
k <- 10
f <- factor(rep(seq_len(k), each = 5))
d <- data.frame(f = f, y = rnorm(k)[f] + rnorm(5 * k, 0, 0.5))
fit <- suppressWarnings(frm(bf(y ~ a + b, a ~ 0 + f, b ~ 0 + 1, nl = TRUE),
                            data = d))
sig <- exp(fit$opt$par[length(fit$opt$par)])
cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  })
  list(v = v, w = w)
}
nm <- variables(fit)
cat("variables:", paste(head(nm, 4), collapse = " "), "...\n")
for (h in c("b_a_f1 + b_b_Intercept = 0", "b_a_f1 = 0",
            "b_a_f1 - b_a_f2 = 0")) {
  r <- cap(hypothesis(fit, h, class = NULL))
  cat(sprintf("%-28s Est.Error %s (sigma_ML/sqrt(5) = %.6g) warnings %d\n",
              h, format(r$v$hypothesis$Est.Error, digits = 7), sig / sqrt(5),
              length(r$w)))
}

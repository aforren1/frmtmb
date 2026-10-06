# Punch round 1, B2: emmeans() and hypothesis() on the spread ridge
# (k = 10): a grid row along a lost direction gets NaN and one warning;
# the identified mu = a_k + b keeps its SE where the route allows it.
#   Rscript dev/nanse-p1-emm.R [lib]
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/wt-nanse-lib"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r5",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
set.seed(1)
k <- 10
f <- factor(rep(seq_len(k), each = 5))
dd <- data.frame(f = f, y = rnorm(k)[f] + rnorm(k * 5, 0, 0.5))
fit <- suppressWarnings(frm(bf(y ~ a + b, a ~ 0 + f, b ~ 1, nl = TRUE),
                            data = dd))
show <- function(lab, expr) {
  w <- character()
  v <- withCallingHandlers(tryCatch(expr, error = function(e) {
    paste("ERROR:", conditionMessage(e))
  }), warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  cat("\n--", lab, "\n")
  if (is.character(v)) cat(v, "\n") else print(utils::head(v, 4))
  cat("   warnings:", length(w), if (length(w)) substr(w[1], 1, 110), "\n")
}
if (requireNamespace("emmeans", quietly = TRUE)) {
  show("emmeans(~ f) on mu", summary(emmeans::emmeans(fit, ~ f)))
  show("emmeans(~ f, dpar = 'a')",
       summary(emmeans::emmeans(fit, ~ f, dpar = "a")))
}
show("hypothesis(a_f1 = 0)", hypothesis(fit, "a_f1 = 0"))
show("hypothesis(a_f1 + b_Intercept = 0)",
     hypothesis(fit, "a_f1 + b_Intercept = 0"))

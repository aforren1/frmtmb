# test-nl-family.R:93 on hurdle_negbinomial: ubuntu-24.04 (R release)
# warned "Large maximum absolute gradient at the optimum (0.00201 at
# shape_(Intercept))" at 0.79.1, run 42116483453, job 127557370985.
# Fits the file's three models per family on its data (seed 20261029)
# and prints the shape, the warnings, and the file's ratio.
#
#   Rscript dev/ci-hnbprobe.R <lib1;lib2;...>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(args[[1L]], ";", fixed = TRUE)[[1L]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(20261029)
n <- 300
x <- round(rnorm(n), 3)
eta <- 0.3 + 0.5 * x
pos <- exp(0.4 + 0.3 * x + rnorm(n, 0, 0.4))
cnt <- rpois(n, exp(eta))
prob <- plogis(eta)
zero <- runif(n) < 0.2
fams <- c("negbinomial", "zero_inflated_negbinomial", "hurdle_negbinomial")
fit1 <- function(form, f, d) {
  w <- character()
  fit <- withCallingHandlers(
    suppressMessages(frm(form, family = get(f, asNamespace("frmtmb"))(),
                         data = d)),
    warning = function(cnd) {
      w <<- c(w, substr(conditionMessage(cnd), 1, 90))
      invokeRestart("muffleWarning")
    })
  p <- fit$opt$par
  list(ll = as.numeric(logLik(fit)), w = w,
       par = fit$opt$par, conv = fit$opt$convergence,
       msg = fit$opt$message)
}
for (f in fams) {
  d <- data.frame(x = x)
  d$y <- if (f == "hurdle_negbinomial") ifelse(zero, 0L, cnt + 1L) else
    if (f == "zero_inflated_negbinomial") ifelse(zero, 0L, cnt) else cnt
  lin <- fit1(bf(y ~ x), f, d)
  int <- fit1(bf(y ~ 1), f, d)
  nl <- fit1(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), f, d)
  cat(sprintf("%s: lin %.10g int %.10g nl %.10g ratio %.3g\n", f, lin$ll,
              int$ll, nl$ll, abs(nl$ll - lin$ll) / (lin$ll - int$ll)))
  for (nm in c("lin", "int", "nl")) {
    r <- get(nm)
    cat(sprintf("  %-3s conv %d (%s) par %s\n", nm, r$conv, r$msg,
                paste(sprintf("%s=%.4g", names(r$par), r$par),
                      collapse = " ")))
    for (m in r$w) cat("    warning:", m, "\n")
  }
}

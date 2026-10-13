# Reviewer of the 0.69.0 consolidation: the precedence run's fit 58 of
# test-diagnostics-ux.R ("a theta in the flat set ..."): the nonlinear
# bump is numerically zero on every row, so pk, lamp, lsig and the sd of
# pk's random effect do not enter the likelihood. What does the user
# get on this build: which warning and message, and what diagnose()
# lists?
#   Rscript dev/relrev069-flatsd.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), find.package("frmtmb"),
    "\n")
flat_peak_data <- function(n_id = 10, seed = 5) {
  set.seed(seed)
  w <- seq(1, 45, by = 0.25)
  d <- do.call(rbind, lapply(seq_len(n_id), function(i) {
    mu_i <- stats::rnorm(1, 10, 1.1)
    logf <- stats::rnorm(1, 1.2, 0.3) - 1.5 * log(w) +
      1.1 * exp(-0.5 * ((w - mu_i) / 1.6)^2)
    data.frame(id = i, w = w, logw = log(w),
               I = stats::rexp(length(w), rate = 1 / exp(logf)))
  }))
  d$id <- factor(d$id)
  d
}
d <- flat_peak_data()
fit <- withCallingHandlers(frm(
  bf(I ~ apo - chi * logw + exp(lamp) * exp(-0.5 * ((w - pk) / exp(lsig))^2),
     apo ~ 1 + (1 | id), chi ~ 1, pk ~ 1 + (1 | id), lamp ~ 1, lsig ~ 1,
     nl = TRUE),
  family = exponential(link = "log"), data = d),
  warning = function(w) { cat("WARNING:", conditionMessage(w), "\n\n")
    invokeRestart("muffleWarning") },
  message = function(m) { cat("MESSAGE:", conditionMessage(m), "\n")
    invokeRestart("muffleMessage") })
s <- withCallingHandlers(summary(fit),
  warning = function(w) { cat("WARNING (summary):", conditionMessage(w), "\n\n")
    invokeRestart("muffleWarning") },
  message = function(m) { cat("MESSAGE (summary):", conditionMessage(m), "\n")
    invokeRestart("muffleMessage") })
dg <- diagnose(fit, quiet = TRUE)
cat("diagnose()$flat:", dg$flat, "\n")
cat("diagnose()$singular:", unlist(dg$singular), "\n")
cat("max over rows of the bump term:",
    with(d, max(exp(fit$estimates$betad[1]) * 0)), "\n")
print(fit$opt$par)

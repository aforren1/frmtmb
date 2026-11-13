# The third case of backlog 3269: predict() of se(s, sigma = TRUE) with
# sigma estimated near 0 (the construction of dev/ndaterm-probe.R,
# n = 200, seed 7). Usage: Rscript dev/ppdraw-probe4.R <lib or "base">
args <- commandArgs(TRUE)
lib <- if (args[1] != "base") args[1] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
n <- 200
set.seed(7)
dg <- data.frame(x = rnorm(n), s = runif(n, 0.2, 1))
dg$y <- 1 + dg$x + rnorm(n, 0, dg$s)
ndg <- data.frame(x = c(0, 0), s = c(0.05, 5))
fg <- withCallingHandlers(
  frm(bf(y | se(s, sigma = TRUE) ~ x), family = gaussian(), data = dg),
  message = function(m) cat("  M:", conditionMessage(m), "\n"),
  warning = function(w) cat("  W:", conditionMessage(w), "\n"))
print(fg$estimates)
print(sqrt(diag(vcov(fg, full = TRUE))))
set.seed(7)
pg <- withCallingHandlers(
  predict(fg, newdata = ndg, summary = FALSE, ndraws = 400),
  message = function(m) cat("  M:", conditionMessage(m), "\n"),
  warning = function(w) cat("  W:", conditionMessage(w), "\n"))
cat("predict column sd:", apply(pg, 2, sd, na.rm = TRUE), "\n")
sgd <- simulate(fg, nsim = 400, newdata = ndg, seed = 7)
cat("simulate row sd:", apply(as.matrix(sgd), 1, sd), "\n")
print(fitted(fg, newdata = ndg))

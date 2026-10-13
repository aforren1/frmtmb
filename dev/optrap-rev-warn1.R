# Reviewer of lane optrap: the full convergence warning on student,
# sigma ~ x, y * 1e-6, seed 1 (dev/optrap-rev-scale.R's generator).
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1)
d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
b <- rnorm(20)[d$g]
d$y <- (1 + 0.5 * d$x + 0.5 * b + exp(0.3 * d$x) * rt(100, 5)) * 1e-6
f <- withCallingHandlers(frm(bf(y ~ x + (1 | g), sigma ~ x),
                             family = student(), data = d),
  warning = function(w) {cat("W:", conditionMessage(w), "\n");
    invokeRestart("muffleWarning")})
cat("code", f$opt$convergence, "ll", format(logLik(f), digits = 15), "\n")
g <- f$obj$gr(f$opt$par)
print(rbind(par = f$opt$par, grad = g, par_units = f$par_units,
            resp_units = f$resp_units %||% NA))

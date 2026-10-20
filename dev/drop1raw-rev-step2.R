# Reviewer: step() moves through update(evaluate = FALSE) and
# eval.parent(), so does a step read the fit's data or the workspace's?
# Arg: lane. Log: dev/drop1raw-rev-out/step2-lane.txt (seed 21, n = 300)
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-drop1raw-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(21)
n <- 300
d <- data.frame(x = rnorm(n), x2 = rexp(n) + 0.1, x3 = rnorm(n),
                x4 = rnorm(n))
d$z <- 1 + 0.5 * d$x + log(d$x2) + rnorm(n)
f4 <- z ~ x + log(x2) + x3 + x4
fit <- q(frm(bf(f4), data = d))
cat("update(fit, '~ . - x4', evaluate = FALSE) data arg:",
    deparse1(update(fit, "~ . - x4", evaluate = FALSE)$data), "\n")
s_ok <- q(step(fit, trace = 0))
cat("step, data in place:", deparse1(formula(s_ok)), "| logLik",
    sprintf("%.8f", as.numeric(logLik(s_ok))), "\n")
fit_gone <- local({
  d_gone <- d
  q(frm(bf(f4), data = d_gone))
})
r <- tryCatch(q(step(fit_gone, trace = 0)), error = function(e) e)
cat("step, data gone:", if (inherits(r, "error")) paste("ERROR",
    conditionMessage(r)) else deparse1(formula(r)), "\n")
cat("drop1, data gone:", paste(sprintf("%.4f", q(drop1(fit_gone))$AIC),
                                collapse = " "), "\n")
d_keep <- d
d$x2 <- d$x2^2 + 1
r2 <- tryCatch(q(step(fit, trace = 0)), error = function(e) e)
cat("step after d$x2 changed:", if (inherits(r2, "error")) paste("ERROR",
    conditionMessage(r2)) else paste(deparse1(formula(r2)), "| logLik",
    sprintf("%.8f", as.numeric(logLik(r2)))), "\n")
d <- d_keep
# add1() through stats::add1.default()
r3 <- tryCatch(q(add1(fit_gone, ~ . + I(x^2))), error = function(e) e)
cat("add1, data gone:", if (inherits(r3, "error")) paste("ERROR",
    conditionMessage(r3)) else "ok", "\n")

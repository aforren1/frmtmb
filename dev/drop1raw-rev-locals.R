# Reviewer: update() of a fit whose call names its formula and family by
# local symbols, the way a fitting helper or a simulation loop writes
# it. drop1() takes both from the fit; does update()?
# Arg: lane. Log: dev/drop1raw-rev-out/locals-lane.txt (seed 5, n = 200)
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-drop1raw-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(5)
n <- 200
d <- data.frame(x = rnorm(n), x2 = rexp(n) + 0.1, x3 = rnorm(n))
d$z <- 1 + d$x + log(d$x2) + rnorm(n)
fit_with <- function(form, fam, data) frm(form, family = fam, data = data)
fit <- q(fit_with(bf(z ~ x + log(x2) + x3), gaussian(), d))
cat("call:", deparse1(fit$call), "\n")
show <- function(lab, expr) {
  r <- tryCatch(q(expr), error = function(e) e)
  cat(sprintf("%-36s %s\n", lab, if (inherits(r, "error"))
    paste("ERROR", substr(conditionMessage(r), 1, 120)) else
      sprintf("logLik %.8f", as.numeric(logLik(r)))))
}
cat("drop1(fit) AIC:", tryCatch(paste(sprintf("%.4f", q(drop1(fit))$AIC), collapse = " "), error = conditionMessage), "
")
show("update(fit)", update(fit))
show("update(fit, ~ . - x3)", update(fit, ~ . - x3))
show("update(fit, family = student())", update(fit, family = student()))
form <- bf(z ~ x)
show("update(fit), global `form` = z ~ x", update(fit))
fam <- gaussian()
u <- q(update(fit))
cat(sprintf("update(fit), global form = z ~ x and fam: formula %s | logLik %.8f | fit %.8f\n",
            deparse1(formula(u)), as.numeric(logLik(u)),
            as.numeric(logLik(fit))))

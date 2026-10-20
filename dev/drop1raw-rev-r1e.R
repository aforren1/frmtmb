# Reviewer, re-check of fix round 1: the stats idiom of editing the call
# update(evaluate = FALSE) returns, and a bound on a kept scale() term
# in drop1(). Arg: lane. Log: dev/drop1raw-rev-out/r1e-lane.txt (seed 5)
.libPaths(c("C:/Users/adf44/source/r/wt-drop1raw-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(5)
n <- 200
d <- data.frame(x = rnorm(n), x2 = rexp(n) + 0.1, x3 = rnorm(n))
d$z <- 1 + d$x + log(d$x2) + rnorm(n)
d$x3[c(3, 17, 40)] <- NA
fit <- q(frm(bf(z ~ splines::ns(x2, 3) + x + x3), data = d))
cl <- update(fit, evaluate = FALSE)
cl$formula <- bf(z ~ x)
f2 <- q(eval(cl))
cat("edited call: fitted formula", deparse1(formula(f2)), "| recorded call",
    substr(deparse1(getCall(f2)), 1, 120), "\n")
cl2 <- update(fit, evaluate = FALSE)
cl2$formula <- bf(z ~ splines::ns(x2, 3))
cl2$data <- d[!is.na(d$x3), ]
f3 <- q(eval(cl2))
ref <- q(frm(bf(z ~ splines::ns(x2, 3)), data = d[!is.na(d$x3), ]))
cat(sprintf("edited call, other data: logLik %.10f | frm on that data %.10f\n",
            as.numeric(logLik(f3)), as.numeric(logLik(ref))))
cat("  its na_action:", paste(unclass(f3$frame$na_action), collapse = ","),
    "| rows", nrow(f3$frame$data_frame), "| length(residuals)",
    NROW(residuals(f3)), "\n")
fb <- q(frm(bf(z ~ scale(x) + log(x2) + x3), data = d,
            prior = set_prior("", class = "b", coef = "scalex", ub = 0.5)))
msgs <- character(0)
tab <- withCallingHandlers(drop1(fb), message = function(m) {
  msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
}, warning = function(w) invokeRestart("muffleWarning"))
cat("bound ub = 0.5 on b_scalex: fit coef", sprintf("%.6f",
    fixef(fb)["scalex", 1]), "| drop1 messages:", length(msgs), "\n")
cat(paste(" ", msgs), sep = "\n")

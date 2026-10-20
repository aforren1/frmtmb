# Two more readers of the stored model frame: a re_formula term that
# re-evaluates a transformed term, and make_conditions() on a fit.
# Arg: base | lane. Log: dev/drop1raw-out/probe-more-<arm>.txt (seed 5)
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-drop1raw-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(5)
n <- 200
d <- data.frame(x = rnorm(n), x2 = rexp(n) + 0.1,
                g = factor(sample(letters[1:8], n, TRUE)))
d$z <- 1 + d$x + log(d$x2) + rnorm(8)[d$g] + rnorm(8)[d$g] * log(d$x2) +
  rnorm(n)
r <- function(lab, e) cat(sprintf("%-46s %s\n", lab, tryCatch(
  paste("ok", paste(q(e), collapse = " | ")),
  error = function(err) paste("ERROR", substr(conditionMessage(err), 1, 90)))))
fit <- q(frm(bf(z ~ x + log(x2) + (1 + log(x2) | g)), data = d))
r("predict(re_formula = ~ (1 + log(x2) | g))",
  format(length(predict(fit, re_formula = ~ (1 + log(x2) | g))[, 1])))
r("predict(re_formula = ~ (1 | g))",
  format(length(predict(fit, re_formula = ~ (1 | g))[, 1])))
r("make_conditions(fit, 'x2')", format(nrow(make_conditions(fit, "x2"))))
r("make_conditions(fit, 'x')", format(nrow(make_conditions(fit, "x"))))

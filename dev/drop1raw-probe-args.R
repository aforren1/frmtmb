# drop1() and step() on models WITHOUT a transformed term: a family
# named by a symbol of a function's frame, a prior on the dropped
# coefficient, and step() that removes a term. Arg: base | lane.
# Log: dev/drop1raw-out/probe-args-<arm>.txt (seed 5, n = 120)
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-drop1raw-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(5)
n <- 120
d <- data.frame(x = rnorm(n), x3 = rnorm(n))
d$z <- 1 + d$x + rnorm(n)
d$cnt <- rpois(n, exp(0.5 + 0.3 * d$x))
r <- function(lab, e) cat(sprintf("%-40s %s\n", lab, tryCatch(
  paste("ok", paste(q(e), collapse = " | ")),
  error = function(err) paste("ERROR", substr(conditionMessage(err), 1, 90)))))
f_local <- function(form, dat) {
  fam <- poisson()
  frm(bf(form), family = fam, data = dat)
}
fit <- q(f_local(cnt ~ x + x3, d))
r("drop1, family = <local symbol>", rownames(drop1(fit)))
fp <- q(frm(bf(z ~ x + x3), data = d,
            prior = set_prior("normal(0, 1)", class = "b", coef = "x3")))
r("drop1, prior on the dropped coef", rownames(drop1(fp)))
fs <- q(frm(bf(z ~ x + x3), data = d))
r("step() removes x3", deparse1(formula(q(step(fs, trace = 0)))))
r("lm step() for reference", deparse1(formula(step(lm(z ~ x + x3, d),
                                                  trace = 0))))

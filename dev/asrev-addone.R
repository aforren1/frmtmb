## REVIEW 2.5e: a base-R code path that passes an argument the new
## refusal rejects. stats::add1.default and stats::step call
## nobs(object, use.fallback = TRUE); nobs.frmtmb_fit now refuses it.
## Usage: Rscript dev/asrev-addone.R <lib>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1],
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("frmtmb from:", dirname(system.file(package = "frmtmb")), "\n")
set.seed(11)
n <- 90
dd <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:9, each = 10)))
dd$y <- 0.3 + 0.5 * dd$x + rnorm(9, 0, .6)[dd$g] + rnorm(n, 0, .5)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

p <- function(lab, expr) {
  r <- tryCatch(force(expr), error = function(e) e)
  cat(sprintf("%-40s %s\n", lab,
              if (inherits(r, "error"))
                paste("ERROR:", gsub("\n", " ", conditionMessage(r)))
              else paste("ok", paste(class(r), collapse = "/"))))
}
p("nobs(fit, use.fallback = TRUE)", stats::nobs(fit, use.fallback = TRUE))
p("add1(fit, ~ . + z)", stats::add1(fit, ~ . + z))
p("drop1(fit)", stats::drop1(fit))
p("step(fit)", utils::capture.output(stats::step(fit)))
p("extractAIC(fit)", stats::extractAIC(fit))
p("print(fit, digits = 3)", utils::capture.output(print(fit, digits = 3)))
p("print(summary(fit), digits = 3)",
  utils::capture.output(print(summary(fit), digits = 3)))
p("summary(fit, prob = 0.9)", summary(fit, prob = 0.9))
p("VarCorr(fit, sigma = 1)", VarCorr(fit, sigma = 1))
p("fixef(fit, pars = 'x')", fixef(fit, pars = "x"))
p("nobs(fit, resp = 'y')", stats::nobs(fit, resp = "y"))
p("family(fit, resp = 'y')", stats::family(fit, resp = "y"))
p("vcov(fit, correlation = TRUE)", stats::vcov(fit, correlation = TRUE))
cat("DONE\n")

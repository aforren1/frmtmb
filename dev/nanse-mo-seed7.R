# What a user sees at seed 7 of brms_monotonic's data code: the fit,
# summary(), vcov(), with every condition printed.
#   Rscript dev/nanse-mo-seed7.R [lib] [seed]
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r5"
s <- if (length(args) > 1) as.integer(args[2]) else 7L
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r5",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "seed", s, "\n")
set.seed(s)
lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
d <- data.frame(income, ls)
d$age <- rnorm(100, mean = 40, sd = 10)
show <- function(tag, expr) {
  cat("\n-----", tag, "\n")
  withCallingHandlers(expr, warning = function(w) {
    cat("WARNING:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }, message = function(m) {
    cat("MESSAGE:", conditionMessage(m))
    invokeRestart("muffleMessage")
  })
}
f <- show("frm", frm(ls ~ mo(income) * age, data = d))
show("summary", print(summary(f)))
show("vcov", print(vcov(f)))
show("confint", print(confint(f)))

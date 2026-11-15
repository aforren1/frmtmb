# Lane famend: test-prior-table-brms.R's "zero-density ends" case on one
# build: beta(2, 2) on a slope the data put below 0, with no bound, with
# lb = 1e-8, with the closed box [0, 1], and with a caller's start on 0.
#   Rscript dev/famend-probe-beta22.R <lib or "base">
a <- commandArgs(trailingOnly = TRUE)
base <- "C:/Users/adf44/source/r/rellib-r21"
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (a[1] != "base") a[1], base, ulib))
suppressPackageStartupMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
set.seed(1)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 12)))
d$y <- 1 + 0.3 * d$z + rnorm(10)[d$g] * 0.3 + rnorm(n, sd = 0.3)
d$yneg <- d$y - 0.8 * d$z
sp <- function(...) set_prior("beta(2, 2)", class = "b", coef = "z", ...)
run <- function(lab, expr) {
  w <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }), error = function(e) structure(conditionMessage(e), class = "err"))
  if (inherits(r, "err")) {
    cat(sprintf("%-26s ERROR %s\n", lab, substr(r, 1, 100)))
  } else {
    cat(sprintf("%-26s code %d objective %.10g z %.6g | %s\n", lab,
                r$opt$convergence, r$opt$objective,
                fixef(r)["z", "Estimate"],
                paste(substr(w, 1, 40), collapse = " || ")))
  }
}
run("lb = 1e-8, ub = 1", frm(yneg ~ z, data = d, prior = sp(lb = 1e-8, ub = 1)))
run("no bound", frm(yneg ~ z, data = d, prior = sp()))
run("lb = 0, ub = 1", frm(yneg ~ z, data = d, prior = sp(lb = 0, ub = 1)))
run("lb = 0, ub = 1, start 0", frm(yneg ~ z, data = d,
                                   prior = sp(lb = 0, ub = 1),
                                   start = list(beta = c(1, 0))))
run("no bound, start 0", frm(yneg ~ z, data = d, prior = sp(),
                             start = list(beta = c(1, 0))))

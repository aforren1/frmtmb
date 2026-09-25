# zero_one_inflated_beta()'s fit-end check on four response shapes,
# with coi estimated and with coi held at 0.5.
#
#   Rscript dev/fams-coi-check.R > dev/fams-coi-check.txt
lib <- Sys.getenv("FRMTMB_LIB", "/opt/rlib/lane-fams")
.libPaths(c(if (lib != "base") lib, "/opt/rlib/base", "/opt/rlib/deps",
            "/opt/r/lib/R/library"))
suppressMessages(library(frmtmb))
set.seed(4)
n <- 300
x <- rnorm(n)
yb <- rbeta(n, 2, 3)
edge <- runif(n) < 0.2
ds <- list(zeros_only = ifelse(edge, 0, yb),
           ones_only = ifelse(edge, 1, yb),
           none = yb,
           both = ifelse(edge, rbinom(n, 1, 0.5), yb))
fit_w <- function(form, d) {
  w <- character(0)
  fit <- withCallingHandlers(
    frm(form + zero_one_inflated_beta(), data = d),
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond))
      invokeRestart("muffleWarning")
    })
  list(fit = fit, w = w)
}
for (nm in names(ds)) {
  d <- data.frame(x = x, y = ds[[nm]])
  a <- fit_w(bf(y ~ x), d)
  b <- fit_w(bf(y ~ x, coi = 0.5), d)
  se <- sqrt(diag(vcov(a$fit)))
  cat(sprintf("== %s: %d zeros, %d ones\n", nm, sum(d$y == 0),
              sum(d$y == 1)))
  cat("  coi estimated: SE(x) ", format(se[["x"]], digits = 4),
      "; warnings: ", if (length(a$w)) a$w else "none", "\n", sep = "")
  cat("  coi = 0.5:     SE(x) ",
      format(sqrt(diag(vcov(b$fit)))[["x"]], digits = 4),
      "; warnings: ", if (length(b$w)) b$w else "none", "\n", sep = "")
}

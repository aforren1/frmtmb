# Review of lane gpgap: does gp_nugget_message() break a fit the base build
# completes? sigma with no estimated coefficient. Seed 5501 (the warn script's
# design).
#
#   Rscript dev/review-gpgap-offset.R <lib stack, ; separated>
lib <- strsplit(commandArgs(trailingOnly = TRUE)[1], ";")[[1]]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(5501)
n <- 60
x <- round(runif(n, 0, 6), 3)
e <- rnorm(n)
d <- data.frame(x = x, y = sin(2 * x) + 0.3 * e, wsm = 0.3, wlow = 1e-3)
d$ylow <- sin(2 * x) + 1e-3 * e
try_one <- function(tag, expr) {
  r <- tryCatch({
    f <- suppressWarnings(expr)
    sprintf("ok, ll %.3f sdgp %.4f conv %d", -f$opt$objective,
            exp(f$estimates$theta[[1]]), f$opt$convergence)
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-46s %s\n", tag, r))
}
try_one("sigma ~ 0 + offset(log(wsm)), ordinary noise",
        frm(bf(y ~ gp(x), sigma ~ 0 + offset(log(wsm))), data = d))
try_one("sigma ~ 0 + offset(log(wlow)), low noise",
        frm(bf(ylow ~ gp(x), sigma ~ 0 + offset(log(wlow))), data = d))
try_one("sigma ~ 1 + offset(log(wsm))",
        frm(bf(y ~ gp(x), sigma ~ 1 + offset(log(wsm))), data = d))
try_one("hsgp, sigma ~ 0 + offset(log(wsm))",
        frm(bf(y ~ gp(x, k = 20), sigma ~ 0 + offset(log(wsm))), data = d))
f <- tryCatch(frm(bf(ylow ~ gp(x), sigma ~ 1 + offset(log(wlow))), data = d),
              warning = function(w) paste("warned:", class(w)[1]))
cat("sigma ~ 1 + offset(log(wlow)) low noise, first warning:",
    if (is.character(f)) f else "none", "\n")

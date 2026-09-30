# The behavioral failures on the base build (rellib-r4, frmtmb 0.66.0)
# that the lane's tests pin, beside the lane build's answer to the same
# call. Seed 20260930. Usage: Rscript dev/ordinal-base-behavior.R
# lane|base. Output: dev/ordinal-log-base-behavior-<arm>.txt
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r4",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) {
  libs <- c("C:/Users/adf44/source/r/wt-ordinal-lib", libs)
}
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
set.seed(20260930)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n))
u <- stats::rlogis(n) / exp(0.4 * d$z) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
d$yh <- ifelse(runif(n) < 0.25, 0L, d$y)
show <- function(lab, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  if (inherits(r, "frmtmb_fit")) {
    r <- paste("fits; logLik", format(as.numeric(logLik(r)), digits = 10))
  }
  cat(sprintf("%-46s %s\n", lab, substr(paste(r, collapse = " "), 1, 150)))
}
show("cumulative(), disc ~ 0 + z",
     frm(bf(y ~ x, disc ~ 0 + z), family = cumulative(), data = d))
show("hurdle_cumulative(threshold = 'equidistant')",
     frm(yh ~ x, family = hurdle_cumulative(threshold = "equidistant"),
         data = d))
show("hurdle_cumulative(threshold = 'sum_to_zero')",
     frm(yh ~ x, family = hurdle_cumulative(threshold = "sum_to_zero"),
         data = d))
show("acat('probit')", frm(y ~ x, family = acat("probit"), data = d))
show("brmsfamily('sratio', threshold = 'equidistant')",
     frm(y ~ x, family = brmsfamily("sratio", threshold = "equidistant"),
         data = d))
show("flexible cumulative (unchanged)",
     frm(y ~ x, family = cumulative(), data = d))
show("flexible acat, cs(z) (unchanged)",
     frm(y ~ x + cs(z), family = acat(), data = d))

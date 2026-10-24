# Called by dev/splitsurv-floor.sh with the scratch library: load
# frmtmb.survival over rellib-r11's core and try what it needs from core.
fl <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(fl, "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
cat("frmtmb", format(packageVersion("frmtmb")), dirname(find.package("frmtmb")),
    "\n")
step <- function(label, expr) {
  r <- tryCatch({ force(expr); "ok" },
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-28s %s\n", label, substr(r, 1, 160)))
}
step("library(frmtmb.survival)", suppressMessages(library(frmtmb.survival)))
set.seed(1)
d <- data.frame(x = rnorm(120))
d$time <- rexp(120, exp(-0.5 + 0.7 * d$x))
step("cox()", f <- cox())
step("frm(family = cox())", fit <- frm(time ~ x, family = cox(), data = d))
step("family = \"cox\"", frm(time ~ x, family = "cox", data = d))
step("frm_curve exists in core", stopifnot(exists("frm_curve")))

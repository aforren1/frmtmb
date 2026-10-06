# Lane fixes' consolidation item 1 on the release build after lane
# nanse's merge: the exact ridge a + b with c0 near 0 (log(c0) undefined
# a step below it), from dev/fixes-rev3-flat.R case 5, and its control
# with c0 near 1. Prints every warning and the fixef standard errors.
#   Rscript dev/rel068-item1.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
run <- function(lab, expr) {
  w <- character()
  fit <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  se <- withCallingHandlers(fixef(fit)[, "Est.Error"], warning = function(x) {
    w <<- c(w, paste("[fixef]", conditionMessage(x))); invokeRestart("muffleWarning") })
  cat(lab, "| SE", format(se, digits = 3), "\n")
  for (x in unique(w)) cat("   W:", substr(x, 1, 220), "\n")
  if (!length(w)) cat("   (no warning)\n")
}
set.seed(45)
d5 <- data.frame(x = rnorm(200))
d5$y <- 1 + 0.5 * d5$x + log(5e-5) + rnorm(200, 0, 0.3)
run("c0 near 0", frm(bf(y ~ a + b + log(c0), a ~ 1 + x, b ~ 1 + x, c0 ~ 1,
                        nl = TRUE), data = d5,
                     start = list(beta = c(0.5, 0.25, 0.5, 0.25, 5e-5))))
run("c0 near 1 (control)",
    frm(bf(y ~ a + b + log(c0), a ~ 1 + x, b ~ 1 + x, c0 ~ 1, nl = TRUE),
        data = transform(d5, y = y - log(5e-5)),
        start = list(beta = c(0.5, 0.25, 0.5, 0.25, 1))))

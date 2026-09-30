# Lane ceplot: what a recorded base-graphics display list holds, to count
# the lines and rug marks plot() draws in a test.
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(5)
dd <- data.frame(x = rnorm(60))
dd$y <- rnorm(60, 1 + 0.5 * dd$x)
fit <- frm(bf(y ~ x), family = gaussian(), data = dd)
ce <- conditional_effects(fit, "x", resolution = 5)
ops <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  grDevices::dev.control("enable")
  force(expr)
  rp <- grDevices::recordPlot()
  vapply(rp[[1]], function(e) {
    f <- e[[2]][[1]]
    if (is.list(f) && !is.null(f$name)) f$name else paste(deparse(f), collapse = "")
  }, "")
}
print(table(ops(plot(ce, ask = FALSE))))
print(table(ops(plot(ce, ask = FALSE, rug = TRUE))))

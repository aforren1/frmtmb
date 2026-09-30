# Reviewer re-check (lane ceplot, punch round 1): ask = NULL, and the
# ignored-parameter warning on each plot entry point, with its guard
# absent (a plain call warns nothing). Data seed 5.
#   Rscript dev/ceplot-rev-misc3.R > dev/ceplot-rev-log/p1/misc3.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(5)
dd <- data.frame(x = rnorm(80), f = factor(rep(c("a", "b"), 40)))
dd$y <- rnorm(80, 1 + 0.5 * dd$x + (dd$f == "b"))
fit <- frm(bf(y ~ x * f), family = gaussian(), data = dd)
ce <- conditional_effects(fit, effects = c("x", "f"), resolution = 5)
h <- hypothesis(fit, c("x > 0", "fb = 0"))
p <- plot(ce, plot = FALSE)
ph <- plot(h, plot = FALSE)
grDevices::pdf(NULL)
run <- function(label, expr) {
  w <- character(0)
  r <- tryCatch(withCallingHandlers({
    force(expr)
    "ok"
  }, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }), error = function(e) paste("ERROR", conditionMessage(e)))
  cat(sprintf("%-34s %s | warnings %d%s\n", label, substr(r, 1, 80),
              length(w), if (length(w)) paste0(": ", substr(w[1], 1, 70))
              else ""))
}
run("plot(ce, ask = NULL)", plot(ce, ask = NULL))
run("plot(h, ask = NULL)", plot(h, ask = NULL))
run("plot(ce, ask = FALSE)", plot(ce, ask = FALSE))
run("plot(h, ask = FALSE)", plot(h, ask = FALSE))
run("print(ce, ask = FALSE)", print(ce, ask = FALSE))
run("print(p$x)", print(p$x))
run("plot(p$x)", plot(p$x))
run("print(ph[[1]])", print(ph[[1]]))
run("plot(h, col = 'red')", plot(h, ask = FALSE, col = "red"))
run("print(p$x, main = 't')", print(p$x, main = "t"))
run("plot(p$x, col = 'red')", plot(p$x, col = "red"))
run("plot(ph[[1]], lwd = 2)", plot(ph[[1]], lwd = 2))
run("plot(p$x, foo = 1)", plot(p$x, foo = 1))
run("plot(ce, ask = 'yes')", plot(ce, ask = "yes"))
run("plot(ce, do_plot = FALSE)", plot(ce, do_plot = FALSE))
invisible(dev.off())

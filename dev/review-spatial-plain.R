.libPaths(c("C:/Users/adf44/source/r/wt-spatial-lib", "C:/Users/adf44/source/r/rellib-r11", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1); d <- data.frame(g = gl(10, 5), t = rep(1:5, 10), x = rnorm(50), g1 = gl(10, 5), g2 = factor(sample(1:10, 50, TRUE)))
d$y <- rnorm(50) + rnorm(10)[d$g]
for (fo in c("y ~ x + ar()", "y ~ x + ar(gr = g)", "y ~ x + ma()", "y ~ x + arma()", "y ~ x + cosy()", "y ~ x + cosy(gr = g)", "y ~ x + (1 | gr(g))", "y ~ x + (1 | mm(g1, g2))")) {
  b <- tryCatch({brms::standata(as.formula(fo), data = d); "ok"}, error = function(e) paste("ERR", substr(conditionMessage(e), 1, 90)))
  f <- tryCatch({suppressWarnings(frm(as.formula(fo), data = d)); "ok"}, error = function(e) paste("ERR", substr(conditionMessage(e), 1, 120)))
  cat(sprintf("%-28s brms: %-10s frmtmb: %s\n", fo, b, f))
}

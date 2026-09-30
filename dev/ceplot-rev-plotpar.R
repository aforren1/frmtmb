# Reviewer check (lane ceplot): par() layout settings after each plot()
# call (usr and friends change with any plot, so only the settings a
# caller would notice are compared: mfrow, mfcol, mar, oma, cex, ask).
#   Rscript dev/ceplot-rev-plotpar.R > dev/ceplot-rev-log/plotpar.txt
# Data seed 5.
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(5)
dd <- data.frame(x = rnorm(80), z = rnorm(80),
                 f = factor(rep(c("a", "b"), 40)))
dd$y <- rnorm(80, 1 + 0.5 * dd$x + (dd$f == "b") + 0.3 * dd$x * dd$z, 1)
fit <- frm(bf(y ~ x * f + x * z), family = gaussian(), data = dd)
ce <- conditional_effects(fit, effects = c("x", "f", "x:f"), resolution = 5)
cc <- conditional_effects(fit, "x", resolution = 5,
                          conditions = data.frame(f = c("a", "b")))
su <- conditional_effects(fit, "x:z", surface = TRUE, resolution = 5,
                          conditions = data.frame(f = c("a", "b")))
h <- hypothesis(fit, c("x > 0", "fb = 0", "x + z = 0"))
keep <- c("mfrow", "mfcol", "mar", "oma", "cex", "ask")
grDevices::pdf(NULL)
p0 <- par(keep)
calls <- list(
  ce = quote(plot(ce, ask = FALSE)),
  facets = quote(plot(cc, ask = FALSE)),
  facets_base = quote(plot(cc, ask = FALSE, rug = TRUE)),
  surface = quote(plot(su, ask = FALSE)),
  hyp = quote(plot(h, ask = FALSE, nvariables = 2)),
  print_obj = quote(print(plot(cc, plot = FALSE, rug = TRUE)[[1]])))
for (nm in names(calls)) {
  invisible(eval(calls[[nm]]))
  p1 <- par(keep)
  ch <- names(Filter(isFALSE, Map(identical, p0, p1)))
  cat(sprintf("%-12s par changed: %s\n", nm,
              if (length(ch)) paste(ch, collapse = ",") else "none"))
  par(p0)
}
dev.off()

# Lane ceplot: render the new plot() arguments to PNG files for a visual
# check. Data seeds 5 and 1.
#   Rscript dev/ceplot-plots.R
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
out <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot/dev/ceplot-log"
set.seed(5)
dd <- data.frame(x = rnorm(120), z = rnorm(120),
                 f = factor(rep(c("a", "b"), 60)))
dd$y <- rnorm(120, 1 + 0.5 * dd$x + (dd$f == "b") + 0.3 * dd$x * dd$z, 1)
fit <- frm(bf(y ~ x * f + x * z), family = gaussian(), data = dd)
ce <- conditional_effects(fit, effects = c("x", "f", "x:f"))
png(file.path(out, "plot-1-default.png"), 900, 300)
par(mfrow = c(1, 3))
p <- plot(ce, ask = FALSE, plot = FALSE)
for (q in p) print(q)
dev.off()
png(file.path(out, "plot-2-args.png"), 900, 300)
par(mfrow = c(1, 3))
p <- plot(ce, points = TRUE, rug = TRUE, ask = FALSE, plot = FALSE,
          line_args = list(colour = "red", linewidth = 3, fill = "orange",
                           alpha = 0.3),
          point_args = list(size = 1, colour = "blue", width = 0.2),
          errorbar_args = list(width = 0.4, colour = "darkgreen"),
          cat_args = list(size = 2, shape = 17),
          rug_args = list(sides = "bl", colour = "purple"))
for (q in p) print(q)
dev.off()
cs <- conditional_effects(fit, "x:z", surface = TRUE, resolution = 20)
png(file.path(out, "plot-3-surface.png"), 600, 300)
par(mfrow = c(1, 2))
plot(cs, ask = FALSE)
plot(cs, stype = "raster", ask = FALSE)
dev.off()
cb <- conditional_effects(fit, "x", band = "boot", boot = 30, seed = 1,
                          spaghetti = TRUE)
png(file.path(out, "plot-4-spaghetti.png"), 600, 300)
par(mfrow = c(1, 2))
plot(cb, ask = FALSE, spaghetti_args = list(colour = "red"))
plot(cb, ask = FALSE, mean = FALSE)
dev.off()
cc <- conditional_effects(fit, "x", conditions = data.frame(
  f = c("a", "b", "a", "b"), z = c(-1, -1, 1, 1)))
png(file.path(out, "plot-5-facets.png"), 600, 600)
plot(cc, ask = FALSE, facet_args = list(nrow = 1, scales = "free"),
     points = TRUE)
dev.off()
h <- hypothesis(fit, c("x > 0", "fb = 0",
                       "x + fb + z + x:fb + x:z = 0.123456789"))
png(file.path(out, "plot-6-hyp.png"), 500, 600)
plot(h, ask = FALSE, chars = 20)
dev.off()
cat("done\n")

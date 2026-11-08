# predictions() standard error against fitted()'s Est.Error on the cases the
# lane leaves to a warning (a latent residual, me(), quadrature) and on car()
# from frmtmb.spatial, one arm.
#   MFX_ARM=before|after Rscript dev/mfxse-ext.R
# Data seeds 9101 (latent), 31 (car, the spatial suite's), 7411 (me, quad).
source("C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/mfxse-lib.R")
suppressPackageStartupMessages(library(frmtmb.spatial))
cat("spatial:", find.package("frmtmb.spatial"), "\n")
show <- function(lab, fit, nd, ...) {
  ws <- character(0)
  p <- withCallingHandlers(predictions(fit, newdata = nd, ...),
                           warning = function(w) {
                             ws <<- c(ws, substr(conditionMessage(w), 1, 90))
                             invokeRestart("muffleWarning")
                           })
  ft <- quiet(fitted(fit, newdata = nd, ...))
  fse <- if (length(dim(ft)) == 3L) as.numeric(ft[, "Est.Error", ]) else
    ft[, "Est.Error"]
  cat(sprintf("%-10s ratio %s\n           warn: %s\n", lab,
              paste(sprintf("%.4f", p$std.error / fse), collapse = " "),
              paste(unique(ws), collapse = " | ")))
}
# a latent residual beside (1 | s)
set.seed(9101)
d <- expand.grid(t = 1:25, s = factor(1:8))
d$x <- rnorm(nrow(d))
re <- rnorm(8, 0, 0.8)
e <- ave(rnorm(nrow(d)), d$s, FUN = function(z) stats::filter(z, 0.5, "recursive"))
d$y <- rpois(nrow(d), exp(1 + 0.3 * d$x + re[d$s] + 0.3 * e))
f <- quiet(frm(y ~ x + (1 | s) + ar(time = t, gr = s, cov = TRUE), data = d,
               family = poisson()))
show("latres", f, data.frame(t = 26, s = factor(c(1, 2)), x = c(0, 1)))
# car(), esicar and icar, on the spatial suite's lattice
lat <- expand.grid(r = 1:3, c = 1:3)
W <- matrix(0L, nrow(lat), nrow(lat))
for (i in seq_len(nrow(lat))) for (j in seq_len(nrow(lat))) {
  if (abs(lat$r[i] - lat$r[j]) + abs(lat$c[i] - lat$c[j]) == 1L) W[i, j] <- 1L
}
dimnames(W) <- list(as.character(seq_len(nrow(lat))),
                    as.character(seq_len(nrow(lat))))
set.seed(31)
m <- 20
dc <- data.frame(g = factor(rep(seq_len(nrow(lat)), each = m)),
                 x = rnorm(nrow(lat) * m))
u <- rnorm(nrow(lat), 0, 0.8)
dc$y <- 0.6 * dc$x + (u - mean(u))[as.integer(dc$g)] + rnorm(nrow(dc))
ndc <- dc[c(1, 25, 45, 170), c("x", "g")]
for (ty in c("esicar", "icar", "escar")) {
  fc <- quiet(frm(bf(y ~ x + car(W, gr = g, type = ty)), data = dc))
  show(ty, fc, ndc)
}
# me(), and quadrature
set.seed(7411)
n <- 300
dm <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(sample(1:12, n, TRUE)))
dm$sx <- 0.5
dm$xo <- dm$x + rnorm(n, 0, 0.5)
dm$cnt <- rpois(n, exp(0.3 + 0.4 * dm$x + 0.3 * dm$z + rnorm(12, 0, 0.6)[dm$g]))
fm <- quiet(frm(bf(cnt ~ me(xo, sx) + z + (1 | g)), family = poisson(),
                data = dm))
show("me+g", fm, data.frame(xo = c(0, 1), sx = 0.5, z = 0, g = factor(1:2)))
fq <- quiet(frm(cnt ~ x + z + (1 | g), family = poisson(), data = dm,
                quadrature = TRUE))
show("quad", fq, dm[1:3, ])
fl <- quiet(frm(cnt ~ x + z + (1 | g), family = poisson(), data = dm))
show("laplace", fl, dm[1:3, ])
show("newlevel", fl, data.frame(x = 0, z = 0, g = factor(c("1", "zz"))),
     allow_new_levels = TRUE)

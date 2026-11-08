# What marginaleffects' predictions() standard error holds against
# fitted(newdata = )'s Est.Error, per kind of term, on one build.
#   MFX_ARM=before|after Rscript dev/mfxse-cmp.R
# Data seed 7411 (dev/predread-mfx-se.R's design) and 8150 for gp().
source("C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/mfxse-lib.R")
set.seed(7411)
n <- 400
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(sample(1:12, n, TRUE)))
ug <- rnorm(12, 0, 0.8)
d$cnt <- rpois(n, exp(0.5 + 0.4 * d$x + ug[d$g]))
d$yg <- 1 + 0.5 * d$x + ug[d$g] + rnorm(n, 0, 0.6)
d$lg <- exp(0.5 + 0.4 * d$x + ug[d$g] + rnorm(n, 0, 0.6))
d$sm <- sin(2 * d$x) + rnorm(n, 0, 0.4)
d$bin <- rbinom(n, 1, plogis(0.3 * d$x + ug[d$g]))
nd <- data.frame(x = c(-1, 0, 1, 2.5), z = 0, g = factor(c(1, 2, 3, 3), 1:12))
fits <- list(
  poisson_g = quiet(frm(cnt ~ x + (1 | g), family = poisson(), data = d)),
  poisson_gx = quiet(frm(cnt ~ x + (1 + x | g), family = poisson(),
                         data = d)),
  gaussian_g = quiet(frm(yg ~ x + (1 | g), data = d)),
  gaussian_reml = quiet(frm(yg ~ x + (1 | g), data = d, REML = TRUE)),
  lognormal_g = quiet(frm(lg ~ x + (1 | g), family = lognormal(), data = d)),
  bernoulli_g = quiet(frm(bin ~ x + (1 | g), family = bernoulli(), data = d)),
  s = quiet(frm(sm ~ s(x), data = d)),
  gp15 = quiet(frm(sm ~ gp(x, k = 15), data = d)),
  poisson_plain = quiet(frm(cnt ~ x, family = poisson(), data = d))
)
set.seed(8150)
m <- 90
dg <- data.frame(x = runif(m, 0, 6), z = runif(m, 0, 6))
dg$y <- sin(dg$x) + 0.3 * dg$z + rnorm(m, sd = 0.3)
fits$gp_exact <- quiet(frm(y ~ z + gp(x), data = dg))
ndg <- data.frame(x = c(dg$x[5], 3.33, 7), z = 3)
show <- function(lab, fit, nd) {
  ws <- character(0)
  p <- withCallingHandlers(predictions(fit, newdata = nd),
                           warning = function(w) {
                             ws <<- c(ws, class(w)[1])
                             invokeRestart("muffleWarning")
                           })
  ft <- quiet(fitted(fit, newdata = nd))
  r <- p$std.error / ft[, "Est.Error"]
  cat(sprintf("%-14s mfx %s\n%-14s fit %s\n%-14s ratio %s  warn: %s\n",
              lab, paste(sprintf("%.4f", p$std.error), collapse = " "), "",
              paste(sprintf("%.4f", ft[, "Est.Error"]), collapse = " "), "",
              paste(sprintf("%.4f", r), collapse = " "),
              paste(unique(ws), collapse = ",")))
}
for (k in names(fits)) {
  show(k, fits[[k]], if (k == "gp_exact") ndg else nd)
}
cat("\nre_formula = NA:\n")
for (k in c("poisson_gx", "gaussian_g", "s")) {
  p <- quiet(predictions(fits[[k]], newdata = nd, re_formula = NA))
  ft <- quiet(fitted(fits[[k]], newdata = nd, re_formula = NA))
  cat(sprintf("%-14s ratio %s\n", k, paste(sprintf("%.4f", p$std.error /
                                                      ft[, "Est.Error"]),
                                           collapse = " ")))
}
cat("\nhypotheses(fit) on poisson_gx:\n")
print(as.data.frame(quiet(hypotheses(fits$poisson_gx)))[, c("term",
                                                            "estimate",
                                                            "std.error")])
print(sqrt(diag(vcov(fits$poisson_gx))))

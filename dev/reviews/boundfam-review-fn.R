# Review of lane boundfam: false negatives, the message/warning split,
# the lane's before/after cases, and what predict() does on a limit fit.
#
#   Rscript dev/reviews/boundfam-review-fn.R <lane|base> [case ...]
arm <- commandArgs(trailingOnly = TRUE)[1]
cases <- commandArgs(trailingOnly = TRUE)[-1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-boundfam-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-boundfam/dev/reviews/boundfam-review-fp.R")
# reuse conds(), show(), ref() from the FP script
eval(parse(text = src[grep("^conds <- ", src):(grep("^ok <- ", src) - 1L)]))
ok <- function(c) !length(cases) || c %in% cases
F <- function(...) conds(frm(...))

# the lane's own before/after cases
if (ok("pois0")) {
  set.seed(2)
  d <- data.frame(x = rnorm(50), y = 0L)
  show("pois0", F(y ~ x, family = poisson(), data = d))
  ref("glm", glm(y ~ x, poisson, d))
}
if (ok("cens")) for (s in c(1L, 3L, 6L)) {
  set.seed(s)
  d <- data.frame(x = rnorm(60))
  d$y <- 1 + d$x + rnorm(60)
  d$cen <- "right"
  show(paste0("cens s", s), F(y | cens(cen) ~ x, data = d))
}
if (ok("maling")) {
  d <- data.frame(k = c(45, 45, 44, 45, 44, 45, 45, 45, 45, 45, 30, 20, 6,
                        44, 44, 27, 25, 17, 14, 27, 35, 30), n = 45L)
  show("maling", F(bf(k | trials(n) ~ 1),
                   family = mixture(beta_binomial, beta_binomial), data = d,
                   start = list(beta = qlogis(c(0.55, 0.95)),
                                betad = c(log(c(15, 15)), 0))))
}
if (ok("nbpois")) for (s in 2:3) {
  set.seed(s)
  d <- data.frame(x = rnorm(200))
  d$y <- rpois(200, exp(1 + 0.4 * d$x))
  show(paste0("nbpois s", s), F(y ~ x, family = negbinomial(), data = d))
}
# nb on poisson data, seed 17 of the FP script: 'nocurv' not 'limit'
if (ok("nb17")) {
  set.seed(17)
  d <- data.frame(x = rnorm(200))
  d$y <- rpois(200, exp(1 + 0.4 * d$x))
  r <- F(y ~ x, family = negbinomial(), data = d)
  show("nb17", r)
  fit <- r$value
  p <- fit$opt$par
  for (t in c(-12, -2, 2, 14)) {
    q <- p; q[3] <- p[3] + t
    cat(sprintf("  nb17 obj change at log shape %+g: %.6g\n", t, fit$obj$fn(q) - fit$obj$fn(p)))
  }
}
# zero-inflation slope separating: zi -> 1 below x = 0, zi -> 0 above
if (ok("zisep")) for (s in 1:2) {
  set.seed(s)
  n <- 200
  d <- data.frame(x = rnorm(n))
  d$y <- ifelse(d$x < 0, 0L, rpois(n, exp(2)))
  show(paste0("zisep s", s), F(bf(y ~ 1, zi ~ x), family = zero_inflated_poisson(), data = d))
  # the bernoulli analogue: the same separation in a mean model
  d$z <- as.integer(d$y == 0)
  show(paste0("zisep-bern s", s), F(z ~ x, family = bernoulli(), data = d))
}
# zero inflation at zero (no excess zeros): the family limit
if (ok("zi0")) for (s in 1:2) {
  set.seed(s)
  d <- data.frame(x = rnorm(300))
  d$y <- rpois(300, exp(1 + 0.3 * d$x))
  show(paste0("zi0 s", s), F(y ~ x, family = zero_inflated_poisson(), data = d))
}
# a variance component at zero, for the 0.81.0 comparison of the report
if (ok("sd0")) {
  set.seed(3)
  d <- data.frame(x = rnorm(200), g = factor(rep(1:20, each = 10)))
  d$y <- rpois(200, exp(1 + 0.3 * d$x))
  show("sd0 poisson (1|g)", F(y ~ x + (1 | g), family = poisson(), data = d))
}
# joint runaway in the mean: one group all zero, poisson y ~ g
if (ok("grp0")) {
  set.seed(5)
  d <- data.frame(g = factor(rep(c("a", "b"), each = 50)))
  d$y <- ifelse(d$g == "a", 0L, rpois(100, 3))
  show("grp0 poisson y ~ g", F(y ~ g, family = poisson(), data = d))
  ref("glm", glm(y ~ g, poisson, d))
  show("grp0 poisson y ~ 0 + g", F(y ~ 0 + g, family = poisson(), data = d))
  show("grp0 zip y ~ g", F(y ~ g, family = zero_inflated_poisson(), data = d))
  show("grp0 negbinomial y ~ g", F(y ~ g, family = negbinomial(), data = d))
}
# mixture weight at zero: two gaussians fit to one, several seeds
if (ok("mix1")) for (s in 1:4) {
  set.seed(s)
  d <- data.frame(y = rnorm(300, 1, 1))
  show(paste0("mix1 s", s), F(y ~ 1, family = mixture(gaussian, gaussian), data = d))
}
# a mixture of two poissons on one poisson, with a covariate
if (ok("mixp")) for (s in 1:3) {
  set.seed(s)
  d <- data.frame(x = rnorm(300))
  d$y <- rpois(300, exp(1 + 0.3 * d$x))
  show(paste0("mixp s", s), F(y ~ x, family = mixture(poisson, poisson), data = d))
}
# family parameter at its limit with a random effect, Laplace vs
# quadrature (quadseam merge)
if (ok("nbquad")) {
  set.seed(2)
  d <- data.frame(x = rnorm(200), g = factor(rep(1:20, each = 10)))
  d$y <- rpois(200, exp(1 + 0.4 * d$x + rnorm(20, 0, 0.5)[d$g]))
  show("nbquad laplace", F(y ~ x + (1 | g), family = negbinomial(), data = d))
  show("nbquad quadrature", F(y ~ x + (1 | g), family = negbinomial(), data = d, quadrature = TRUE))
}
# what predict(), fitted() and confint() give on limit fits
if (ok("pred")) {
  set.seed(2)
  d <- data.frame(x = rnorm(200))
  d$y <- rpois(200, exp(1 + 0.4 * d$x))
  r <- F(y ~ x, family = negbinomial(), data = d)
  fit <- r$value
  nd <- d[1:3, , drop = FALSE]
  pr <- conds(predict(fit, newdata = nd, ndraws = 400))
  cat("  pred nb: W", length(pr$w), "M", length(pr$m), "\n")
  for (x in c(pr$w, pr$m)) cat("    ", substr(gsub("\n", " ", x), 1, 300), "\n")
  print(pr$value)
  fi <- conds(fitted(fit, newdata = nd))
  cat("  fitted nb:\n"); print(fi$value)
  for (x in c(fi$w, fi$m)) cat("    ", substr(gsub("\n", " ", x), 1, 300), "\n")
  ci <- conds(confint(fit))
  print(ci$value)
  for (x in c(ci$w, ci$m)) cat("    ", substr(gsub("\n", " ", x), 1, 300), "\n")
  set.seed(2)
  d0 <- data.frame(x = rnorm(50), y = 0L)
  r0 <- F(y ~ x, family = poisson(), data = d0)
  fi <- conds(fitted(r0$value, newdata = d0[1:3, , drop = FALSE]))
  cat("  fitted pois0:\n"); print(fi$value)
  for (x in c(fi$w, fi$m)) cat("    ", substr(gsub("\n", " ", x), 1, 300), "\n")
  pr <- conds(predict(r0$value, newdata = d0[1:3, , drop = FALSE], ndraws = 400))
  cat("  predict pois0:\n"); print(pr$value)
  for (x in c(pr$w, pr$m)) cat("    ", substr(gsub("\n", " ", x), 1, 300), "\n")
}

# Review of lane sampfront: what brms 2.23.0 does with
# drop_unused_levels = FALSE and a grouping factor with an unused level,
# against frm() and frm_sample(). No Stan compile: make_standata() and
# make_stancode().
#   Rscript dev/review-sampfront-dul.R
.libPaths(c("C:/Users/adf44/source/r/wt-sampfront-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(11)
d <- data.frame(g = factor(rep(1:6, each = 8), levels = 1:7), x = rnorm(48),
                f = factor(rep(c("a", "b"), 24), levels = c("a", "b", "c")))
d$y <- 1 + 0.5 * d$x + rnorm(6)[d$g] * 0.7 + rnorm(48)
for (dul in c(TRUE, FALSE)) {
  sd <- make_standata(y ~ x + (1 | g), data = d, drop_unused_levels = dul)
  cat("\nbrms drop_unused_levels =", dul, ": N_1 (levels of g) =", sd$N_1,
      " max J_1 =", max(sd$J_1), "\n")
  sf <- make_standata(y ~ f, data = d, drop_unused_levels = dul)
  cat("  y ~ f: K =", sf$K, " colnames X:", colnames(sf$X), "\n")
}
suppressPackageStartupMessages(library(frmtmb))
for (dul in c(TRUE, FALSE)) {
  f <- suppressWarnings(suppressMessages(
    frm(y ~ x + (1 | g), data = d, drop_unused_levels = dul)))
  cat("\nfrm drop_unused_levels =", dul, ": nlevels(fit$data$g) =",
      nlevels(f$data$g), " length(b) =", length(f$estimates$b),
      " ranef rows =", nrow(ranef(f)$g), "\n")
  ff <- suppressWarnings(suppressMessages(
    frm(y ~ f, data = d, drop_unused_levels = dul)))
  cat("  y ~ f: coefficients:", names(fixef(ff)[, 1]), "\n")
}
cat("\nDONE\n")

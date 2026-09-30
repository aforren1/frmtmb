# Punch round 1, user decision: the conditional_effects() bootstrap
# holds every smooth. Boot width over Wald width on y ~ s(x) and on
# y ~ s(x) + (1 | g) at re_formula NA, NULL and an observed g, on the
# lane build or the base build.
#   Rscript dev/postfit2-p1-smooth.R <lane|base>
# Seeds: data 24 (s(x)) and 26 (s(x) + (1 | g)), boot 60 seed 5.
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/wt-postfit2-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "base")) libs <- libs[-1]
.libPaths(libs)
suppressMessages(library(frmtmb))
say <- function(...) cat(sprintf(...), "\n", sep = "")
f4 <- function(v) paste(sprintf("%.3f", v), collapse = " ")
say("ARM %s: frmtmb from %s", arm, find.package("frmtmb"))
wd <- function(d) d$upper__ - d$lower__
one <- function(label, fit, ...) {
  w <- suppressMessages(conditional_effects(fit, "x", resolution = 5,
                                            ...))$x
  b <- suppressMessages(conditional_effects(fit, "x", resolution = 5,
                                            band = "boot", boot = 60,
                                            seed = 5, ...))$x
  say("%s: boot/wald width %s; estimate inside the boot band %s", label,
      f4(wd(b) / wd(w)),
      all(b$lower__ <= w$estimate__ & w$estimate__ <= b$upper__))
}
set.seed(24)
d <- data.frame(x = runif(100))
d$y <- sin(2 * pi * d$x) + rnorm(100, 0, 0.3)
f1 <- frm(bf(y ~ s(x)), family = gaussian(), data = d)
one("s(x), re_formula NA", f1)
one("s(x), re_formula NULL", f1, re_formula = NULL)
set.seed(26)
d2 <- data.frame(x = runif(160), g = factor(rep(1:8, 20)))
d2$y <- sin(2 * pi * d2$x) + rnorm(8, 0, 1)[d2$g] + rnorm(160, 0, 0.3)
f2 <- frm(bf(y ~ s(x) + (1 | g)), family = gaussian(), data = d2)
one("s(x) + (1 | g), re_formula NA", f2)
one("s(x) + (1 | g), re_formula NULL (new g)", f2, re_formula = NULL)
one("s(x) + (1 | g), observed g = 3", f2, re_formula = NULL,
    conditions = list(g = factor("3", levels = 1:8)))

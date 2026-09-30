# conditional_effects(re_formula = NULL, conditions = <an OBSERVED level>)
# must read that level. Construction from lane sampfix's
# dev/sampfix-08-ce-level.R (its tree, read only), extended to the fit
# method's band = "boot" and to the spaghetti path.
#   POSTFIT2_LIB=base Rscript dev/postfit2-celevel.R   (the base build)
#   Rscript dev/postfit2-celevel.R                      (this lane)
lib <- Sys.getenv("POSTFIT2_LIB", "lane")
libs <- c("C:/Users/adf44/source/r/wt-postfit2-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (lib == "base") libs <- libs[-1]
.libPaths(libs)
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
say <- function(...) cat(sprintf(...), "\n", sep = "")
say("arm %s: frmtmb %s, frmtmb.sample %s", lib,
    format(packageVersion("frmtmb")), format(packageVersion("frmtmb.sample")))
set.seed(9)
dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
lab <- frmtmb::brms_par_labels(fit)
tpl <- fit$frame$par_template
est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
set.seed(1)
M <- matrix(rep(est, each = 5) + rnorm(5 * length(est), 0, 0.05), 5,
            dimnames = list(NULL, lab))
M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
full <- structure(list(stanfit = NULL, draws = M, fit = fit),
                  class = "frmtmb_draws")
ce <- function(o, lev, ...) {
  conditional_effects(o, effects = "x", resolution = 3, re_formula = NULL,
                      conditions = data.frame(g = factor(lev, levels = 1:6)),
                      seed = 4, ...)$x
}
for (lev in 1:2) {
  shifted <- full
  col <- sprintf("r_g[%d,Intercept]", lev)
  shifted$draws[, col] <- shifted$draws[, col] + 10
  a <- ce(full, lev)$estimate__
  b <- ce(shifted, lev)$estimate__
  say("draws, g = %d: curve moves by %s when r_g[%d,] moves by 10",
      lev, paste(format(b - a, digits = 6), collapse = " "), lev)
}
say("draws: level 1 minus level 2: %s (the median of r_g[1] - r_g[2]: %.5f)",
    paste(format(ce(full, 1)$estimate__ - ce(full, 2)$estimate__,
                 digits = 6), collapse = " "),
    median(M[, "r_g[1,Intercept]"] - M[, "r_g[2,Intercept]"]))
sp <- tryCatch({
  s1 <- attr(ce(full, 2, spaghetti = TRUE), "spaghetti")
  sh <- full
  sh$draws[, "r_g[2,Intercept]"] <- sh$draws[, "r_g[2,Intercept]"] + 10
  s2 <- attr(ce(sh, 2, spaghetti = TRUE), "spaghetti")
  sprintf("spaghetti moves by %s", paste(format(range(s2$estimate__ -
                                                        s1$estimate__),
                                                 digits = 6),
                                          collapse = " to "))
}, error = function(e) paste("spaghetti:", conditionMessage(e)))
say("draws, g = 2: %s", sp)
# the fit method's bootstrap: the refits' curves at an observed level
# centre on that level's curve, not on the population one
cf <- function(lev, band, ...) {
  conditional_effects(fit, effects = "x", resolution = 3, re_formula = NULL,
                      conditions = data.frame(g = factor(lev, levels = 1:6)),
                      band = band, ...)
}
rg <- ranef(fit)$g
re <- if (length(dim(rg)) == 3L) rg[, 1, 1] else rg[, 1]
lev <- which.max(abs(re))
w <- cf(lev, "wald")$x
bt <- cf(lev, "boot", boot = 40, seed = 3)
m <- colMeans(attr(bt, "boot")$t, na.rm = TRUE)
say("fit, g = %d (mode %.4f): wald estimate %s", lev, re[lev],
    paste(format(w$estimate__, digits = 5), collapse = " "))
say("fit, boot mean of the refit curves %s; its offset from the estimate over the mode: %.3f",
    paste(format(m, digits = 5), collapse = " "),
    max(abs(m - w$estimate__)) / abs(re[lev]))
say("fit, boot band width / wald band width: %.3f",
    mean(bt$x$upper__ - bt$x$lower__) / mean(w$upper__ - w$lower__))
b1 <- cf(1, "boot", boot = 40, seed = 3)$x
b2 <- cf(2, "boot", boot = 40, seed = 3)$x
say("fit, boot band at g = 1 minus at g = 2: lower %s (estimates differ by %s)",
    paste(format(b1$lower__ - b2$lower__, digits = 4), collapse = " "),
    paste(format(b1$estimate__ - b2$estimate__, digits = 4),
          collapse = " "))
# the same with group effects large enough to see (group SD 2)
set.seed(21)
d2 <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
d2$y <- rnorm(80, 1 + 0.5 * d2$x + rnorm(8, 0, 2)[d2$g], 1)
f2 <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d2)
rg2 <- ranef(f2)$g
re2 <- if (length(dim(rg2)) == 3L) rg2[, 1, 1] else rg2[, 1]
lv <- which.max(abs(re2))
cnd <- data.frame(g = factor(lv, levels = 1:8))
w2 <- conditional_effects(f2, "x", resolution = 3, re_formula = NULL,
                          conditions = cnd)$x
bb <- conditional_effects(f2, "x", resolution = 3, re_formula = NULL,
                          conditions = cnd, band = "boot", boot = 60,
                          seed = 5)
m2 <- colMeans(attr(bb, "boot")$t, na.rm = TRUE)
say("fit (group SD 2, seed 21), g = %d, mode %.3f: boot mean minus estimate %s; over the mode: %.3f",
    lv, re2[lv], paste(format(m2 - w2$estimate__, digits = 4), collapse = " "),
    max(abs(m2 - w2$estimate__)) / abs(re2[lv]))
say("  estimate inside the boot band at every point: %s; band width / wald width %.3f",
    all(bb$x$lower__ <= w2$estimate__ & w2$estimate__ <= bb$x$upper__),
    mean(bb$x$upper__ - bb$x$lower__) / mean(w2$upper__ - w2$lower__))

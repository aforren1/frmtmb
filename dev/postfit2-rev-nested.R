# Reviewer: the +/-1410 Wald band of the nested fit in
# dev/postfit2-rev-adv-fit.R (section C), on both arms. Is it the fit
# (g and g:h confounded in that design) or the new-level variance?
# A second, properly nested design (every g crossed with 4 h) is the
# control.
#   Rscript dev/postfit2-rev-nested.R <base|lane>
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "lane"
libs <- c("C:/Users/adf44/source/r/wt-postfit2-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
suppressMessages(library(frmtmb))
say <- function(...) cat(sprintf(...), "\n", sep = "")
f3 <- function(v) paste(sprintf("%.4f", v), collapse = " ")
say("ARM %s", arm)
one <- function(label, dn) {
  fn <- frm(bf(y ~ x + (1 | g / h)), family = gaussian(), data = dn)
  vc <- VarCorr(fn)
  for (nm in names(vc)) say("%s VarCorr %s sd %s", label, nm,
                            f3(vc[[nm]]$sd[, 1]))
  say("%s sdr of theta: %s", label, f3(sqrt(diag(vcov(fn, full = TRUE)))[
    grep("theta", rownames(vcov(fn, full = TRUE)))]))
  rg <- ranef(fn)$g
  rg <- if (length(dim(rg)) == 3L) rg[, 1, 1] else rg[, 1]
  lv <- which.max(abs(rg))
  cn <- data.frame(g = factor(lv, levels = levels(dn$g)))
  w <- conditional_effects(fn, "x", resolution = 3, re_formula = NULL,
                           conditions = cn)$x
  say("%s g=%d set, g:h new: wald est %s lo %s hi %s", label, lv,
      f3(w$estimate__), f3(w$lower__), f3(w$upper__))
  w0 <- conditional_effects(fn, "x", resolution = 3, re_formula = NULL)$x
  say("%s nothing set: wald lo %s hi %s", label, f3(w0$lower__),
      f3(w0$upper__))
  bo <- suppressMessages(conditional_effects(
    fn, "x", resolution = 3, re_formula = NULL, conditions = cn,
    band = "boot", boot = 60, seed = 5))$x
  say("%s g=%d set, g:h new: boot60 seed5 lo %s hi %s; width/wald %s",
      label, lv, f3(bo$lower__), f3(bo$upper__),
      f3((bo$upper__ - bo$lower__) / (w$upper__ - w$lower__)))
  cn2 <- data.frame(g = factor(lv, levels = levels(dn$g)),
                    h = factor(dn$h[dn$g == lv][1], levels = levels(dn$h)))
  w2 <- conditional_effects(fn, "x", resolution = 3, re_formula = NULL,
                            conditions = cn2)$x
  say("%s g and h both set (observed g:h): wald lo %s hi %s", label,
      f3(w2$lower__), f3(w2$upper__))
}
set.seed(23)
dn <- data.frame(x = rnorm(160), g = factor(rep(1:8, 20)),
                 h = factor(rep(1:4, each = 2, length.out = 160)))
dn$y <- rnorm(160, 1 + 0.5 * dn$x + rnorm(8, 0, 2)[dn$g] +
                rnorm(32, 0, 0.7)[interaction(dn$g, dn$h)], 1)
say("confounded design: g:h levels %d, g levels %d",
    nlevels(droplevels(interaction(dn$g, dn$h))), nlevels(dn$g))
one("C-confounded", dn)
set.seed(25)
dm <- expand.grid(rep = 1:5, h = factor(1:4), g = factor(1:8))
dm$x <- rnorm(nrow(dm))
dm$y <- rnorm(nrow(dm), 1 + 0.5 * dm$x + rnorm(8, 0, 2)[dm$g] +
                rnorm(32, 0, 0.7)[interaction(dm$g, dm$h)], 1)
say("crossed-within design: g:h levels %d",
    nlevels(droplevels(interaction(dm$g, dm$h))))
one("C-nested", dm)

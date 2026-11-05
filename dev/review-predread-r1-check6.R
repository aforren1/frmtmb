# Short check: se__ of conditional_effects() on families with no me() term: the lane's
# fix moved lognormal and hurdle. Against brms's sampled draws. Seed 7801's recipe.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7801)
n <- 240
d <- data.frame(x = rnorm(n))
eta <- 0.4 + 0.6 * d$x
d$cnt <- rpois(n, exp(eta))
d$pos <- exp(eta + rnorm(n, 0, 0.4))
d$hp <- ifelse(runif(n) < 0.25, 0, d$pos)
for (sp in list(c("bf(cnt ~ x)", "poisson()"), c("bf(pos ~ x)", "lognormal()"), c("bf(hp ~ x)", "hurdle_lognormal()"))) {
  cat("\n===", sp, "\n")
  ff <- quiet(frm(rev_ev(sp[1], "frmtmb"), family = rev_ev(sp[2], "frmtmb"), data = d))
  bb <- quiet(brms::brm(rev_ev(sp[1], "brms"), family = rev_ev(sp[2], "brms"), data = d, chains = 2, iter = 1500,
                        refresh = 0, seed = 3, silent = 2))
  cf <- conditional_effects(ff, "x", resolution = 4)[[1]]
  cb <- brms::conditional_effects(bb, "x", resolution = 4)[[1]]
  print(signif(data.frame(x = cf$x, frm_est = cf$estimate__, frm_se = cf$se__, frm_lo = cf$lower__, frm_hi = cf$upper__,
                          brms_est = cb$estimate__, brms_se = cb$se__, brms_lo = cb$lower__, brms_hi = cb$upper__), 4))
  fv <- fitted(ff, newdata = data.frame(x = cf$x))
  fb <- fitted(bb, newdata = data.frame(x = cf$x))
  print(signif(cbind(fv, fb), 4))
}

# Item 4: on the fitted rows nothing of a me() fit may change: digit-for-digit against
# rellib-r16, one process per build. Also the difference curve on both builds.
#   REV_LIB=before Rscript dev/review-predread-regress-me.R <out.rds>
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
outf <- commandArgs(TRUE)[1]
set.seed(7402)
n <- 200
d <- data.frame(z = rnorm(n), lat = rnorm(n), sx = 0.4, g = factor(rep(sprintf("g%02d", 1:20), each = 10)))
d$xo <- d$lat + rnorm(n, 0, 0.4)
gl <- rnorm(20)[d$g]
d$sg <- 0.4
d$xg <- gl + rnorm(20, 0, 0.4)[d$g]
d$y <- 1 + 1.2 * d$lat - 0.5 * d$z + rnorm(n, 0, 0.5)
d$yg <- 1 + 1.2 * gl - 0.5 * d$z + rnorm(n, 0, 0.5)
d$cnt <- rpois(n, exp(0.5 + 0.8 * d$lat))
fits <- list(
  plain = quiet(frm(bf(y ~ me(xo, sx) + z), family = gaussian(), data = d)),
  gr = quiet(frm(bf(yg ~ me(xg, sg, gr = g) + z), family = gaussian(), data = d)),
  log = quiet(frm(bf(cnt ~ me(xo, sx)), family = poisson(), data = d)))
out <- list()
for (nm in names(fits)) {
  f <- fits[[nm]]
  ds <- rev_draws_at(f, 3L)
  out[[paste0(nm, "|fitted")]] <- fitted(f)
  out[[paste0(nm, "|linpred")]] <- frm_linpred(f, se.fit = TRUE)
  out[[paste0(nm, "|predict")]] <- {set.seed(3); predict(f, ndraws = 40)}
  out[[paste0(nm, "|predict_nop")]] <- {set.seed(3); predict(f, ndraws = 40, propagate_error = FALSE)}
  out[[paste0(nm, "|simulate")]] <- as.matrix(simulate(f, nsim = 3, seed = 3))
  out[[paste0(nm, "|residuals")]] <- residuals(f)
  out[[paste0(nm, "|basis")]] <- frm_lp_basis(f)[c("eta", "extra_var")]
  out[[paste0(nm, "|draws_epred")]] <- frmtmb.sample::posterior_epred(ds)
  out[[paste0(nm, "|draws_pp")]] <- {set.seed(4); frmtmb.sample::posterior_predict(ds)}
}
saveRDS(out, outf)
g <- data.frame(xo = c(1, 1), sx = 0.6, z = c(0, 1))
h <- data.frame(xo = c(1, 1), sx = 0.6, z = c(1, 0))
cu <- frm_curve(fits$plain, newdata = g, contrast = h, simultaneous = FALSE)
cat("frm_curve contrast in z, same xo, sx: est", signif(cu$.estimate, 5), " se", signif(cu$.se, 5), "\n")
cat("coef z se:", signif(summary(fits$plain)$fixed["z", 2], 5), "\n")
em <- as.data.frame(emmeans::contrast(emmeans::emmeans(fits$plain, ~ z, at = list(z = c(0, 1))), "revpairwise"))
cat("emmeans contrast in z:", signif(em$estimate, 5), " SE", signif(em$SE, 5), "\n")
gg <- data.frame(xg = c(1, 1), sg = 0.6, z = c(0, 1), g = "g01")
hh <- data.frame(xg = c(1, 1), sg = 0.6, z = c(1, 0), g = "g01")
cg <- frm_curve(fits$gr, newdata = gg, contrast = hh, simultaneous = FALSE)
cat("me(gr =), same level: est", signif(cg$.estimate, 5), " se", signif(cg$.se, 5), " coef z se:",
    signif(summary(fits$gr)$fixed["z", 2], 5), "\n")

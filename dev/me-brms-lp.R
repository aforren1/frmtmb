# me() against brms 2.23.0's own log density (lane me, 2026-09-25).
#
# For each shape: fit with frmtmb, translate the estimates and the
# latent modes into brms's Stan parameters (stan_pars_from_fit() in
# tests/testthat/helper-brms.R), and compare rstan::log_prob() of brms's
# program with flat priors against frmtmb's joint density plus the
# log-Jacobian of the map. Then the gradient of brms's log density on
# the zme block, which vanishes when the latent modes are brms's modes.
#
# Run from the worktree root, one process, Stan compiles each shape once:
#   FRMTMB_STAN_CACHE=/tmp/lanes/me/stan-cache \
#     Rscript dev/me-brms-lp.R
# The output is dev/me-brms-lp.txt.

lib <- Sys.getenv("FRMTMB_LIB", "/opt/rlib/lane-me")
.libPaths(c(lib, "/opt/rlib/base", "/opt/rlib/deps", "/opt/r/lib/R/library"))
suppressMessages({
  library(frmtmb)
  library(testthat)
})
# the helpers call package internals, as they do inside test_file()
env <- new.env(parent = asNamespace("frmtmb"))
sys.source("tests/testthat/helper-brms.R", envir = env)
options(frmtmb.brms_lp_report = TRUE)

set.seed(23)
n <- 80
tx <- rnorm(n, 1, 0.8)
tz <- 0.5 * tx + rnorm(n, 0, 0.7)
d <- data.frame(x = tx + rnorm(n, 0, 0.3), sx = runif(n, 0.2, 0.4),
                z = tz + rnorm(n, 0, 0.3), sz = 0.3, w = rnorm(n),
                g = factor(rep(1:16, each = 5)))
d$y <- 2 + 0.7 * tx - 0.4 * tz + 0.2 * d$w + rnorm(n, 0, 0.5)
d$xg <- rep(rnorm(16), each = 5)
d$sxg <- rep(runif(16, 0.2, 0.4), each = 5)
d$cnt <- rpois(n, exp(0.2 + 0.4 * tx))

cases <- list(
  cor_interaction = list(brm = brms::bf(y ~ me(x, sx) * me(z, sz) + w),
                         frm = bf(y ~ me(x, sx) * me(z, sz) + w),
                         fam = gaussian()),
  mecor_false = list(brm = brms::bf(y ~ me(x, sx) + me(z, sz)) +
                       brms::set_mecor(FALSE),
                     frm = bf(y ~ me(x, sx) + me(z, sz)) + set_mecor(FALSE),
                     fam = gaussian()),
  sigma_dpar = list(brm = brms::bf(y ~ w, sigma ~ me(x, sx)),
                    frm = bf(y ~ w, sigma ~ me(x, sx)), fam = gaussian()),
  gr_level = list(brm = brms::bf(y ~ me(x, sx) + me(xg, sxg, gr = g)),
                  frm = bf(y ~ me(x, sx) + me(xg, sxg, gr = g)),
                  fam = gaussian()),
  poisson = list(brm = brms::bf(cnt ~ me(x, sx)), frm = bf(cnt ~ me(x, sx)),
                 fam = poisson()))

for (nm in names(cases)) {
  cs <- cases[[nm]]
  fit <- frm(cs$frm + cs$fam, data = d)
  res <- env$with_brms_me(function() {
    env$brms_lp_check(cs$brm, cs$fam, d, fit, joint = TRUE)
  })
  cat(sprintf(paste0("ME-LP %-16s brms lp %.10f  frmtmb joint + logJ ",
                     "%.10f  diff %.3e  rel %.3e  max|grad zme| %.3e\n"),
              nm, res$lp, res$ours, res$measured_const,
              abs(res$measured_const) / abs(res$ours), res$max_grad))
}

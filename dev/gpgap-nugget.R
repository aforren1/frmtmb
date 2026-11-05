# Once the positions are in brms's order, is the rest of the exact gp()
# gap the one line where the nuggets differ?
#
#   Rscript dev/gpgap-nugget.R "<spec key>"
#
# brms's program is compiled twice: as brms writes it, and with its
# nugget line, cov[n, n] += 1e-12, replaced by frmtmb's relative nugget,
# cov[n, n] += 1e-6 * square(sdgp). Both are read at the same z, taken
# through frmtmb's kernel in brms's order. At frmtmb's nugget the kernel
# is conditioned about 1e7, so double precision reads it: a constant of
# 0 says every other piece of the two models (scaling, kernel, latent
# map, prior of z, likelihood) is the same. Base build unless GPGAP_LIB.
args <- commandArgs(trailingOnly = TRUE)
key <- if (length(args)) args[1] else "lognormal x gp"
source("dev/gpgap-lib.R")
E$key <- key
local(envir = E, {
  s <- gpgap_setup(eq_specs[[key]]())
  sf0 <- gpgap_sf(s$code, s$sdat)
  sf1 <- gpgap_sf(gpgap_frm_nugget(s$code), s$sdat)
  L <- t(chol(s$kb))
  z <- as.numeric(solve(L, s$fb))
  p <- s$pars
  p[[s$zn]] <- array(z, length(z))
  logJ <- s$logJ - sum(log(diag(t(chol(s$kf))))) + sum(log(diag(L)))
  ours <- s$ours0 + logJ
  c0 <- gpgap_lp(sf0, p) - ours
  c1 <- gpgap_lp(sf1, p) - ours
  u1 <- rstan::unconstrain_pars(sf1, p)
  g1 <- rstan::grad_log_prob(sf1, u1, adjust_transform = FALSE)
  idx <- brms_inner_index(sf1, p)
  ev <- eigen(s$kb / s$sdgp^2, symmetric = TRUE, only.values = TRUE)$values
  cat(sprintf(paste0("GPNUG %s: ours %.4f, tol 1e-6*|ours| %.3g, cond(K+nug) ",
                     "%.3g\n",
                     "  brms's program, 1e-12        constant %.6f\n",
                     "  frmtmb's nugget in that line constant %.3e\n",
                     "  gradient in zgp at frmtmb's nugget, max abs %.3e\n"),
              key, ours, 1e-6 * abs(ours), max(ev) / min(ev), c0, c1,
              max(abs(g1[idx]))))
})

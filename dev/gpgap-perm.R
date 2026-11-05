# Is the exact gp() L1 gap the order of the latent positions?
#
#   Rscript dev/gpgap-perm.R "<spec key>"
#
# frmtmb sorts a term's distinct positions (R/frame.R, posdf[order]);
# brms keeps them in order of first appearance (Jgp_1 = 1 2 3 ...). The
# suite's translator takes zgp as L_frm^-1 f in frmtmb's order and hands
# it to brms, which multiplies by L_brms in brms's order. On x that is
# not already sorted that puts each latent value at another position.
#
# Reads at frmtmb's estimate, base build rellib-r16 (no model change):
#  A   the instrument as it is, stan_pars_from_fit();
#  A1  z through frmtmb's own kernel, permuted to brms's order;
#  and for A1 the likelihood at Stan's rebuilt field against R's reading
#  of it (lp - std_normal(z) - ll(f)), which is the whole of A1 if the
#  prior of z cancels the frmtmb prior of f.
args <- commandArgs(trailingOnly = TRUE)
key <- if (length(args)) args[1] else "lognormal x gp"
.libPaths(c("C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
E$key <- key
local(envir = E, {
  spec <- eq_specs[[key]]()
  fit <- eq_frm_fit(spec)
  bform <- eq_eval(spec$bf, "brms")
  ft <- eq_fam_text(spec, "brms")
  family <- if (is.na(ft)) NULL else eq_eval(ft, "brms")
  prior <- eq_brms_prior(spec, bform, family, spec$data)
  args <- eq_brms_args(spec, spec$data, prior)
  code <- do.call(brms::make_stancode, c(list(bform, family = family), args))
  sdat <- suppressMessages(do.call(brms::make_standata,
                                   c(list(bform, family = family), args)))
  sf <- suppressMessages(rstan::sampling(brms_stan_model(code), data = sdat,
                                         chains = 0))
  pars <- stan_pars_from_fit(fit, sdat, code)
  lp <- function(p) {
    rstan::log_prob(sf, rstan::unconstrain_pars(sf, p),
                    adjust_transform = FALSE)
  }
  e <- fit$obj$env
  ours0 <- -e$f(e$last.par.best)
  A <- lp(pars) - (ours0 + attr(pars, "logJ"))

  zn <- grep("^zgp_", names(pars), value = TRUE)
  stopifnot(length(zn) == 1L)
  ip <- brms_idx_parts(zn, "zgp")
  dp <- ip$dpar
  bk <- brms_gp_block(fit, dp, ip$idx[[1]])
  gsfx <- sub("^zgp_", "", zn)
  xb <- as.matrix(sdat[[paste0("Xgp_", gsfx)]])
  jb <- as.integer(sdat[[paste0("Jgp_", gsfx)]])
  kf <- brms_block_cov(fit, bk)
  fv <- as.numeric(brms_block_b(fit, bk))
  # brms's position of frmtmb's i-th sorted position; dividing by dmax
  # keeps the order, so sorting brms's scaled inputs lands on frmtmb's
  ob <- do.call(order, as.data.frame(xb))
  fb <- numeric(length(fv))
  fb[ob] <- fv
  kb <- matrix(0, nrow(kf), ncol(kf))
  kb[ob, ob] <- kf
  lp1 <- t(chol(kb))
  z1 <- as.numeric(solve(lp1, fb))
  p1 <- pars
  p1[[zn]] <- array(z1, length(z1))
  logJ0 <- sum(log(diag(t(chol(kf)))))
  logJ1 <- sum(log(diag(lp1)))
  A1 <- lp(p1) - (ours0 + attr(pars, "logJ") - logJ0 + logJ1)
  th <- fit$estimates$theta[bk$theta_idx]
  sdgp <- exp(th[[1]])
  ev <- eigen(kf / sdgp^2, symmetric = TRUE, only.values = TRUE)$values
  # the kernel in brms's order and scale against frmtmb's, at the fit's
  # lscale, off the diagonal
  rho <- as.numeric(pars[[paste0("lscale_", gsfx)]])
  d2 <- as.matrix(stats::dist(xb))^2
  kstan <- sdgp^2 * exp(-d2 / (2 * rho^2))
  offd <- max(abs((kstan - kb)[upper.tri(kb)]))
  cat(sprintf(paste0("GPPERM %s: positions %d, sdgp %.4g, rho %.4g; ",
                     "eig<1e-6 %d, eig<1e-12 %d\n",
                     "  sorted x already: %s\n",
                     "  A  translator as is        %.6f\n",
                     "  A1 z through frm kernel, brms order  %.6f\n",
                     "  logJ0 - logJ1 (rounding of permuted det)  %.3g\n",
                     "  kernel off-diagonal, brms scale and order: %.3g\n"),
              key, length(fv), sdgp, rho, sum(ev < 1e-6), sum(ev < 1e-12),
              !is.unsorted(ob), A, A1, logJ0 - logJ1, offd))
})

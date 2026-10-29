# Where does the exact gp() L1 gap come from?
#
#   Rscript dev/equiv-gp-gap.R ["<feature_a> x <feature_b>"]
#
# The suite's translator standardizes frmtmb's field f through frmtmb's
# OWN kernel, z = L_frm^-1 f, and brms rebuilds the field through ITS
# kernel, L_brms z. On a kernel with eigenvalues near both nuggets the
# rebuilt field is not f, so the likelihood term moves too. This script
# reads the pair three ways at frmtmb's estimate:
#   A  the instrument's L1 constant (z through L_frm);
#   B  the constant with z taken through brms's kernel, z = L_brms^-1 f,
#      and the map's Jacobian log|L_brms|, so brms rebuilds f exactly;
#   C  the gap the two kernels alone predict for that same f,
#      log N(f | 0, K_brms) - log N(f | 0, K_frm).
# If B equals C, the whole gap is the two nuggets acting on the prior of
# the field and the model is otherwise the same function. Writes nothing.
args <- commandArgs(trailingOnly = TRUE)
key <- if (length(args)) args[1] else "lognormal x gp"
.libPaths(c("C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/review-equiv-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
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

  dp <- if (grepl("mvbf", spec$bf)) "y1_mu" else "mu"
  bk <- brms_gp_block(fit, dp, 1)
  kf <- brms_block_cov(fit, bk)
  fv <- as.numeric(brms_block_b(fit, bk))
  th <- fit$estimates$theta[bk$theta_idx]
  sdgp <- exp(th[[1]])
  kb <- kf
  diag(kb) <- diag(kb) - frmtmb:::gp_nugget * sdgp^2 + 1e-12
  lf <- t(chol(kf))
  lb <- t(chol(kb))
  zn <- grep("^zgp_", names(pars), value = TRUE)
  p2 <- pars
  p2[[zn]] <- array(as.numeric(solve(lb, fv)), length(fv))
  logJ2 <- attr(pars, "logJ") - sum(log(diag(lf))) + sum(log(diag(lb)))
  B <- lp(p2) - (ours0 + logJ2)
  ldm <- function(f, L) {
    -0.5 * sum(solve(L, f)^2) - sum(log(diag(L))) - 0.5 * length(f) *
      log(2 * pi)
  }
  C <- ldm(fv, lb) - ldm(fv, lf)
  ev <- eigen(kf / sdgp^2, symmetric = TRUE, only.values = TRUE)$values
  cat(sprintf(paste0("GPGAP %s: positions %d, sdgp %.4g; eigenvalues of ",
                     "the correlation below 1e-6: %d, below 1e-12: %d\n",
                     "  A instrument constant  %.6f\n",
                     "  B z through brms kernel %.6f\n",
                     "  C nugget-only gap      %.6f\n",
                     "  B - C                  %.3g\n"),
              key, length(fv), sdgp, sum(ev < 1e-6), sum(ev < 1e-12),
              A, B, C, B - C))
})

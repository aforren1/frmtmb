# Shared setup of the gpgap scripts: one pair's fit, brms's program and
# data, and the block read in brms's order of positions. The libraries
# come from GPGAP_LIB (`;`-separated), so a script can read the base
# build or the lane's.
gpgap_libs <- strsplit(Sys.getenv("GPGAP_LIB",
                                  "C:/Users/adf44/source/r/rellib-r16"),
                       ";")[[1]]
.libPaths(c(gpgap_libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
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

# The Stan program as brms writes it, or with one line replaced, on the
# same data.
local(envir = E, {
  gpgap_setup <- function(spec) {
    fit <- eq_frm_fit(spec)
    bform <- eq_eval(spec$bf, "brms")
    ft <- eq_fam_text(spec, "brms")
    family <- if (is.na(ft)) NULL else eq_eval(ft, "brms")
    prior <- eq_brms_prior(spec, bform, family, spec$data)
    args <- eq_brms_args(spec, spec$data, prior)
    code <- do.call(brms::make_stancode,
                    c(list(bform, family = family), args))
    sdat <- suppressMessages(do.call(brms::make_standata,
                                     c(list(bform, family = family), args)))
    pars <- stan_pars_from_fit(fit, sdat, code)
    zn <- grep("^zgp_", names(pars), value = TRUE)
    stopifnot(length(zn) == 1L)
    ip <- brms_idx_parts(zn, "zgp")
    bk <- brms_gp_block(fit, ip$dpar, ip$idx[[1]])
    gsfx <- sub("^zgp_", "", zn)
    xb <- as.matrix(sdat[[paste0("Xgp_", gsfx)]])
    ob <- do.call(order, as.data.frame(xb))
    kf <- brms_block_cov(fit, bk)
    fv <- as.numeric(brms_block_b(fit, bk))
    fb <- numeric(length(fv))
    fb[ob] <- fv
    kb <- matrix(0, nrow(kf), ncol(kf))
    kb[ob, ob] <- kf
    th <- fit$estimates$theta[bk$theta_idx]
    e <- fit$obj$env
    list(fit = fit, code = code, sdat = sdat, pars = pars, zn = zn, bk = bk,
         gsfx = gsfx, xb = xb, ob = ob, kf = kf, kb = kb, fv = fv, fb = fb,
         sdgp = exp(th[[1]]),
         rho = as.numeric(pars[[paste0("lscale_", gsfx)]]),
         ours0 = -e$f(e$last.par.best), logJ = attr(pars, "logJ"))
  }
  gpgap_sf <- function(code, sdat) {
    suppressMessages(rstan::sampling(brms_stan_model(code), data = sdat,
                                     chains = 0))
  }
  gpgap_lp <- function(sf, p) {
    rstan::log_prob(sf, rstan::unconstrain_pars(sf, p),
                    adjust_transform = FALSE)
  }
  # brms's program with frmtmb's nugget in brms's one line
  gpgap_frm_nugget <- function(code) {
    new <- gsub("cov[n, n] += 1e-12;", "cov[n, n] += 1e-6 * square(sdgp);",
                code, fixed = TRUE)
    stopifnot(!identical(new, code))
    new
  }
})

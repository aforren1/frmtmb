# The five survey pairs with an exact gp(): brms's own line against the line
# rewritten to frmtmb's relative one. The rewrite is attribution only. Also the
# closed form's cost where one exists: the marginal at brms's absolute line at frm()'s
# maximum, against its own maximum. Spec data, seeds in dev/equiv-specs.R.
#   Rscript dev/ncgp-l1eq.R <key>
key <- commandArgs(TRUE)[1]
LIB <- "C:/Users/adf44/source/r/wt-ncgp-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) sys.source(f, envir = E)
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
  sdat <- suppressMessages(do.call(brms::make_standata, c(list(bform, family = family), args)))
  pars <- stan_pars_from_fit(fit, sdat, code)
  e <- fit$obj$env
  ours <- -as.numeric(e$f(e$last.par.best)) + attr(pars, "logJ")
  read <- function(cd) {
    sf <- suppressMessages(rstan::sampling(brms_stan_model(cd), data = sdat, chains = 0))
    up <- rstan::unconstrain_pars(sf, pars)
    lp <- rstan::log_prob(sf, up, adjust_transform = FALSE)
    g <- rstan::grad_log_prob(sf, up, adjust_transform = FALSE)
    idx <- brms_inner_index(sf, pars)
    c(lp - ours, max(abs(g[idx])))
  }
  a <- read(code)
  b <- read(gsub("cov[n, n] += 1e-12;", "cov[n, n] += 1e-12 * square(sdgp);", code, fixed = TRUE))
  bk <- Filter(function(k) identical(k$covstruct, "gp"), fit$frame$re_blocks)[[1L]]
  cat(sprintf("NCGPE %s sdgp %.3f ours %.3f tol %.2e | brms's line const %.3e grad %.2e | rewritten const %.3e grad %.2e\n",
              key, exp(fit$estimates$theta[bk$theta_idx][1L]), ours, 1e-6 * abs(ours), a[1], a[2], b[1], b[2]))
})

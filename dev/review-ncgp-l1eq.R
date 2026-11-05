# Item 6: the five survey pairs on brms's own unedited program, and the same pairs with
# y multiplied so that sdgp is far from 1. Args: key scale response
# Rscript dev/review-ncgp-l1eq.R "student x gp" 100 y
args <- commandArgs(TRUE)
key <- args[1]; sc <- as.numeric(args[2]); resp <- args[3]
setwd("C:/Users/adf44/source/r/frmtmb-wt-ncgp")
source("dev/review-ncgp-lib.R")
rev_libs("lane")
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"))
suppressPackageStartupMessages(library(testthat))
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) sys.source(f, envir = E)
E$key <- key; E$sc <- sc; E$resp <- resp
local(envir = E, {
  spec <- eq_specs[[key]]()
  if (sc != 1) spec$data[[resp]] <- spec$data[[resp]] * sc
  fit <- suppressWarnings(eq_frm_fit(spec))
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
    up <- tryCatch(rstan::unconstrain_pars(sf, pars), error = function(err) err)
    if (inherits(up, "error")) return(c(NA, NA))
    lp <- tryCatch(rstan::log_prob(sf, up, adjust_transform = FALSE), error = function(err) NA)
    g <- tryCatch(rstan::grad_log_prob(sf, up, adjust_transform = FALSE), error = function(err) NA)
    idx <- brms_inner_index(sf, pars)
    c(lp - ours, max(abs(g[idx])))
  }
  a <- tryCatch(read(code), error = function(err) { cat("brms line error:", conditionMessage(err), "\n"); c(NA, NA) })
  b <- tryCatch(read(gsub("cov[n, n] += 1e-12;", "cov[n, n] += 1e-12 * square(sdgp);", code, fixed = TRUE)),
                error = function(err) { cat("rewrite error:", conditionMessage(err), "\n"); c(NA, NA) })
  bk <- Filter(function(k) identical(k$covstruct, "gp"), fit$frame$re_blocks)[[1L]]
  cat(sprintf("REV1 %s scale %g sdgp %.3f conv %d ours %.3f tol %.2e | brms's line const %.3e over %.1f grad %.2e | rewritten const %.3e grad %.2e\n",
              key, sc, exp(fit$estimates$theta[bk$theta_idx][1L]), fit$opt$convergence, ours, 1e-6 * abs(ours), a[1], abs(a[1]) / (1e-6 * abs(ours)), a[2], b[1], b[2]))
})

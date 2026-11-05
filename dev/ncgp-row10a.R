# Row 10a's gradient read on the BLAS it runs under: OpenBLAS 0.3.32 read 1.29e-3
# against the helper's 1e-3 (dev/ncgp-out/ob32/test-brms-likelihood.txt). What
# the number is on each, and whether the fit moved or only the rounding did.
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressPackageStartupMessages({library(testthat); library(frmtmb)})
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
cat("BLAS:", extSoftVersion()[["BLAS"]], "\n")
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("helper-warnings.R", "helper-brms.R")) sys.source(file.path("tests/testthat", f), envir = E)
local(envir = E, {
  set.seed(23)
  d10 <- data.frame(x = sort(runif(80, 0, 10)))
  d10$y <- sin(d10$x) + rnorm(80, 0, 0.3)
  fit <- frm(bf(y ~ gp(x)) + gaussian(), data = d10)
  bform <- brms::bf(y ~ gp(x))
  prior <- brms_flat_prior(bform, data = d10, family = gaussian())
  code <- brms::make_stancode(bform, data = d10, family = gaussian(), prior = prior)
  sdat <- brms_standata(bform, data = d10, family = gaussian(), prior = prior)
  rel <- gsub("cov[n, n] += 1e-12;", "cov[n, n] += 1e-12 * square(sdgp);", code, fixed = TRUE)
  sf <- suppressMessages(rstan::sampling(brms_stan_model(rel), data = sdat, chains = 0))
  pars <- stan_pars_from_fit(fit, sdat, rel)
  up <- rstan::unconstrain_pars(sf, pars)
  g <- rstan::grad_log_prob(sf, up, adjust_transform = FALSE)
  idx <- brms_inner_index(sf, pars)
  e <- fit$obj$env
  ours <- -as.numeric(e$f(e$last.par.best))
  lp <- rstan::log_prob(sf, up, adjust_transform = FALSE, gradient = FALSE)
  r <- e$random
  own <- -as.numeric(e$f(e$last.par.best, order = 1))[r]
  cat(sprintf("R10 ll %.10f conv %d iter %d sdgp %.8f rho %.8f sigma %.8f\n",
              as.numeric(logLik(fit)), fit$opt$convergence, fit$opt$iterations,
              exp(fit$estimates$theta[1]), exp(fit$estimates$theta[2]),
              as.numeric(sigma(fit))))
  cat(sprintf("R10 const %.3e stan inner grad %.3e, largest at z %d of %d | tape's inner grad %.3e\n",
              lp - ours, max(abs(g[idx])), which.max(abs(g[idx])), length(idx), max(abs(own))))
  cat("R10 top 5 of |stan grad|:", formatC(sort(abs(g[idx]), decreasing = TRUE)[1:5], digits = 2, format = "e"), "\n")
  cat("R10 |z| at those:", formatC(abs(up[idx])[order(abs(g[idx]), decreasing = TRUE)[1:5]], digits = 2, format = "e"), "\n")
})

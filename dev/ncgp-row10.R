# Row 10a's numbers: the design of test-brms-likelihood.R, brms's own line against
# the rewritten one. Seed 23.
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressPackageStartupMessages({library(testthat); library(frmtmb)})
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("helper-warnings.R", "helper-brms.R")) sys.source(file.path("tests/testthat", f), envir = E)
local(envir = E, {
  set.seed(23)
  n <- 80
  d <- data.frame(x = sort(runif(n, 0, 10)))
  d$y <- sin(d$x) + rnorm(n, 0, 0.3)
  fit <- frm(bf(y ~ gp(x)) + gaussian(), data = d)
  th <- fit$estimates$theta
  cat("sdgp", exp(th[1]), "rho", exp(th[2]), "\n")
  bform <- brms::bf(y ~ gp(x))
  prior <- brms_flat_prior(bform, data = d, family = gaussian())
  code <- brms::make_stancode(bform, data = d, family = gaussian(), prior = prior)
  sdat <- brms_standata(bform, data = d, family = gaussian(), prior = prior)
  pars <- stan_pars_from_fit(fit, sdat, code)
  e <- fit$obj$env
  ours <- -as.numeric(e$f(e$last.par.best)) + attr(pars, "logJ")
  for (cd in list(code, gsub("cov[n, n] += 1e-12;", "cov[n, n] += 1e-12 * square(sdgp);", code, fixed = TRUE))) {
    sf <- suppressMessages(rstan::sampling(brms_stan_model(cd), data = sdat, chains = 0))
    up <- rstan::unconstrain_pars(sf, pars)
    lp <- rstan::log_prob(sf, up, adjust_transform = FALSE)
    g <- rstan::grad_log_prob(sf, up, adjust_transform = FALSE)
    idx <- brms_inner_index(sf, pars)
    cat(sprintf("NCGPT ours %.5f tol %.2e const %.4e grad %.3e\n", ours, 1e-6 * abs(ours), lp - ours, max(abs(g[idx]))))
  }
})

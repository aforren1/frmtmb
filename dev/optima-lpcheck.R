# Lane optima: test-brms-likelihood.R's check C row 3 variant by hand,
# with brms_lp_check()'s report on, on the arm given.
#   Rscript dev/optima-lpcheck.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache"))
suppressMessages({library(testthat); library(frmtmb)})
env <- new.env(parent = asNamespace("frmtmb"))
testthat::source_test_helpers("tests/testthat", env = env)
options(frmtmb.brms_lp_report = TRUE)
evalq({
  set.seed(3)
  dm <- data.frame(inc = sample(0:3, 300, TRUE), z = rnorm(300),
                   g = factor(rep(1:20, 15)))
  dm$y <- 1 + c(0, 1, 1.6, 2)[dm$inc + 1] + 0.3 * dm$z + rnorm(300)
  fit <- frm(bf(y ~ mo(inc):z + (1 | g)) + gaussian(), data = dm)
  print(fit$opt$par)
  print(frmtmb:::mo_simplex(fit$estimates$zeta1))
  r <- tryCatch(brms_lp_check(brms::bf(y ~ mo(inc):z + (1 | g)), gaussian(),
                              dm, fit, joint = TRUE, const = lgamma(3)),
                error = function(e) conditionMessage(e))
  if (is.list(r)) print(r$pars$simo_1) else print(r)
}, env)
evalq({
  e <- fit$obj$env
  a <- -e$f(e$last.par.best)
  fit$obj$fn(fit$opt$par)
  b <- -e$f(e$last.par)
  cat("joint at last.par.best", format(a, digits = 12),
      "; joint after re-solving at opt$par", format(b, digits = 12), "\n")
  r <- e$random
  cat("max |random effects difference|",
      max(abs(e$last.par.best[r] - e$last.par[r])), "\n")
}, env)
evalq({
  bform <- brms::bf(y ~ mo(inc):z + (1 | g))
  prior <- brms_flat_prior(bform, data = dm, family = gaussian())
  code <- brms::make_stancode(bform, data = dm, family = gaussian(),
                              prior = prior)
  sdat <- brms_standata(bform, data = dm, family = gaussian(), prior = prior)
  rtab <- brms_ranef_table(bform, dm, gaussian(), prior)
  sf <- suppressMessages(rstan::sampling(brms_stan_model(code), data = sdat,
                                         chains = 0))
  pars <- stan_pars_from_fit(fit, sdat, code, rtab)
  lp_at <- function(p) rstan::log_prob(sf, rstan::unconstrain_pars(sf, p),
                                       adjust_transform = FALSE,
                                       gradient = FALSE)
  cat("lp at the lane's pars:", format(lp_at(pars), digits = 12), "\n")
  for (w2 in c(1e-12, 1e-10, 1e-8, 1e-6)) {
    q <- pars
    s <- q$simo_1
    s[2] <- w2
    s <- s / sum(s)
    q$simo_1 <- s
    cat("simo_1[2] =", w2, ": lp", format(lp_at(q), digits = 12), "\n")
  }
  cat(grep("simo|Dirichlet|dirichlet", strsplit(code, "\n")[[1]],
           value = TRUE), sep = "\n")
}, env)

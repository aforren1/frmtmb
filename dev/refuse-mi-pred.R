# D8, the predictor side: mi(w) in a formula under each family brms takes
# it for, the variable imputed by a gaussian model. frm() fit and, where the
# translator has the family, L1 against brms's log_prob. Seed 170.
.libPaths(c("C:/Users/adf44/source/r/wt-refuse-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           NOT_CRAN = "true")
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
fams <- list(
  list("poisson", "poisson()", "poisson()", ""),
  list("bernoulli", "bernoulli()", "bernoulli()", ""),
  list("binomial", "binomial()", "binomial()", " | trials(10)"),
  list("negbinomial", "negbinomial()", "negbinomial()", ""),
  list("cumulative", "cumulative()", "cumulative()", ""),
  list("sratio", "sratio()", "sratio()", ""),
  list("categorical", "categorical()", "categorical()", ""),
  list("Gamma", "stats::Gamma(link = 'log')", "brmsfamily('gamma', 'log')", ""),
  list("beta", "Beta()", "brmsfamily('beta')", ""),
  list("student", "student()", "student()", ""),
  list("hurdle_gamma", "hurdle_gamma()", "hurdle_gamma()", ""),
  list("zero_inflated_poisson", "zero_inflated_poisson()",
       "zero_inflated_poisson()", ""),
  list("von_mises", "von_mises()", "von_mises()", ""),
  list("asym_laplace", "asym_laplace()", "asym_laplace()", ""))
set.seed(170)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$wf <- d$z + rnorm(n, 0, 0.5)
d$w <- d$wf
d$w[sample(n, 15)] <- NA
eta <- 0.3 + 0.4 * d$x + 0.4 * d$wf
gen <- function(nm) switch(nm,
  poisson = rpois(n, exp(eta)), bernoulli = rbinom(n, 1, plogis(eta)),
  binomial = rbinom(n, 10, plogis(eta)), negbinomial = rnbinom(n, 3, mu = exp(eta)),
  cumulative = , sratio = , categorical =
    findInterval(eta + rlogis(n), c(-1, 0.5, 2)) + 1L,
  Gamma = rgamma(n, 3, 3 / exp(eta)),
  beta = rbeta(n, plogis(eta) * 8, (1 - plogis(eta)) * 8),
  student = eta + rt(n, 4),
  hurdle_gamma = ifelse(runif(n) < 0.2, 0, rgamma(n, 3, 3 / exp(eta))),
  zero_inflated_poisson = ifelse(runif(n) < 0.2, 0, rpois(n, exp(eta))),
  von_mises = ((2 * atan(eta) + rnorm(n, 0, 0.5) + pi) %% (2 * pi)) - pi,
  asym_laplace = eta + rexp(n) - rexp(n))
for (fm in fams) {
  nm <- fm[[1L]]
  dd <- d
  dd$v <- gen(nm)
  rhs <- "mi(w) + x"
  ff <- tryCatch({
    frm_fam <- eval(parse(text = fm[[2L]]), asNamespace("frmtmb"))
    form <- as.formula(paste0("v", fm[[4L]], " ~ ", rhs))
    suppressWarnings(frm(mvbf(bf(form, family = frm_fam),
                              bf(w | mi() ~ z), rescor = FALSE), data = dd))
  }, error = function(e) e)
  if (inherits(ff, "error")) {
    cat(sprintf("%-22s frm: REFUSED %s\n", nm,
                substr(conditionMessage(ff), 1, 120)))
    next
  }
  l1 <- tryCatch({
    bfam <- eval(parse(text = fm[[3L]]), asNamespace("brms"))
    fb <- brms::bf(as.formula(paste0("v", fm[[4L]], " ~ ", rhs)), family = bfam) +
      brms::bf(w | mi() ~ z) + brms::set_rescor(FALSE)
    pr <- E$brms_flat_prior(fb, data = dd, family = NULL)
    code <- brms::make_stancode(fb, data = dd, prior = pr)
    sdat <- E$brms_standata(fb, data = dd, family = NULL, prior = pr)
    mod <- E$brms_stan_model(code)
    sf <- suppressMessages(rstan::sampling(mod, data = sdat, chains = 0))
    pars <- E$stan_pars_from_fit(ff, sdat, code, NULL)
    up <- rstan::unconstrain_pars(sf, pars)
    lp <- rstan::log_prob(sf, up, adjust_transform = FALSE, gradient = FALSE)
    ours <- -ff$obj$env$f(ff$obj$env$last.par.best) + attr(pars, "logJ")
    g <- E$brms_lp_grad(sf, up, ff, pars)
    idx <- E$brms_inner_index(sf, pars)
    sprintf("const %.3g inner grad %.3g", lp - ours, max(abs(g[idx])))
  }, error = function(e) paste("L1 error:", substr(conditionMessage(e), 1, 90)))
  cat(sprintf("%-22s frm: logLik %.4f  %s\n", nm, as.numeric(logLik(ff)), l1))
}

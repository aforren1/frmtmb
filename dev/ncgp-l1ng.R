# The non-gaussian designs of dev/ncgp-nongauss.R on brms's own program at frm()'s
# estimate, and rewritten for the attribution, because brms's chains do not mix at a
# large log mean. The joint and the inner gradient say the two tapes have one mode;
# the Laplace term on top is checked by hand in test-gp-noncentered.R.
#   Rscript dev/ncgp-l1ng.R <case>
case <- commandArgs(TRUE)[1]
LIB <- "C:/Users/adf44/source/r/wt-ncgp-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressPackageStartupMessages({library(testthat); library(frmtmb)})
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("helper-warnings.R", "helper-brms.R")) sys.source(file.path("tests/testthat", f), envir = E)
n <- 40L
set.seed(5601)
x <- sort(round(stats::runif(n, 0, 6), 3))
f <- sin(2 * x)
spec <- switch(
  case,
  pois4 = list(y = stats::rpois(n, exp(4 + f)), fam = "poisson"),
  pois12 = list(y = stats::rpois(n, exp(12 + f)), fam = "poisson"),
  pois16 = list(y = stats::rpois(n, exp(16 + f)), fam = "poisson"),
  binom1e5 = list(y = stats::rbinom(n, 1e5, plogis(f)), fam = "binomial", N = 1e5),
  binom1e7 = list(y = stats::rbinom(n, 1e7, plogis(f)), fam = "binomial", N = 1e7),
  bern = list(y = stats::rbinom(200L, 1L, plogis(rep(f, 5L))), fam = "bernoulli"))
if (case == "bern") x <- rep(x, 5L)
d <- data.frame(x = x, y = spec$y)
if (!is.null(spec$N)) d$N <- spec$N
form <- function(ns) {
  bfn <- get("bf", ns)
  if (!is.null(spec$N)) bfn(y | trials(N) ~ gp(x)) else bfn(y ~ gp(x))
}
famf <- function() switch(spec$fam, poisson = poisson(), binomial = binomial(), bernoulli = bernoulli())
fit <- suppressWarnings(frm(form(asNamespace("frmtmb")), data = d, family = famf()))
bform <- form(asNamespace("brms"))
prior <- E$brms_flat_prior(bform, data = d, family = famf())
code <- brms::make_stancode(bform, data = d, family = famf(), prior = prior)
sdat <- E$brms_standata(bform, data = d, family = famf(), prior = prior)
pars <- E$stan_pars_from_fit(fit, sdat, code)
e <- fit$obj$env
ours <- -as.numeric(e$f(e$last.par.best)) + attr(pars, "logJ")
read <- function(cd) {
  sf <- suppressMessages(rstan::sampling(E$brms_stan_model(cd), data = sdat, chains = 0))
  up <- rstan::unconstrain_pars(sf, pars)
  lp <- rstan::log_prob(sf, up, adjust_transform = FALSE)
  g <- rstan::grad_log_prob(sf, up, adjust_transform = FALSE)
  idx <- E$brms_inner_index(sf, pars)
  c(lp - ours, max(abs(g[idx])))
}
a <- read(code)
b <- read(gsub("cov[n, n] += 1e-12;", "cov[n, n] += 1e-12 * square(sdgp);", code, fixed = TRUE))
cat(sprintf("NCGPLN %s sdgp %.3f ours %.3f tol %.2e | brms's line const %.3e grad %.2e | rewritten const %.3e grad %.2e\n",
            case, exp(fit$estimates$theta[1]), ours, 1e-6 * abs(ours), a[1], a[2], b[1], b[2]))

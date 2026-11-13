# Review of lane nlstop, item 4850: is the lane's point the posterior mode
# that brms 2.23.0's program reaches? rstan::optimizing() (LBFGS, no
# Jacobian: the mode in sigma, as frmtmb's penalized objective) on brms's
# program with the same priors and a flat sigma prior (frmtmb sets none),
# from 6 random inits, from the lane's point and from 0.81.0's point
# ("off": prior_conflict_escape() as the identity).
# Usage: Rscript dev/reviews/nlstop-review-brms.R
source("dev/nlstop-common.R")
suppressPackageStartupMessages({library(brms); library(rstan)})
ns <- asNamespace("frmtmb")
orig <- ns$prior_conflict_escape
set_arm <- function(arm) {
  utils::assignInNamespace("prior_conflict_escape",
    if (arm == "on") orig else function(obj, opt, ...) opt, "frmtmb")
}
source("tests/testthat/helper-brms.R")
cat("brms", as.character(packageVersion("brms")), "rstan",
    as.character(packageVersion("rstan")), "\n")
mk <- function(seed, psd, s = 1e6) {
  set.seed(seed); n <- 200
  d <- data.frame(x = rnorm(n))
  d$y <- 0.5 * d$x + rnorm(n)
  d$yl <- s * (1 + 0.05 * d$x) + rnorm(n, 0, s * 1e-3)
  d
}
cases <- list(
  v6201 = list(seed = 6201, psd = 1000, cls = "b"),
  v6202 = list(seed = 6202, psd = 1000, cls = "b"),
  v6201p300 = list(seed = 6201, psd = 300, cls = "b"),
  v6203p3000 = list(seed = 6203, psd = 3000, cls = "b"),
  int1e4 = list(seed = 6201, psd = 1e4, cls = "Intercept"),
  int3e3 = list(seed = 6201, psd = 3e3, cls = "Intercept"),
  int2e4 = list(seed = 6201, psd = 2e4, cls = "Intercept"))
pt_of <- function(r, d) {
  fe <- fixef(r)[, "Estimate"]
  sg <- exp(unname(r$estimates$betad[1]))
  list(b_Intercept = unname(fe[1]), b = unname(fe[2]), sigma = sg,
       Intercept = unname(fe[1] + fe[2] * mean(d$x)))
}
for (cn in names(cases)) {
  cs <- cases[[cn]]
  d <- mk(cs$seed, cs$psd)
  pr <- if (cs$cls == "b") {
    c(set_prior(sprintf("normal(0, %g)", cs$psd), class = "b"),
      set_prior("normal(1e6, 1e6)", class = "Intercept"))
  } else {
    set_prior(sprintf("normal(0, %g)", cs$psd), class = "Intercept")
  }
  set_arm("on")
  r_on <- tt(frm(yl ~ x, data = d, prior = pr))
  set_arm("off")
  r_off <- tt(frm(yl ~ x, data = d, prior = pr))
  set_arm("on")
  p_on <- pt_of(r_on, d); p_off <- pt_of(r_off, d)
  cat(sprintf("CASE %s frm on: llp %.6f b_Int %.6g b %.6g sigma %.6g msg %s warn %s\n",
              cn, llp(r_on), p_on$b_Intercept, p_on$b, p_on$sigma,
              any(attr(r_on, "mc") == "frmtmb_prior_multistart"), wsum(r_on)))
  cat(sprintf("CASE %s frm off: llp %.6f b_Int %.6g b %.6g sigma %.6g\n",
              cn, llp(r_off), p_off$b_Intercept, p_off$b, p_off$sigma))
  bpr <- c(pr, set_prior("", class = "sigma"))
  code <- make_stancode(yl ~ x, data = d, prior = bpr)
  sdat <- make_standata(yl ~ x, data = d, prior = bpr)
  mod <- brms_stan_model(code)
  inits <- list()
  for (k in 1:6) inits[[paste0("rand", k)]] <- "random"
  inits$at_on <- list(Intercept = p_on$Intercept, b = array(p_on$b),
                      sigma = p_on$sigma)
  inits$at_off <- list(Intercept = p_off$Intercept, b = array(p_off$b),
                       sigma = p_off$sigma)
  for (nm in names(inits)) {
    seed <- match(nm, names(inits))
    o <- tryCatch(suppressWarnings(optimizing(mod, data = sdat,
      init = inits[[nm]], seed = 100 + seed, algorithm = "LBFGS",
      iter = 20000, as_vector = FALSE, verbose = FALSE, hessian = FALSE)),
      error = function(e) NULL)
    if (is.null(o)) { cat(sprintf("  stan %s: error\n", nm)); next }
    bi <- o$par$Intercept - o$par$b * mean(d$x)
    cat(sprintf("  stan %s: lp %.6f rc %d b_Int %.6g b %.6g sigma %.6g\n",
                nm, o$value, o$return_code, bi, o$par$b[1], o$par$sigma))
  }
  flush.console()
}
cat("DONE\n")

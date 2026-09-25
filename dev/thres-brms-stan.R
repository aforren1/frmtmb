# Lane thres: brms's COMPILED log density against frmtmb, grouped
# thresholds, at frmtmb's estimate.
#
#   Rscript dev/thres-brms-stan.R > dev/thres-brms-stan-log.txt 2>&1
#
# For each family, brms 2.23.0 writes the Stan program for
# `y | thres(gr = g) ~ x` with every prior flat; rstan compiles it
# (about a minute each) and evaluates log_prob and its gradient at
# frmtmb's maximum likelihood estimate, with adjust_transform = FALSE so
# that no Jacobian enters. The log density must equal frmtmb's logLik
# (brms keeps every constant, `lpmf` not `lupmf`), and its gradient
# must vanish there, which says frmtmb's estimate is brms's mode too.
# A last run puts brms's normal(0, 2) on class Intercept, with the
# default rows per group, and compares against frmtmb's penalized
# objective; frmtmb's density on an ORDERED vector carries the
# log-Jacobian of (first, log increments) -> thresholds, sum of the log
# increments, which log_prob(adjust_transform = FALSE) leaves out, so
# that sum is added back on frmtmb's side.
# Data seed 11, n = 240, three levels with 4, 3 and 5 categories.
.libPaths(c("/opt/rlib/stan", "/opt/rlib/lane-thres", "/opt/rlib/base",
            "/opt/rlib/deps", "/opt/r/lib/R/library"))
suppressMessages({
  library(frmtmb)
  library(brms)
  library(rstan)
})
cat("frmtmb from", find.package("frmtmb"), "\n")

set.seed(11)
n <- 240
d <- data.frame(x = rnorm(n), g = sample(c("a", "b", "c"), n, TRUE))
u <- rlogis(n, 0.7 * d$x)
tau <- list(a = c(-1, 0.2, 1.3), b = c(-0.5, 0.6),
            c = c(-1.5, -0.4, 0.5, 1.6))
d$y <- vapply(seq_len(n), function(i) 1L + sum(u[i] > tau[[d$g[i]]]), 1L)

stan_at <- function(bform, family, prior, fit) {
  code <- make_stancode(bform, data = d, family = family, prior = prior)
  sdat <- make_standata(bform, data = d, family = family, prior = prior)
  mod <- stan_model(model_code = code)
  sf <- suppressMessages(sampling(mod, data = sdat, chains = 0))
  fe <- fixef(fit)[, "Estimate"]
  th <- fit$spec$responses$y$family[["thres"]]
  pars <- list(b = array(fe[["x"]], 1L))
  for (k in seq_along(th$groups)) {
    nm <- paste0("Intercept[", th$groups[k], ",", seq_len(th$nthres[k]),
                 "]")
    pars[[paste0("Intercept_", k)]] <- as.array(unname(fe[nm]))
  }
  up <- unconstrain_pars(sf, pars)
  list(lp = log_prob(sf, up, adjust_transform = FALSE),
       grad = grad_log_prob(sf, up, adjust_transform = FALSE))
}

cat("\n== flat priors: brms log_prob vs frmtmb logLik at the MLE ==\n")
for (fam in c("cumulative", "sratio", "cratio", "acat")) {
  ff <- get(fam, envir = asNamespace("frmtmb"))
  fit <- frm(y | thres(gr = g) ~ x, data = d, family = ff())
  bf0 <- bf(y | thres(gr = g) ~ x)
  bfam <- get(fam, envir = asNamespace("brms"))()
  pr <- get_prior(bf0, data = d, family = bfam)
  pr$prior <- ""
  r <- stan_at(bf0, bfam, pr, fit)
  ll <- as.numeric(logLik(fit))
  cat(sprintf("%-10s brms %.10f frmtmb %.10f rel %.2e max|grad| %.2e\n",
              fam, r$lp, ll, abs(r$lp - ll) / abs(ll), max(abs(r$grad))))
}

cat("\n== normal(0, 2) on class Intercept, cumulative ==\n")
pf <- set_prior("normal(0, 2)", class = "Intercept")
fit <- frm(y | thres(gr = g) ~ x, data = d, family = cumulative(),
           prior = pf)
bf0 <- bf(y | thres(gr = g) ~ x)
pr <- get_prior(bf0, data = d, family = brms::cumulative())
pr$prior <- ""
pr$prior[pr$class == "Intercept" & pr$coef == "" & pr$group == ""] <-
  "normal(0, 2)"
r <- stan_at(bf0, brms::cumulative(), pr, fit)
raw <- fit$estimates$tau_raw
th <- fit$spec$responses$y$family[["thres"]]
lay <- frmtmb:::thres_layout(th$nthres)
logjac <- sum(unlist(lapply(seq_len(lay$G), function(g) {
  if (lay$nthres[g] > 1L) raw[(lay$start[g] + 1L):lay$end[g]] else 0
})))
ours <- -fit$opt$objective - logjac
cat(sprintf("brms %.10f frmtmb %.10f rel %.2e max|grad| %.2e\n", r$lp,
            ours, abs(r$lp - ours) / abs(ours), max(abs(r$grad))))

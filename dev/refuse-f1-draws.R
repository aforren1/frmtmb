# M5 and M3: frm_sample()'s missing values against brms's Ymi, by name and by
# KS, and a mixture on a bounded support. Usage: Rscript
# dev/refuse-f1-draws.R <Gamma|von_mises|beta|mix>. Seed 9700.
which_case <- commandArgs(trailingOnly = TRUE)[1L]
.libPaths(c("C:/Users/adf44/source/r/wt-refuse-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb); library(frmtmb.sample); library(testthat)
})
E <- new.env(parent = asNamespace("frmtmb"))
sys.source("tests/testthat/helper-brms.R", envir = E)
set.seed(9700)
n <- 120
d <- data.frame(x = rnorm(n))
eta <- 0.4 + 0.5 * d$x
d$y <- switch(which_case,
  Gamma = rgamma(n, 2.5, 2.5 / exp(eta)),
  beta = rbeta(n, plogis(eta) * 4, (1 - plogis(eta)) * 4),
  von_mises = ((2 * atan(0.5 * eta) + rnorm(n, 0, 0.5) + pi) %% (2 * pi)) - pi,
  mix = {
    k <- rbinom(n, 1, 0.5)
    exp(ifelse(k == 1, 0.6, -0.6) + 0.3 * d$x + rnorm(n, 0, 0.6))
  })
d$y[sample(n, 15)] <- NA
frm_f <- switch(which_case,
  Gamma = stats::Gamma(link = "log"), beta = Beta(), von_mises = von_mises(),
  mix = mixture(lognormal(), lognormal()))
brm_f <- switch(which_case,
  Gamma = brms::brmsfamily("gamma", "log"), beta = brms::brmsfamily("beta"),
  von_mises = brms::brmsfamily("von_mises"),
  mix = brms::mixture(brms::lognormal(), brms::lognormal(), order = "none"))
fb <- if (which_case == "mix") {
  brms::bf(y | mi() ~ 1, mu1 ~ x, mu2 ~ x)
} else {
  brms::bf(y | mi() ~ x)
}
ff <- if (which_case == "mix") {
  bf(y | mi() ~ 1, mu1 ~ x, mu2 ~ x)
} else {
  bf(y | mi() ~ x)
}
# two overlapping components, so that a missing value has one mode to mix
# over, each held at its side by a prior, so that
# neither program's chains can swap them
mix_pr <- c("set_prior('normal(0.6, 0.1)', class = 'Intercept', dpar = 'mu1')",
            "set_prior('normal(-0.6, 0.1)', class = 'Intercept', dpar = 'mu2')",
            "set_prior('normal(0, 1)', class = 'b', dpar = 'mu1')",
            "set_prior('normal(0, 1)', class = 'b', dpar = 'mu2')")
ours <- suppressWarnings(suppressMessages(frm_sample(
  ff, family = frm_f, data = d, chains = 4, iter = 2000, seed = 11,
  refresh = 0, cores = 1,
  prior = if (which_case == "mix") {
    eval(parse(text = paste(mix_pr, collapse = " + ")))
  })))
pr <- brms::default_prior(fb, family = brm_f, data = d)
if (which_case == "mix") {
  pr <- eval(parse(text = paste(mix_pr, collapse = " + ")),
             envir = asNamespace("brms"))
}
code <- brms::stancode(fb, family = brm_f, data = d, prior = pr)
sdat <- brms::standata(fb, family = brm_f, data = d, prior = pr)
mod <- E$brms_stan_model(code)
ini <- if (which_case == "mix") {
  list(list(Ymi = array(rep(1, 15), 15L)))[rep(1L, 4L)]
} else {
  "random"
}
sf <- rstan::sampling(mod, data = sdat, chains = 4, iter = 2000, seed = 12,
                      refresh = 0, init = ini)
ex <- rstan::extract(sf, permuted = TRUE)
ym <- as.matrix(ex$Ymi)
dm <- as.matrix(as_draws_matrix(ours))
cn <- grep("^Ymi", colnames(dm), value = TRUE)
cat("columns:", head(cn, 3), "... n =", length(cn), "\n")
ks <- vapply(seq_along(cn), function(k) {
  suppressWarnings(stats::ks.test(dm[, cn[k]], ym[, k])$statistic)
}, 0)
cat(sprintf("%s KS of Ymi: median %.3f max %.3f; our range [%.2f, %.2f], brms [%.2f, %.2f]\n",
            which_case, median(ks), max(ks), min(dm[, cn]), max(dm[, cn]),
            min(ym), max(ym)))
for (nm in c("shape", "phi", "kappa", "b")) {
  if (!is.null(ex[[nm]])) {
    cat(nm, "brms mean", format(mean(ex[[nm]]), digits = 4), "\n")
  }
}
cn2 <- intersect(c("b_x", "shape", "phi", "kappa", "b_mu1_x", "b_mu2_x",
                   "sigma1", "sigma2"), colnames(dm))
print(round(colMeans(dm[, cn2, drop = FALSE]), 3))
rh <- rstan::summary(sf)$summary
cat("brms max Rhat", format(max(rh[, "Rhat"], na.rm = TRUE), digits = 4), "\n")
ep <- posterior_epred(ours)
pp <- posterior_predict(ours)
cat("readers: epred", dim(ep), "finite", all(is.finite(ep)), "; predict",
    dim(pp), "finite", all(is.finite(pp)), "\n")
ll <- tryCatch(dim(log_lik(ours)), error = function(e) conditionMessage(e))
cat("log_lik:", substr(paste(ll, collapse = " "), 1, 100), "\n")
s <- posterior_summary(ours)
cat("summary has Ymi:", any(grepl("^Ymi", rownames(s))), "\n")

# Item 2: dpar readers against brms 2.23.0 at the same parameter values, per
# family class. Seed 7101. Writes dev/review-predread-out/dpar-<tag>.tsv
#   REV_LIB=before Rscript dev/review-predread-dpar.R before
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
tag <- commandArgs(TRUE)[1]
only <- commandArgs(TRUE)[2]
set.seed(7101)
n <- 240
d <- data.frame(x = rnorm(n), z = rnorm(n))
eta <- 0.3 + 0.6 * d$x - 0.4 * d$z
d$y <- eta + rnorm(n, 0, exp(-0.2 + 0.3 * d$z))
d$pos <- exp(eta + rnorm(n, 0, 0.4))
d$cnt <- rpois(n, exp(eta))
d$bin <- rbinom(n, 1, plogis(eta))
d$frac <- rbeta(n, plogis(eta) * 8, (1 - plogis(eta)) * 8)
d$zcnt <- ifelse(runif(n) < 0.3, 0L, d$cnt)
d$hcnt <- ifelse(d$cnt == 0, 1L, d$cnt)
d$hcnt[runif(n) < 0.3] <- 0L
d$hpos <- ifelse(runif(n) < plogis(-1 + d$z), 0, d$pos)
d$ord <- as.integer(cut(eta + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf),
                        labels = FALSE))
d$cat <- as.integer(cut(eta + rlogis(n), c(-Inf, 0, 1, Inf), labels = FALSE))
d$tr <- 10L
d$succ <- rbinom(n, 10, plogis(eta))
comp <- rbinom(n, 1, 0.4)
d$mixy <- ifelse(comp == 1, 2 + 0.8 * d$x, -1 - 0.5 * d$x) + rnorm(n, 0, 0.6)
d$mixc <- rpois(n, ifelse(comp == 1, exp(1.5 + 0.4 * d$x), exp(0.2 - 0.3 * d$x)))
d$mixp <- exp(d$mixy / 2)
d$hord <- d$ord
d$hord[runif(n) < 0.2] <- 0L
d$ang <- atan2(sin(0.5 * d$x + rnorm(n, 0, 0.5)), cos(0.5 * d$x + rnorm(n, 0, 0.5)))
d$wei <- rweibull(n, 1.5, exp(eta))
nd <- data.frame(x = c(-1, 0.3, 1.2), z = c(0.5, -0.4, 1), tr = 10L)
specs <- list(
  identity = list(f = "bf(y ~ x, sigma ~ z)", fam = "gaussian()",
                  dp = c("mu", "sigma")),
  student = list(f = "bf(y ~ x, sigma ~ z, nu ~ 1)", fam = "student()",
                 dp = c("mu", "sigma", "nu")),
  log = list(f = "bf(cnt ~ x)", fam = "poisson()", dp = "mu"),
  logit = list(f = "bf(bin ~ x)", fam = "bernoulli()", dp = "mu"),
  binom = list(f = "bf(succ | trials(tr) ~ x)", fam = "binomial()", dp = "mu"),
  probit = list(f = "bf(bin ~ x)", fam = "bernoulli('probit')", dp = "mu"),
  inverse = list(f = "bf(pos ~ x, shape ~ z)", fam = "Gamma('inverse')",
                 dp = c("mu", "shape")),
  gamma_log = list(f = "bf(pos ~ x, shape ~ z)", fam = "Gamma('log')",
                   dp = c("mu", "shape")),
  invgauss = list(f = "bf(pos ~ x)", fam = "inverse.gaussian()",
                  dp = c("mu", "shape")),
  lognormal = list(f = "bf(pos ~ x, sigma ~ z)", fam = "lognormal()",
                   dp = c("mu", "sigma")),
  beta = list(f = "bf(frac ~ x, phi ~ z)", fam = "Beta()",
              dp = c("mu", "phi")),
  negbin = list(f = "bf(cnt ~ x, shape ~ z)", fam = "negbinomial()",
                dp = c("mu", "shape")),
  weibull = list(f = "bf(wei ~ x, shape ~ z)", fam = "weibull()",
                 dp = c("mu", "shape")),
  vonmises = list(f = "bf(ang ~ x)", fam = "von_mises()", dp = c("mu", "kappa")),
  cumulative = list(f = "bf(ord ~ x)", fam = "cumulative()", dp = c("mu", "disc")),
  cum_probit = list(f = "bf(ord ~ x)", fam = "cumulative('probit')",
                    dp = c("mu", "disc")),
  sratio_cs = list(f = "bf(ord ~ x + cs(z))", fam = "sratio('cloglog')",
                   dp = c("mu", "disc")),
  acat = list(f = "bf(ord ~ x, disc ~ z)", fam = "acat('logit')",
              dp = c("mu", "disc")),
  cratio = list(f = "bf(ord ~ x)", fam = "cratio('cauchit')", dp = "mu"),
  hurdle_pois = list(f = "bf(hcnt ~ x, hu ~ z)", fam = "hurdle_poisson()",
                     dp = c("mu", "hu")),
  hurdle_ln = list(f = "bf(hpos ~ x, hu ~ z)", fam = "hurdle_lognormal()",
                   dp = c("mu", "hu", "sigma")),
  hurdle_gam = list(f = "bf(hpos ~ x, hu ~ z)", fam = "hurdle_gamma()",
                    dp = c("mu", "hu", "shape")),
  zip = list(f = "bf(zcnt ~ x, zi ~ z)", fam = "zero_inflated_poisson()",
             dp = c("mu", "zi")),
  zinb = list(f = "bf(zcnt ~ x, zi ~ z)", fam = "zero_inflated_negbinomial()",
              dp = c("mu", "zi", "shape")),
  hurdle_cum = list(f = "bf(hord ~ x, hu ~ z)", fam = "hurdle_cumulative()",
                    dp = c("mu", "hu")),
  mix_gauss = list(f = "bf(mixy ~ x)", fam = "mixture(gaussian(), gaussian())",
                   dp = c("mu1", "mu2", "sigma1", "theta1")),
  mix_pois = list(f = "bf(mixc ~ x)", fam = "mixture(poisson(), poisson())",
                  dp = c("mu1", "mu2", "theta1")),
  mix_mixed = list(f = "bf(mixp ~ x)", fam = "mixture(lognormal(), Gamma('log'))",
                   dp = c("mu1", "mu2", "sigma1", "shape2")),
  mix_ord = list(f = "bf(ord ~ x)",
                 fam = "mixture(cumulative(), cumulative('probit'))",
                 dp = c("mu1", "mu2")),
  categorical = list(f = "bf(cat ~ x)", fam = "categorical()",
                     dp = c("mu2", "mu3")))
sel <- names(specs)
if (!is.na(only)) sel <- strsplit(only, ",")[[1L]]
out <- list()
add <- function(...) out[[length(out) + 1L]] <<- data.frame(...)
num <- function(v) signif(as.numeric(v), 6)
for (nm in sel) {
  sp <- specs[[nm]]
  cat("== ", nm, "\n")
  pr <- tryCatch(rev_pair(sp$f, sp$fam, d), error = function(e) e)
  if (inherits(pr, "error")) {
    add(spec = nm, dpar = "", reader = "pair", rel = NA, a = NA, b = NA,
        note = substr(conditionMessage(pr), 1, 150))
    next
  }
  fit <- pr$fit
  fx <- pr$fx
  ds <- rev_draws_at(fit, 2L)
  # the draws' own estimates, not zeros: read at the estimate
  for (dp in sp$dp) {
    bm <- function(e) {
      v <- try_num(e)
      if (!is.null(attr(v, "err"))) return(v)
      if (length(dim(v)) == 2L) v[1, ] else v
    }
    pairs <- list(
      "fitted response" = list(
        function() fitted(fit, newdata = nd, dpar = dp)[, "Estimate"],
        function() fitted(fx, newdata = nd, dpar = dp, summary = FALSE)),
      "fitted linear" = list(
        function() fitted(fit, newdata = nd, dpar = dp, scale = "linear")[, "Estimate"],
        function() fitted(fx, newdata = nd, dpar = dp, scale = "linear", summary = FALSE)),
      "frm_linpred response" = list(
        function() frm_linpred(fit, newdata = nd, dpar = dp, type = "response"),
        function() posterior_linpred(fx, newdata = nd, dpar = dp, transform = TRUE)),
      "frm_linpred link" = list(
        function() frm_linpred(fit, newdata = nd, dpar = dp, type = "link"),
        function() posterior_linpred(fx, newdata = nd, dpar = dp)),
      "frm_curve transform" = list(
        function() frm_curve(fit, newdata = nd, dpar = dp, simultaneous = FALSE,
                             transform = TRUE)$.estimate,
        function() posterior_epred(fx, newdata = nd, dpar = dp)),
      "frm_curve link" = list(
        function() frm_curve(fit, newdata = nd, dpar = dp, simultaneous = FALSE)$.estimate,
        function() posterior_linpred(fx, newdata = nd, dpar = dp)),
      "CE dpar" = list(
        function() conditional_effects(fit, "x", dpar = dp, resolution = 3,
                                       conditions = data.frame(z = 0.1))[[1L]]$estimate__,
        function() brms::conditional_effects(fx, "x", dpar = dp, resolution = 3,
                                              conditions = data.frame(z = 0.1))[[1L]]$estimate__),
      "posterior_epred" = list(
        function() frmtmb.sample::posterior_epred(ds, newdata = nd, dpar = dp)[1, ],
        function() posterior_epred(fx, newdata = nd, dpar = dp)),
      "posterior_linpred" = list(
        function() frmtmb.sample::posterior_linpred(ds, newdata = nd, dpar = dp)[1, ],
        function() posterior_linpred(fx, newdata = nd, dpar = dp)),
      "posterior_linpred transform" = list(
        function() frmtmb.sample::posterior_linpred(ds, newdata = nd, dpar = dp, transform = TRUE)[1, ],
        function() posterior_linpred(fx, newdata = nd, dpar = dp, transform = TRUE)),
      "draws fitted linear" = list(
        function() frmtmb.sample:::fitted.frmtmb_draws(ds, newdata = nd, dpar = dp, scale = "linear", summary = FALSE)[1, ],
        function() fitted(fx, newdata = nd, dpar = dp, scale = "linear", summary = FALSE)))
    for (r in names(pairs)) {
      a <- try_num(pairs[[r]][[1L]]())
      b <- bm(pairs[[r]][[2L]]())
      add(spec = nm, dpar = dp, reader = r, rel = signif(rd(a, b), 3),
          a = paste(num(a), collapse = " "), b = paste(num(b), collapse = " "),
          note = paste0(attr(a, "err") %||% "", " | ", attr(b, "err") %||% ""))
    }
  }
}
res <- do.call(rbind, out)
dir.create("dev/review-predread-out", showWarnings = FALSE)
write.table(res, sprintf("dev/review-predread-out/dpar-%s.tsv", tag), sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("rows", nrow(res), "agree", sum(res$rel < 1e-8, na.rm = TRUE), "differ",
    sum(res$rel >= 1e-8, na.rm = TRUE), "no reading", sum(is.na(res$rel)), "\n")
print(res[is.na(res$rel) | res$rel >= 1e-8, ], right = FALSE)

# Item 3: every reader over specs with no cs(), no me(), no layered family, on the
# lane's build and rellib-r16. Writes a list keyed spec|reader.
#   REV_LIB=before Rscript dev/review-predread-regress.R <out.rds>
# Data seed 7801. Fits are made in this process; the comparison script reports
# bitwise equality and the largest relative difference.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
outf <- commandArgs(TRUE)[1]
set.seed(7801)
n <- 240
d <- data.frame(x = rnorm(n), z = rnorm(n), f = factor(sample(c("a", "b", "c"), n, TRUE)),
                g = factor(rep(sprintf("g%02d", 1:12), each = 20)))
re <- rnorm(12, 0, 0.6)[d$g]
eta <- 0.4 + 0.6 * d$x - 0.3 * d$z + re
d$y <- eta + rnorm(n, 0, 0.7)
d$y2 <- -0.3 + 0.4 * d$x + rnorm(n, 0, 0.6)
d$cnt <- rpois(n, exp(eta))
d$bin <- rbinom(n, 1, plogis(eta))
d$pos <- exp(eta + rnorm(n, 0, 0.4))
d$frac <- rbeta(n, plogis(eta) * 8, (1 - plogis(eta)) * 8)
d$succ <- rbinom(n, 10, plogis(eta))
d$tr <- 10L
d$zcnt <- ifelse(runif(n) < 0.3, 0L, d$cnt)
d$hpos <- ifelse(runif(n) < 0.25, 0, d$pos)
d$ord <- as.integer(cut(eta + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf), labels = FALSE))
d$m <- sample(0:3, n, TRUE)
d$off <- runif(n, 0, 0.5)
d$t <- rep(1:20, 12)
d$w <- runif(n, 0.5, 1.5)
d$xm <- d$x
d$xm[seq(5, n, by = 17)] <- NA
nd <- data.frame(x = c(-1, 0, 0.7, 1.5), z = c(0.2, -0.5, 1, 0),
                 f = factor(c("a", "b", "c", "a"), levels = levels(d$f)),
                 g = factor(c("g01", "g02", "g03", "g01"), levels = levels(d$g)),
                 m = c(0L, 1L, 2L, 3L), off = 0.25, t = c(1L, 2L, 3L, 4L), tr = 10L,
                 y = 0.1, w = 1, xm = c(-1, 0, 0.7, 1.5))
specs <- list(
  gaussian = list(bf(y ~ x + z), gaussian()),
  gaussian_w = list(bf(y | weights(w) ~ x), gaussian()),
  group = list(bf(y ~ x + (1 + x | g)), gaussian()),
  factor = list(bf(y ~ x * f), gaussian()),
  offset = list(bf(y ~ x + offset(off)), gaussian()),
  poly = list(bf(y ~ poly(x, 2) + scale(z)), gaussian()),
  mo = list(bf(y ~ mo(m) + x), gaussian()),
  smooth = list(bf(y ~ s(x, k = 5) + z), gaussian()),
  gp = list(bf(y ~ gp(x, k = 8, c = 5/4)), gaussian()),
  sigma = list(bf(y ~ x, sigma ~ z), gaussian()),
  nl = list(bf(y ~ a + b * x, a ~ 1 + z, b ~ 1, nl = TRUE), gaussian()),
  ar_cov = list(bf(y ~ x + ar(time = t, gr = g, cov = TRUE)), gaussian()),
  ar_cond = list(bf(y ~ x + ar(time = t, gr = g, cov = FALSE)), gaussian()),
  student = list(bf(y ~ x), student()),
  poisson = list(bf(cnt ~ x + (1 | g)), poisson()),
  negbin = list(bf(cnt ~ x, shape ~ z), negbinomial()),
  bernoulli = list(bf(bin ~ x + z), bernoulli()),
  binomial = list(bf(succ | trials(tr) ~ x), binomial()),
  gamma = list(bf(pos ~ x), Gamma(link = "log")),
  lognormal = list(bf(pos ~ x), lognormal()),
  beta = list(bf(frac ~ x, phi ~ z), Beta()),
  zip = list(bf(zcnt ~ x, zi ~ z), zero_inflated_poisson()),
  hurdle = list(bf(hpos ~ x, hu ~ z), hurdle_lognormal()),
  cumulative = list(bf(ord ~ x + z), cumulative()),
  sratio_nocs = list(bf(ord ~ x), sratio()),
  mv = list(mvbf(bf(y ~ x), bf(y2 ~ x + z), rescor = FALSE), gaussian()),
  mv_rescor = list(mvbf(bf(y ~ x), bf(y2 ~ x + z), rescor = TRUE), gaussian()),
  mi = list(bf(y ~ mi(xm) + z) + bf(xm | mi() ~ z) + set_rescor(FALSE), gaussian()))
out <- list()
num <- function(e) {
  v <- tryCatch(quiet(e), error = function(err) {
    paste("ERROR:", substr(conditionMessage(err), 1, 200))
  })
  v
}
for (nm in names(specs)) {
  cat("== ", nm, "\n")
  fit <- tryCatch(quiet(frm(specs[[nm]][[1L]], family = specs[[nm]][[2L]], data = d)),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    out[[paste0(nm, "|frm")]] <- paste("ERROR:", conditionMessage(fit))
    next
  }
  mv <- length(fit$spec$responses) > 1L
  ord <- nm %in% c("cumulative", "sratio_nocs")
  R <- list(
    linpred_link = function() frm_linpred(fit, type = "link"),
    linpred_link_nd = function() frm_linpred(fit, newdata = nd, type = "link"),
    linpred_se = function() frm_linpred(fit, newdata = nd, se.fit = TRUE),
    linpred_resp = function() frm_linpred(fit, newdata = nd, type = "response"),
    fitted = function() fitted(fit),
    fitted_nd = function() fitted(fit, newdata = nd),
    fitted_linear = function() fitted(fit, newdata = nd, scale = "linear"),
    predict = function() {set.seed(3); predict(fit, ndraws = 30)},
    predict_nd = function() {set.seed(3); predict(fit, newdata = nd, ndraws = 30)},
    predict_nop = function() {set.seed(3); predict(fit, newdata = nd, ndraws = 30, propagate_error = FALSE)},
    simulate = function() as.matrix(simulate(fit, nsim = 3, seed = 3)),
    simulate_nd = function() as.matrix(simulate(fit, nsim = 3, seed = 3, newdata = nd)),
    residuals = function() residuals(fit),
    basis = function() frm_lp_basis(fit, newdata = nd, extra_cov = TRUE)[c("eta", "A", "extra_var")],
    curve = function() as.data.frame(frm_curve(fit, newdata = nd, seed = 1)),
    curve_tr = function() as.data.frame(frm_curve(fit, newdata = nd, seed = 1, transform = TRUE)),
    curve_contrast = function() as.data.frame(frm_curve(fit, newdata = nd[1:2, ], contrast = nd[3:4, ], seed = 1)),
    curve_deriv = function() as.data.frame(frm_curve_deriv(fit, "x", newdata = nd, seed = 1)),
    emmeans = function() as.data.frame(emmeans::emmeans(fit, ~ x, at = list(x = c(-1, 1)))),
    get_predict = function() marginaleffects::get_predict(fit, newdata = nd, re_formula = NA)$estimate,
    get_predict_link = function() marginaleffects::get_predict(fit, newdata = nd, type = "link", re_formula = NA)$estimate,
    avg_pred = function() {
      r <- marginaleffects::avg_predictions(fit, newdata = nd, re_formula = NA)
      c(r$estimate, r$std.error)
    },
    CE = function() {
      lapply(conditional_effects(fit, resolution = 4), function(p) {
        as.data.frame(p)[, c("estimate__", "se__", "lower__", "upper__")]
      })
    },
    CE_predict = function() {
      lapply(conditional_effects(fit, method = "predict", resolution = 4, seed = 4), function(p) {
        as.data.frame(p)[, c("estimate__", "lower__", "upper__")]
      })
    },
    log_lik = function() log_lik(fit),
    insight = function() as.numeric(insight::get_predicted(fit, data = nd)))
  if (!ord) {
    R$linpred_mu <- function() frm_linpred(fit, newdata = nd, dpar = if (mv) NULL else "mu", type = "response")
    R$fitted_mu <- function() fitted(fit, newdata = nd, dpar = if (mv) NULL else "mu")
  }
  if (mv) {
    for (r in c("y", "y2")) {
      local({
        rr <- r
        R[[paste0("linpred_resp_", rr)]] <<- function() frm_linpred(fit, newdata = nd, resp = rr, se.fit = TRUE)
        R[[paste0("fitted_", rr)]] <<- function() fitted(fit, newdata = nd, resp = rr)
        R[[paste0("curve_", rr)]] <<- function() as.data.frame(frm_curve(fit, newdata = nd, resp = rr, seed = 1))
      })
    }
  }
  ds <- rev_draws_at(fit, 5L)
  R$draws_linpred <- function() frmtmb.sample::posterior_linpred(ds)
  R$draws_linpred_nd <- function() frmtmb.sample::posterior_linpred(ds, newdata = nd)
  R$draws_epred <- function() frmtmb.sample::posterior_epred(ds, newdata = nd)
  R$draws_pp <- function() {
    set.seed(5)
    frmtmb.sample::posterior_predict(ds, newdata = nd)
  }
  R$draws_fitted <- function() frmtmb.sample:::fitted.frmtmb_draws(ds, newdata = nd, summary = FALSE)
  R$draws_loglik <- function() frmtmb.sample::log_lik(ds)
  for (r in names(R)) {
    out[[paste0(nm, "|", r)]] <- num(R[[r]]())
  }
}
saveRDS(out, outf)
cat("keys", length(out), "\n")

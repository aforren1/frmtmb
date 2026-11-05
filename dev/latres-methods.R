# Post-fit methods on latent-residual fits against brms 2.23.0 at the same
# parameter values (Stan's Fixed_param, brms_fixed_fit()). Prints, per shape,
# what each brms method does with the latent residual and whether frmtmb's
# method agrees. Data seed 167.
#
# Usage: Rscript dev/latres-methods.R [shapes]
args <- commandArgs(TRUE)
setwd("C:/Users/adf44/source/r/frmtmb-wt-latres")
.libPaths(c("C:/Users/adf44/source/r/wt-latres-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-latres/dev/stan-cache")
suppressPackageStartupMessages({
  library(frmtmb)
  library(rstan)
})
hel <- new.env(parent = asNamespace("frmtmb"))
hel$`%||%` <- function(a, b) if (is.null(a)) b else a
sys.source("dev/equiv-lib.R", envir = hel)

set.seed(167)
d <- expand.grid(t = 1:5, g = factor(1:30))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
e <- as.vector(replicate(30, arima.sim(list(ar = 0.5), 5, sd = 0.5)))
eta <- 0.8 + 0.4 * d$x + e
d$count <- rpois(n, exp(eta))
d$bin <- rbinom(n, 1, plogis(eta - 0.8))
d$pos <- rgamma(n, 4, 4 / exp(eta))

shapes <- list(
  ar1t = list("count ~ x + ar(t, g, p = 1, cov = TRUE)", "poisson()"),
  ar1f = list("count ~ x + ar(t, g, p = 1)", "poisson()"),
  ma1t = list("bin ~ x + ma(t, g, q = 1, cov = TRUE)", "bernoulli()"),
  gamf = list("pos ~ x + ar(t, g, p = 1)", "Gamma(link = 'log')"),
  cosy = list("count ~ x + cosy(t, g)", "poisson()"),
  zre = list("count ~ x + (1 | g) + ar(t, g, cov = TRUE)", "poisson()"))
if (length(args)) shapes <- shapes[strsplit(args[[1]], ",")[[1]]]
rel <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  if (length(a) != length(b)) return(sprintf("lengths %d %d", length(a), length(b)))
  sprintf("%.3g", max(abs(a - b) / pmax(1, abs(a), abs(b))))
}
tr <- function(expr) tryCatch(expr, error = function(e) {
  structure(conditionMessage(e), class = "err")
})
say <- function(lab, a, b) {
  if (inherits(a, "err") || inherits(b, "err")) {
    cat(sprintf("%-28s frm: %s | brms: %s\n", lab,
                if (inherits(a, "err")) substr(a, 1, 110) else "ok",
                if (inherits(b, "err")) substr(b, 1, 110) else "ok"))
  } else {
    cat(sprintf("%-28s max rel diff %s\n", lab, rel(a, b)))
  }
}
est <- function(x) if (length(dim(x)) == 2L) x[, 1L] else x
mz <- function(a, b) {
  # per-row mean z over two sets of draws, nrow = draws
  se <- sqrt(apply(a, 2L, var) / nrow(a) + apply(b, 2L, var) / nrow(b))
  z <- (colMeans(a) - colMeans(b)) / se
  sprintf("max |z| %.2f over %d rows", max(abs(z[se > 0])), ncol(a))
}
for (nm in names(shapes)) {
  s <- shapes[[nm]]
  cat("\n=====", nm, s[[1]], s[[2]], "\n")
  f <- as.formula(s[[1]])
  fit <- suppressWarnings(frm(frmtmb::bf(f),
                              family = hel$eq_eval(s[[2]], "frmtmb"),
                              data = d))
  fb <- brms::bf(f)
  famb <- hel$eq_eval(s[[2]], "brms")
  bfit <- suppressWarnings(hel$brms_fixed_fit(fb, famb, d, fit, ndraws = 10))
  bfit4 <- suppressWarnings(hel$brms_fixed_fit(fb, famb, d, fit, ndraws = 2000))
  print(utils::head(brms::variables(bfit), 12))
  print(tail(brms::variables(bfit), 4))
  ep <- brms::posterior_epred(bfit)[1L, ]
  say("fitted", tr(est(fitted(fit))), ep)
  say("fitted re_formula = NA", tr(est(fitted(fit, re_formula = NA))),
      tr(brms::posterior_epred(bfit, re_formula = NA)[1L, ]))
  say("linpred", tr(est(frm_linpred(fit))),
      tr(brms::posterior_linpred(bfit)[1L, ]))
  say("residuals", tr(est(residuals(fit, type = "response"))),
      tr(residuals(bfit, method = "posterior_epred", type = "ordinary",
                   summary = FALSE)[1L, ]))
  say("log_lik", tr(hel$frm_row_loglik(fit)), tr(brms::log_lik(bfit)[1L, ]))
  say("log_lik()", tr(as.numeric(log_lik(fit))), tr(brms::log_lik(bfit)[1L, ]))
  pp <- tr(predict(fit, summary = FALSE, ndraws = 2000,
                   propagate_error = FALSE))
  ppb <- tr(brms::posterior_predict(bfit4))
  if (!inherits(pp, "err") && !inherits(ppb, "err")) {
    cat("predict in sample: ", mz(pp, ppb), "\n")
  } else say("predict", pp, ppb)
  ss <- tr({
    sim <- simulate(fit, nsim = 2000, seed = 7)
    t(vapply(sim, as.numeric, numeric(nrow(sim))))
  })
  if (!inherits(ss, "err") && !inherits(ppb, "err")) {
    cat("simulate in sample:", mz(ss, ppb), "\n")
  } else say("simulate", ss, ppb)
  nd <- d[d$g %in% levels(d$g)[1:2], ]
  nd$g <- factor(paste0("n", nd$g))
  fn <- tr(est(fitted(fit, newdata = nd, allow_new_levels = TRUE)))
  fb_nd <- tr(brms::posterior_epred(bfit4, newdata = nd,
                                    allow_new_levels = TRUE))
  cat("fitted newdata frm: ",
      if (inherits(fn, "err")) substr(fn, 1, 150) else "ok",
      " brms: ", if (inherits(fb_nd, "err")) substr(fb_nd, 1, 150) else {
        sprintf("draws differ: sd over draws of row 1 = %.4g",
                sd(fb_nd[, 1L]))
      }, "\n")
  pn <- tr(predict(fit, newdata = nd, allow_new_levels = TRUE, summary = FALSE,
                   ndraws = 2000, propagate_error = FALSE))
  pb_nd <- tr(brms::posterior_predict(bfit4, newdata = nd,
                                      allow_new_levels = TRUE))
  if (!inherits(pn, "err") && !inherits(pb_nd, "err")) {
    cat("predict newdata:", mz(pn, pb_nd), "\n")
  } else say("predict newdata", pn, pb_nd)
  sn <- tr({
    sim <- simulate(fit, nsim = 2000, seed = 7, newdata = nd,
                    allow_new_levels = TRUE)
    t(vapply(sim, as.numeric, numeric(nrow(sim))))
  })
  if (!inherits(sn, "err") && !inherits(pb_nd, "err")) {
    cat("simulate newdata:", mz(sn, pb_nd), "\n")
  } else say("simulate newdata", sn, pb_nd)
  say("emmeans", tr(as.data.frame(emmeans::emmeans(fit, ~ x))[[2L]]),
      tr(as.data.frame(emmeans::emmeans(bfit, ~ x))[[2L]]))
  say("cond_effects",
      tr(conditional_effects(fit, "x", method = "posterior_epred")[[1L]]$estimate__),
      tr(brms::conditional_effects(bfit, "x", method = "posterior_epred")[[1L]]$estimate__))
  cat("\nsummary:\n")
  print(tr(summary(fit)))
  print(tr(confint_varcorr(fit)))
  print(tr(VarCorr(fit)))
  print(tr(names(ranef(fit))))
  print(tr(prior_summary(fit)))
  print(tr(default_prior(fit)))
  print(tr(variables(fit)))
  print(utils::capture.output(print(brms::prior_summary(bfit)))[1:8])
}

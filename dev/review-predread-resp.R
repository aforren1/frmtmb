# Item 6: resp readers. Data seed 7501.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
library(ggplot2)
set.seed(7501)
n <- 150
d <- data.frame(x = rnorm(n), z = rnorm(n), f = factor(sample(c("a", "b"), n, TRUE)))
d$y <- 1 + d$x + rnorm(n)
d$y2 <- -d$x + 0.5 * d$z + rnorm(n)
comp <- rbinom(n, 1, 0.5)
d$ym <- ifelse(comp == 1, 3 + d$x, -3 - d$x) + rnorm(n, 0, 0.7)
U <- rev_pair("bf(y ~ x)", "gaussian()", d)
ds <- rev_draws_at(U$fit, 8L)
vals <- list(empty = character(0), twice = c("y", "y"), unknown = "nosuch",
             two_unknown = c("a", "b"), na = NA_character_, num = 1L, null = NULL)
nd <- data.frame(x = c(-1, 1))
fmt <- function(e) {
  r <- tryCatch(quiet(e), error = function(err) structure(list(), m = conditionMessage(err)))
  if (!is.null(attr(r, "m"))) return(paste("ERR:", substr(attr(r, "m"), 1, 90)))
  if (is.list(r)) return(paste("list", paste(names(r), collapse = ",")))
  paste("OK", paste(dim(r) %||% length(r), collapse = "x"))
}
cat("=== univariate: every method, every resp value\n")
meths <- list(
  brms_epred = function(r) posterior_epred(U$fx, newdata = nd, resp = r),
  brms_linpred = function(r) posterior_linpred(U$fx, newdata = nd, resp = r),
  brms_pp = function(r) posterior_predict(U$fx, newdata = nd, resp = r),
  brms_fitted = function(r) fitted(U$fx, newdata = nd, resp = r),
  brms_predict = function(r) predict(U$fx, newdata = nd, resp = r),
  brms_resid = function(r) residuals(U$fx, resp = r),
  brms_loglik = function(r) log_lik(U$fx, resp = r),
  brms_R2 = function(r) bayes_R2(U$fx, resp = r),
  brms_CE = function(r) brms::conditional_effects(U$fx, "x", resp = r),
  brms_pperr = function(r) predictive_error(U$fx, resp = r),
  draws_epred = function(r) frmtmb.sample::posterior_epred(ds, newdata = nd, resp = r),
  draws_linpred = function(r) frmtmb.sample::posterior_linpred(ds, newdata = nd, resp = r),
  draws_pp = function(r) frmtmb.sample::posterior_predict(ds, newdata = nd, resp = r),
  draws_fitted = function(r) frmtmb.sample:::fitted.frmtmb_draws(ds, newdata = nd, resp = r),
  draws_predict = function(r) frmtmb.sample:::predict.frmtmb_draws(ds, newdata = nd, resp = r),
  draws_resid = function(r) frmtmb.sample:::residuals.frmtmb_draws(ds, resp = r),
  draws_loglik = function(r) frmtmb.sample::log_lik(ds, resp = r),
  draws_R2 = function(r) frmtmb.sample::bayes_R2(ds, resp = r),
  draws_CE = function(r) conditional_effects(ds, "x", resp = r),
  draws_pperr = function(r) frmtmb.sample::predictive_error(ds, resp = r),
  draws_loo_R2 = function(r) frmtmb.sample::loo_R2(ds, resp = r),
  fit_fitted = function(r) fitted(U$fit, newdata = nd, resp = r),
  fit_predict = function(r) predict(U$fit, newdata = nd, resp = r, ndraws = 5),
  fit_resid = function(r) residuals(U$fit, resp = r),
  fit_linpred = function(r) frm_linpred(U$fit, newdata = nd, resp = r),
  fit_basis = function(r) frm_lp_basis(U$fit, newdata = nd, resp = r)$eta,
  fit_curve = function(r) frm_curve(U$fit, newdata = nd, resp = r, simultaneous = FALSE)$.se,
  fit_curve_deriv = function(r) frm_curve_deriv(U$fit, "x", newdata = nd, resp = r)$.se,
  fit_CE = function(r) conditional_effects(U$fit, "x", resp = r),
  fit_emm = function(r) as.data.frame(emmeans::emmeans(U$fit, ~ x, at = list(x = 1), resp = r))$SE,
  fit_ppcheck = function(r) pp_check(U$fit, resp = r, ndraws = 3),
  fit_simulate = function(r) as.matrix(simulate(U$fit, nsim = 1, seed = 1, resp = r)),
  fit_loglik = function(r) log_lik(U$fit, resp = r),
  fit_pperr = function(r) predictive_error(U$fit, resp = r, ndraws = 3),
  fit_predint = function(r) predictive_interval(U$fit, resp = r, ndraws = 3),
  fit_R2 = function(r) bayes_R2(U$fit, resp = r)
)
for (m in names(meths)) {
  cat(sprintf("%-16s", m))
  for (v in names(vals)) cat(" |", v, ":", fmt(meths[[m]](vals[[v]])))
  cat("\n")
}
cat("\n=== pp_mixture univariate\n")
Mx <- tryCatch(rev_pair("bf(ym ~ x)", "mixture(gaussian(), gaussian())", d), error = function(e) e)
if (!inherits(Mx, "error")) {
  dm <- rev_draws_at(Mx$fit, 4L)
  for (v in names(vals)) {
    cat(v, ": brms", fmt(pp_mixture(Mx$fx, resp = vals[[v]])), "| fit",
        fmt(pp_mixture(Mx$fit, resp = vals[[v]])), "| draws",
        fmt(frmtmb.sample::pp_mixture(dm, resp = vals[[v]])), "\n")
  }
}
cat("\n=== multivariate\n")
MV <- rev_pair("mvbf(bf(y ~ x), bf(y2 ~ x + z + f), rescor = FALSE)", "gaussian()", d)
ms <- rev_draws_at(MV$fit, 8L)
fit <- MV$fit
fx <- MV$fx
nms <- function(e) {
  r <- tryCatch(quiet(e), error = function(err) structure(list(), m = conditionMessage(err)))
  if (!is.null(attr(r, "m"))) return(paste("ERR:", substr(attr(r, "m"), 1, 100)))
  paste(names(r), collapse = " ")
}
for (lab in c("NULL", "empty", "y2", "both", "rev", "twice", "unknown")) {
  r <- switch(lab, "NULL" = NULL, empty = character(0), y2 = "y2", both = c("y", "y2"),
              rev = c("y2", "y"), twice = c("y2", "y2"), unknown = "nosuch")
  cat(sprintf("%-8s brms  : %s\n", lab, nms(brms::conditional_effects(fx, "x", resp = r))))
  cat(sprintf("%-8s fit   : %s\n", lab, nms(conditional_effects(fit, "x", resp = r, resolution = 5))))
  cat(sprintf("%-8s draws : %s\n", lab, nms(conditional_effects(ms, "x", resp = r, resolution = 5))))
}
cat("\n--- effects = NULL, all terms\n")
cat("brms  :", nms(brms::conditional_effects(fx)), "\n")
cat("fit   :", nms(conditional_effects(fit, resolution = 4)), "\n")
cat("draws :", nms(conditional_effects(ms, resolution = 4)), "\n")
cat("\n--- effects = 'z' read by one response only\n")
cat("brms  :", nms(brms::conditional_effects(fx, "z")), "\n")
cat("fit   :", nms(conditional_effects(fit, "z", resolution = 4)), "\n")
cat("draws :", nms(conditional_effects(ms, "z", resolution = 4)), "\n")
cat("\n--- effects = 'x:f'\n")
cat("brms  :", nms(brms::conditional_effects(fx, "x:f")), "\n")
cat("fit   :", nms(conditional_effects(fit, "x:f", resolution = 4)), "\n")
cat("draws :", nms(conditional_effects(ms, "x:f", resolution = 4)), "\n")
cat("\n--- conditions\n")
cond <- data.frame(z = c(-1, 1))
cb <- brms::conditional_effects(fx, "x", conditions = cond, resolution = 4)
cf <- conditional_effects(fit, "x", conditions = cond, resolution = 4)
cd <- conditional_effects(ms, "x", conditions = cond, resolution = 4)
cat("brms:", names(cb), " rows ", sapply(cb, nrow), " cond__ ",
    paste(levels(factor(cb[[2]]$cond__)), collapse = ","), "\n")
cat("fit :", names(cf), " rows ", sapply(cf, nrow), " cond__ ",
    paste(levels(factor(cf[[2]]$cond__)), collapse = ","), "\n")
cat("draws:", names(cd), " rows ", sapply(cd, nrow), "\n")
cat("cols brms :", paste(names(cb[[2]]), collapse = " "), "\n")
cat("cols fit  :", paste(names(cf[[2]]), collapse = " "), "\n")
cat("cols draws:", paste(names(cd[[2]]), collapse = " "), "\n")
cat("\n--- estimates by panel against brms, each panel holds its response\n")
pb <- brms::conditional_effects(fx, "x", resolution = 5)
pf <- conditional_effects(fit, "x", resolution = 5)
pd <- conditional_effects(ms, "x", resolution = 5)
for (k in names(pb)) {
  cat(k, " brms", round(pb[[k]]$estimate__, 4), "\n", k, " fit ",
      round(pf[[k]]$estimate__, 4), "\n", k, " draw", round(pd[[k]]$estimate__, 4), "\n")
  cat("  attr effects brms:", paste(attr(pb[[k]], "effects"), collapse = ","),
      " fit:", paste(attr(pf[[k]], "effects"), collapse = ","), "\n")
  cat("  response brms:", attr(pb[[k]], "response"), " fit:", attr(pf[[k]], "response"),
      " draws:", attr(pd[[k]], "response"), "\n")
}
cat("\n--- method = predict / posterior_predict\n")
cat("brms posterior_predict:", nms(brms::conditional_effects(fx, "x", method = "posterior_predict", resolution = 4)), "\n")
cat("fit predict           :", nms(conditional_effects(fit, "x", method = "predict", resolution = 4, seed = 1)), "\n")
cat("fit posterior_predict :", nms(conditional_effects(fit, "x", method = "posterior_predict", resolution = 4, seed = 1)), "\n")
cat("draws predict         :", nms(conditional_effects(ms, "x", method = "predict", resolution = 4)), "\n")
cat("draws posterior_predict:", nms(conditional_effects(ms, "x", method = "posterior_predict", resolution = 4)), "\n")
cat("\n--- plots: points = TRUE; ggplot label and points by panel\n")
pl <- function(ce, label, ...) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  ps <- tryCatch(quiet(plot(ce, plot = FALSE, points = TRUE, ...)), error = function(e) e)
  if (inherits(ps, "error")) {
    cat(label, "plot ERR:", conditionMessage(ps), "\n")
    return(invisible())
  }
  cat(label, "plots:", paste(names(ps), collapse = " "), "\n")
  for (k in names(ps)) {
    g <- ps[[k]]
    b <- tryCatch(ggplot2::ggplot_build(g), error = function(e) e)
    if (inherits(b, "error")) {
      cat("  ", k, " build ERR:", conditionMessage(b), "\n")
      next
    }
    isp <- which(sapply(g$layers, function(l) inherits(l$geom, "GeomPoint")))
    np <- if (length(isp)) nrow(b$data[[isp[1L]]]) else 0
    py <- if (length(isp)) round(mean(b$data[[isp[1L]]]$y), 4) else NA
    cat("  ", k, ": ylab ", as.character(g$labels$y %||% ""), " points ", np,
        " mean y of points ", py, "\n")
  }
}
cat("means of data: y", round(mean(d$y), 4), " y2", round(mean(d$y2), 4), "\n")
pl(brms::conditional_effects(fx, "x", resolution = 4), "brms")
pl(conditional_effects(fit, "x", resolution = 4), "fit")
pl(conditional_effects(ms, "x", resolution = 4), "draws")
cat("\n--- one named response, a variable it does not read\n")
cat("brms:", nms(brms::conditional_effects(fx, "z", resp = "y")), "\n")
cat("fit :", nms(conditional_effects(fit, "z", resp = "y", resolution = 4)), "\n")
cat("draws:", nms(conditional_effects(ms, "z", resp = "y", resolution = 4)), "\n")
cat("\n--- spaghetti on mv\n")
cat("brms spaghetti:", nms(brms::conditional_effects(fx, "x", spaghetti = TRUE, ndraws = 5)), "\n")
cat("fit spaghetti:", nms(conditional_effects(fit, "x", spaghetti = TRUE, ndraws = 5, resolution = 4)), "\n")
cat("draws spaghetti:", nms(conditional_effects(ms, "x", spaghetti = TRUE, ndraws = 5, resolution = 4)), "\n")
cat("\n--- pp_mixture multivariate\n")
d$ym2 <- ifelse(comp == 1, 1 + d$z, -2 - d$z) + rnorm(n, 0, 0.7)
ftxt <- paste0("mvbf(bf(ym ~ x, family = mixture(gaussian(), gaussian())), ",
               "bf(ym2 ~ z, family = mixture(gaussian(), gaussian())), rescor = FALSE)")
fb <- rev_ev(ftxt, "brms")
ff <- rev_ev(ftxt, "frmtmb")
f0 <- tryCatch(quiet(frm(ff, data = d)), error = function(e) e)
if (inherits(f0, "error")) {
  cat("mv mixture fit ERR:", conditionMessage(f0), "\n")
} else {
  b0 <- tryCatch(quiet(rev_E$brms_fixed_fit(fb, NULL, d, f0, ndraws = 4)), error = function(e) e)
  dmm <- rev_draws_at(f0, 4L)
  for (lab in c("NULL", "empty", "ym2", "both", "twice", "unknown")) {
    r <- switch(lab, "NULL" = NULL, empty = character(0), ym2 = "ym2",
                both = c("ym", "ym2"), twice = c("ym2", "ym2"), unknown = "nosuch")
    if (!inherits(b0, "error")) cat(lab, ": brms", fmt(pp_mixture(b0, resp = r)), "\n")
    cat(lab, ": fit", fmt(pp_mixture(f0, resp = r)), "| draws",
        fmt(frmtmb.sample::pp_mixture(dmm, resp = r)), "\n")
  }
  if (!inherits(b0, "error")) {
    for (r in c("ym", "ym2")) {
      a <- pp_mixture(b0, resp = r)
      f <- pp_mixture(f0, resp = r)
      g <- frmtmb.sample::pp_mixture(dmm, resp = r)
      cat(r, "first-row probs: brms", round(a[1, "Estimate", ], 4), " fit",
          round(f[1, "Estimate", ], 4), " draws", round(g[1, "Estimate", ], 4), "\n")
    }
  } else cat("mv mixture brms ERR:", conditionMessage(b0), "\n")
}
cat("\n--- bayes_R2 and loo_R2, multivariate, stored criteria\n")
sampled <- quiet(frm_sample(frmtmb::mvbf(bf(y ~ x), bf(y2 ~ x), rescor = FALSE),
                            family = gaussian(), data = d, chains = 1, iter = 300,
                            seed = 2, refresh = 0))
bsamp <- quiet(brms::brm(brms::mvbf(brms::bf(y ~ x), brms::bf(y2 ~ x), rescor = FALSE),
                         data = d, chains = 1, iter = 300, seed = 2, refresh = 0, silent = 2))
stor <- function(obj, fn, ...) {
  for (lab in c("NULL", "empty", "y2", "both", "twice")) {
    r <- switch(lab, "NULL" = NULL, empty = character(0), y2 = "y2",
                both = c("y", "y2"), twice = c("y2", "y2"))
    v <- tryCatch(quiet(fn(obj, resp = r)),
                  error = function(e) paste("ERR:", substr(conditionMessage(e), 1, 90)))
    cat("  ", lab, ":", if (is.matrix(v)) paste(paste(rownames(v), collapse = ","), paste(dim(v), collapse = "x")) else v, "\n")
  }
}
cat("brms bayes_R2, not stored\n")
stor(bsamp, brms::bayes_R2)
cat("draws bayes_R2, not stored\n")
stor(sampled, frmtmb.sample::bayes_R2)
cat("brms loo_R2, not stored\n")
stor(bsamp, brms::loo_R2)
cat("draws loo_R2, not stored\n")
stor(sampled, frmtmb.sample::loo_R2)
bs2 <- quiet(brms::add_criterion(bsamp, c("bayes_R2", "loo_R2")))
fs2 <- tryCatch(quiet(frmtmb.sample::add_criterion(sampled, c("bayes_R2", "loo_R2"))),
                error = function(e) e)
if (inherits(fs2, "error")) {
  cat("add_criterion draws ERR:", conditionMessage(fs2), "\n")
} else {
  cat("brms bayes_R2, stored\n")
  stor(bs2, brms::bayes_R2)
  cat("draws bayes_R2, stored\n")
  stor(fs2, frmtmb.sample::bayes_R2)
  cat("brms loo_R2, stored\n")
  stor(bs2, brms::loo_R2)
  cat("draws loo_R2, stored\n")
  stor(fs2, frmtmb.sample::loo_R2)
  cat("stored and fresh on empty: equal?",
      isTRUE(all.equal(frmtmb.sample::bayes_R2(fs2, resp = character(0)),
                       frmtmb.sample::bayes_R2(fs2))), "\n")
}
cat("\n--- univariate draws, stored, empty resp\n")
us <- quiet(frm_sample(frmtmb::bf(y ~ x), family = gaussian(), data = d, chains = 1,
                       iter = 300, seed = 2, refresh = 0))
us2 <- tryCatch(quiet(frmtmb.sample::add_criterion(us, c("bayes_R2", "loo_R2"))),
                error = function(e) e)
if (!inherits(us2, "error")) {
  for (r in list(NULL, character(0), "nosuch", c("y", "y"))) {
    cat(paste(r, collapse = ","), ": bayes_R2", fmt(frmtmb.sample::bayes_R2(us2, resp = r)),
        " loo_R2", fmt(frmtmb.sample::loo_R2(us2, resp = r)), "\n")
  }
}

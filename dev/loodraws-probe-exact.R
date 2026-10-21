# Exactness of the loo_*() functions on draws: each against loo::E_loo()
# and loo::psis() called by hand on frmtmb's own log_lik() and
# prediction matrices, and against brms 2.23.0's own internal functions
# (E_loo_value(), .loo_R2()) on the same matrices.
#   Rscript dev/loodraws-probe-exact.R [lib]
# Log: dev/loodraws-log/probe-exact.txt (lane lib) and
# dev/loodraws-log/probe-exact-base.txt (rellib-r10, the before arm).
a <- commandArgs(trailingOnly = TRUE)
lib <- if (length(a)) a[1] else "C:/Users/adf44/source/r/wt-loodraws-lib"
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-loodraws/dev/stan-cache")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", dirname(system.file(package = "frmtmb.sample")),
    format(packageVersion("frmtmb.sample")), "\n")
set.seed(1)
n <- 60
dd <- data.frame(x = rnorm(n), g = factor(rep(1:10, 6)))
dd$y <- rnorm(n, 1 + 0.5 * dd$x + rnorm(10, 0, 0.7)[dd$g])
dd$y2 <- rnorm(n, -1 + 0.3 * dd$x)
dd$k <- rpois(n, exp(0.2 + 0.4 * dd$x))
smp <- function(...) suppressWarnings(suppressMessages(
  frm_sample(..., chains = 2, iter = 600, warmup = 300, seed = 11,
             refresh = 0)))
ds <- list(
  gaussian = smp(bf(y ~ x), family = gaussian(), data = dd),
  poisson = smp(bf(k ~ x), family = poisson(), data = dd),
  ranef = smp(bf(y ~ x + (1 | g)), family = gaussian(), data = dd),
  mvbf = smp(mvbf(bf(y ~ x), bf(y2 ~ x), rescor = FALSE),
             family = gaussian(), data = dd))
q <- function(e) suppressWarnings(suppressMessages(e))
res <- function(lab, ok) cat(sprintf("%-46s %s\n", lab, ok))
for (nm in names(ds)) {
  d <- ds[[nm]]
  out <- tryCatch({
    resp <- if (nm == "mvbf") "y2" else NULL
    ll <- log_lik(d, resp = resp)
    ps <- q(loo::psis(-ll, r_eff = loo::relative_eff(
      exp(ll), chain_id = attr(ll, "chain_id"))))
    ep <- posterior_epred(d, resp = resp)
    lp <- posterior_linpred(d, resp = resp)
    res(paste(nm, "psis(ds) == by hand"),
        identical(stats::weights(q(psis(d, resp = resp))),
                  stats::weights(ps)))
    res(paste(nm, "loo(save_psis)$psis_object == psis(ds)"),
        identical(stats::weights(q(loo(d, resp = resp, save_psis = TRUE))
                                 $psis_object),
                  stats::weights(ps)))
    for (ty in c("mean", "var", "quantile")) {
      pr <- if (ty == "quantile") c(0.1, 0.5, 0.9) else 0.5
      lty <- if (ty == "var") "variance" else ty
      ref <- loo::E_loo(ep, ps, type = lty, probs = pr)$value
      ref <- if (is.matrix(ref)) t(ref) else matrix(ref)
      got <- q(loo_epred(d, type = ty, probs = pr, resp = resp))
      res(paste(nm, "loo_epred", ty), identical(unname(got), unname(ref)))
      refl <- loo::E_loo(lp, ps, type = lty, probs = pr)$value
      refl <- if (is.matrix(refl)) t(refl) else matrix(refl)
      gotl <- q(loo_linpred(d, type = ty, probs = pr, resp = resp))
      res(paste(nm, "loo_linpred", ty), identical(unname(gotl), unname(refl)))
      if (requireNamespace("brms", quietly = TRUE)) {
        bv <- brms:::E_loo_value(ep, ps, type = ty, probs = pr)
        res(paste(nm, "loo_epred", ty, "== brms E_loo_value"),
            identical(got, bv))
      }
    }
    # posterior_predict is random: one seed for both arms
    set.seed(3); pp <- posterior_predict(d, resp = resp)
    set.seed(3); gp <- q(loo_predict(d, type = "quantile",
                                     probs = c(0.05, 0.95), resp = resp))
    refp <- t(loo::E_loo(pp, ps, type = "quantile",
                         probs = c(0.05, 0.95))$value)
    res(paste(nm, "loo_predict quantile"), identical(unname(gp), refp))
    set.seed(3); gi <- q(loo_predictive_interval(d, prob = 0.9, resp = resp))
    # brms forms the probabilities as (1 - prob) / 2, which is not 0.05
    set.seed(3); gq <- q(loo_predict(d, type = "quantile", resp = resp,
                                     probs = c((1 - 0.9) / 2,
                                               1 - (1 - 0.9) / 2)))
    res(paste(nm, "loo_predictive_interval == loo_predict q"),
        identical(gi, gq))
    res(paste(nm, "interval colnames"), paste(colnames(gi), collapse = ","))
    # loo_R2: by hand with brms's arithmetic, then brms's own .loo_R2
    y <- frmtmb.sample:::draws_observed_y(d$fit,
                                          resp %||% names(d$fit$spec$responses)[1])
    r2 <- q(loo_R2(d, resp = if (nm == "mvbf") "y2", summary = FALSE,
                   seed = 7))
    if (requireNamespace("brms", quietly = TRUE)) {
      set.seed(7)
      br <- q(brms:::.loo_R2(as.numeric(y), ep, ll,
                             loo::relative_eff(exp(ll),
                                               chain_id = attr(ll, "chain_id"))))
      rel <- max(abs(r2[, 1] - br[, 1])) / max(abs(br[, 1]))
      res(paste(nm, "loo_R2 vs brms .loo_R2 max rel diff"),
          sprintf("%.3g (eps %.3g)", rel, .Machine$double.eps))
    }
    TRUE
  }, error = function(e) {
    cat(nm, "ERROR:", conditionMessage(e), "\n")
  })
}
# multivariate without resp: one layer per response, joint weights
d <- ds$mvbf
ll <- log_lik(d)
ps <- q(loo::psis(-ll, r_eff = loo::relative_eff(
  exp(ll), chain_id = attr(ll, "chain_id"))))
got <- tryCatch(q(loo_epred(d)), error = function(e) conditionMessage(e))
if (is.character(got)) cat("mv epred ERROR:", got, "\n") else {
  res("mv loo_epred dim", paste(dim(got), collapse = " x "))
  res("mv loo_epred dimnames[[3]]", paste(dimnames(got)[[3]], collapse = ","))
  for (r in c("y", "y2")) {
    ref <- matrix(loo::E_loo(posterior_epred(d, resp = r), ps)$value)
    res(paste("mv layer", r, "== E_loo joint weights"),
        identical(unname(got[, , r, drop = TRUE]), as.numeric(ref)))
  }
}
# loo_model_weights against loo's own on the same matrices
mw <- tryCatch(q(loo_model_weights(ds$gaussian, ds$ranef)),
               error = function(e) conditionMessage(e))
if (is.character(mw)) cat("lmw ERROR:", mw, "\n") else {
  lls <- list(log_lik(ds$gaussian), log_lik(ds$ranef))
  ref <- q(loo::loo_model_weights(lls, r_eff_list = lapply(lls, function(l)
    loo::relative_eff(exp(l), chain_id = attr(l, "chain_id")))))
  res("loo_model_weights stacking == loo's", identical(as.numeric(mw),
                                                       as.numeric(ref)))
  res("loo_model_weights names", paste(names(mw), collapse = ","))
  set.seed(5); mp <- q(loo_model_weights(ds$gaussian, ds$ranef,
                                         method = "pseudobma"))
  set.seed(5); rp <- q(loo::loo_model_weights(lls, method = "pseudobma",
    r_eff_list = lapply(lls, function(l)
      loo::relative_eff(exp(l), chain_id = attr(l, "chain_id")))))
  res("loo_model_weights pseudobma == loo's", identical(as.numeric(mp),
                                                        as.numeric(rp)))
  res("class", paste(class(mw), collapse = ","))
}
# pp_check's loo weights and psis(): same draws, same weights
if (requireNamespace("bayesplot", quietly = TRUE)) {
  d <- ds$ranef
  lw <- local({
    testthat::local_mocked_bindings(
      draws_bayesplot_fun = function(nm, what) function(y, yrep, lw, ...) lw,
      .package = "frmtmb.sample")
    q(pp_check(d, type = "loo_pit_overlay"))
  })
  res("pp_check lw == psis(ds) log weights",
      identical(lw, stats::weights(q(psis(d)), log = TRUE)))
  ids <- seq(1, ndraws(d), by = 3)
  lws <- local({
    testthat::local_mocked_bindings(
      draws_bayesplot_fun = function(nm, what) function(y, yrep, lw, ...) lw,
      .package = "frmtmb.sample")
    q(pp_check(d, type = "loo_pit_overlay", draw_ids = ids))
  })
  pss <- stats::weights(q(loo::psis(-log_lik(d, draw_ids = ids))),
                        log = TRUE)
  res("thinned: pp_check lw == psis with r_eff 1",
      identical(lws, pss))
  res("thinned: max |pp_check lw - r_eff-1 lw|", sprintf("%.3g",
      max(abs(lws - pss))))
}
# the refusals the brief names
fit <- frm(bf(y ~ x), family = gaussian(), data = dd)
for (f in c("loo_epred", "loo_linpred", "loo_predict",
            "loo_predictive_interval", "loo_R2", "loo_model_weights")) {
  m <- tryCatch({get(f)(fit); "no error"}, error = function(e) conditionMessage(e))
  res(paste("ML", f), substr(m, 1, 120))
}
m <- tryCatch({q(loo_epred(ds$gaussian, newdata = dd[1:3, ])); "no error"},
              error = function(e) conditionMessage(e))
res("newdata", substr(m, 1, 120))
m <- tryCatch({q(loo_linpred(ds$gaussian, scale = "response")); "no error"},
              error = function(e) conditionMessage(e))
res("scale", substr(m, 1, 120))
m <- tryCatch({x <- q(add_criterion(ds$gaussian, "loo_R2")); dim(x$criteria$loo_R2)},
              error = function(e) conditionMessage(e))
res("add_criterion loo_R2", paste(m, collapse = " x "))

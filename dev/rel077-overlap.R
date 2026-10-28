# Cross-lane probes of the 0.77.0 merge: the pairs of lanes whose changes
# meet only in the merged tree. Each section prints what it measured;
# dev/round-20261028.md quotes the log. Seeds are in the code.
#
#   Rscript dev/rel077-overlap.R <section> [lib]
#   <section>: sr (sigmainn x reslot), sp (sigmainn x priordraws),
#   rp (reslot x priordraws), mp (mointerop x priordraws),
#   ms (mointerop x sigmainn), mr (mointerop x reslot), all
#   [lib]: the library in front (default rellib-r15)
a <- commandArgs(TRUE)
sec <- if (length(a)) a[1] else "all"
lib <- if (length(a) > 1L) a[2] else "C:/Users/adf44/source/r/rellib-r15"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-rel077/dev/stan-cache")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
options(mc.cores = 1)
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")),
    find.package("frmtmb.sample"), format(packageVersion("frmtmb.sample")),
    "\n")
q <- function(expr) suppressWarnings(suppressMessages(expr))
msg <- function(expr) tryCatch({ force(expr); "no error" },
                               error = function(e) conditionMessage(e))
rp <- function(fit, prior) q(frmtmb.sample:::sample_resolve_priors(fit, prior))
ent <- function(es) vapply(es, function(e) {
  paste0(e$comp, "[", paste(e$idx, collapse = ","), "] ", e$dist$kind,
         if (isTRUE(e$dist$flat)) " (flat lkj)" else "",
         if (!is.null(e$innov)) " INNOV" else "")
}, "")
nlp <- function(fit, es, set) {
  p <- fit$frame$par_template
  for (cmp in names(set)) p[[cmp]][seq_along(set[[cmp]])] <- set[[cmp]]
  frmtmb:::neg_log_prior_fn(es)(p)
}
set.seed(1)
d <- data.frame(x = rnorm(60), h = factor(rep(1:5, each = 12)),
                t = rep(1:12, 5))
d$y <- 1 + 0.5 * d$x + rnorm(60)
far <- y ~ x + ar(time = t, gr = h, cov = TRUE)

if (sec %in% c("sr", "all")) {
  cat("\n===== sigmainn x reslot: sigma under ar(cov = TRUE) =====\n")
  fit <- frm(far, data = d, dry_run = "objective")
  for (pr in list("flat", set_prior("", class = "sigma", lb = 0.1),
                  set_prior("", class = "sigma"))) {
    r <- rp(fit, pr)
    sg <- Filter(function(e) identical(e$comp, "betad"), r$ri$entries)
    lab <- if (is.character(pr)) "\"flat\"" else
      paste0("set_prior(\"\", class = \"sigma\", lb = ", pr[[1]]$lb, ")")
    cat(lab, ": sigma entries:", ent(sg), "\n")
    v <- vapply(c(-0.8, 0, 0.9), function(a) nlp(fit, sg,
      list(betad = 0.3, thetaac = a)), 0)
    cat("  -log prior at log sigma 0.3, thetaac -0.8, 0, 0.9:",
        format(v, digits = 17), "\n")
    if (!is.null(r$ri$lower)) cat("  lower:", paste(names(r$ri$lower),
      format(r$ri$lower), collapse = "; "), "\n")
  }
  cat("-- prior_summary() round trip of a prior = \"flat\" run\n")
  ds <- q(frm_sample(fit, prior = "flat", chains = 1, iter = 300,
                     seed = 3, refresh = 0))
  ps <- prior_summary(ds)
  print(as.data.frame(ps)[, c("prior", "class", "lb", "ub")])
  e1 <- rp(fit, "flat")$ri$entries
  for (back in list(ps, as.brmsprior(as.data.frame(ps)))) {
    e2 <- rp(fit, back)$ri$entries
    set.seed(99)
    pts <- replicate(20, list(betad = rnorm(1, 0, 0.7),
                              thetaac = rnorm(1, 0, 0.7),
                              beta = rnorm(2)), simplify = FALSE)
    dv <- vapply(pts, function(p) nlp(fit, e1, p) - nlp(fit, e2, p), 0)
    cat("  passed back (", class(back)[1], "): max |difference| of -log ",
        "prior at 20 points: ", max(abs(dv)), "\n", sep = "")
  }
  cat("-- the fit-route table on the ar model\n")
  p <- default_prior(far, data = d)
  cat("  ", msg(rp(fit, p)), "\n")
  cat("  coerced: ", msg(rp(fit, as.brmsprior(as.data.frame(p)))), "\n")
  ps2 <- default_prior(far, data = d, route = "sample")
  cat("  sample route coerced:",
      paste(ent(rp(fit, as.brmsprior(as.data.frame(ps2)))$ri$entries),
            collapse = " | "), "\n")
  cat("  sample route as is:  ",
      paste(ent(rp(fit, ps2)$ri$entries), collapse = " | "), "\n")
}

if (sec %in% c("sp", "all")) {
  cat("\n===== sigmainn x priordraws: \"only\" sigma under ar(cov = TRUE) =====\n")
  pr <- set_prior("normal(0, 1)", class = "b") +
    set_prior("normal(1, 2)", class = "Intercept") +
    set_prior("exponential(1)", class = "sigma") +
    set_prior("normal(0, 0.5)", class = "ar")
  do <- q(frm_sample(bf(far), data = d, family = gaussian(), prior = pr,
                     sample_prior = "only", chains = 4, iter = 3000,
                     warmup = 1000, seed = 9, refresh = 0))
  da <- posterior::as_draws_array(do)
  s <- posterior::subset_draws(da, "sigma")
  cat("sigma: mean", mean(s), "mcse", posterior::mcse_mean(s),
      "z against exponential(1)", (mean(s) - 1) / posterior::mcse_mean(s),
      "\n")
  qs <- c(0.05, 0.5, 0.95)
  for (k in qs) {
    cat("  quantile", k, ":", stats::quantile(s, k), "exact",
        stats::qexp(k), "z", (stats::quantile(s, k) - stats::qexp(k)) /
          posterior::mcse_quantile(s, k), "\n")
  }
  ys <- q(frm_sample(bf(far), data = d, family = gaussian(), prior = pr,
                     sample_prior = "yes", chains = 2, iter = 1000,
                     seed = 9, refresh = 0))
  pd <- prior_draws(ys)
  cat("\"yes\" prior columns:", names(pd), "\n")
  h <- q(hypothesis(ys, "sigma = 1", class = ""))
  print(h$hypothesis)
  cat("  the numerator reads the draws' sigma (innovation sd): posterior",
      "mean", mean(as.matrix(ys)[, "sigma"]), "; the denominator",
      "prior_sigma: mean", mean(pd$sigma), "\n")
}

if (sec %in% c("rp", "all")) {
  cat("\n===== reslot x priordraws: prior = \"flat\" =====\n")
  set.seed(2)
  dg <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
  dg$y <- 1 + 0.5 * dg$x + rnorm(6)[dg$g] + rnorm(60)
  ff <- bf(y ~ x + (1 + x | g))
  ds <- q(frm_sample(ff, data = dg, family = gaussian(), prior = "flat",
                     sample_prior = "yes", chains = 1, iter = 400,
                     seed = 5, refresh = 0))
  cat("effective prior kinds:", paste(ent(rp(frm(ff, data = dg,
      dry_run = "objective"), "flat")$ri$entries), collapse = " | "), "\n")
  cat("prior_draws() columns:", ncol(prior_draws(ds)), names(prior_draws(ds)),
      "\n")
  h <- q(hypothesis(ds, "x = 0"))
  cat("hypothesis(x = 0) Evid.Ratio:", h$hypothesis$Evid.Ratio, "\n")
  h2 <- q(hypothesis(ds, "sd_g__Intercept = 0.5", class = ""))
  cat("hypothesis(sd = 0.5) Evid.Ratio:", h2$hypothesis$Evid.Ratio, "\n")
  cat("\"only\":", msg(q(frm_sample(ff, data = dg, family = gaussian(),
      prior = "flat", sample_prior = "only", chains = 1, iter = 100,
      refresh = 0))), "\n")
  cat("-- a mixture's weights under prior = \"flat\" keep dirichlet(1)\n")
  set.seed(4)
  dm <- data.frame(y = c(rnorm(40, -2), rnorm(40, 2)))
  mpr <- set_prior("normal(-2, 1)", class = "Intercept", dpar = "mu1") +
    set_prior("normal(2, 1)", class = "Intercept", dpar = "mu2")
  dmx <- q(frm_sample(bf(y ~ 1), data = dm,
                      family = mixture(gaussian(), gaussian()),
                      prior = "flat", sample_prior = "yes", chains = 1,
                      iter = 400, seed = 6, refresh = 0))
  cat("  mixture prior columns:", names(prior_draws(dmx)), "\n")
}

if (sec %in% c("mp", "all")) {
  cat("\n===== mointerop x priordraws: a two-value mo() under \"yes\" =====\n")
  set.seed(3)
  dm <- data.frame(m2 = rep(1:2, 30), x = rnorm(60))
  dm$y <- 0.5 * dm$m2 + rnorm(60)
  ds <- q(frm_sample(bf(y ~ mo(m2)), data = dm, family = gaussian(),
                     prior = set_prior("normal(0, 1)", class = "b"),
                     sample_prior = "yes", chains = 2, iter = 400, seed = 7,
                     refresh = 0))
  pd <- prior_draws(ds)
  cat("columns:", names(pd), "\n")
  w <- pd[[grep("^simo", names(pd))]]
  cat("prior simo weight: draws", length(w), "all equal to 1:",
      all(w == 1), "\n")
  cat("posterior simo weight all 1:",
      all(as.matrix(ds)[, grep("^simo", variables(ds))] == 1), "\n")
}

if (sec %in% c("ms", "all")) {
  cat("\n===== mointerop x sigmainn: marginaleffects on ar(cov = TRUE) =====\n")
  fit <- q(frm(bf(far, sigma ~ x), data = d))
  nd <- data.frame(x = c(-1, 0, 1), t = 1, h = factor(1, levels = 1:5))
  me <- q(marginaleffects::predictions(fit, newdata = nd, dpar = "sigma"))
  ft <- fitted(fit, newdata = nd, dpar = "sigma")
  cat("predictions(dpar = sigma):", me$estimate, "\n")
  cat("fitted(dpar = sigma):     ", ft[, "Estimate"], "\n")
  cat("SE:", me$std.error, "against fitted Est.Error", ft[, "Est.Error"], "\n")
  fit0 <- q(frm(far, data = d))
  me0 <- q(marginaleffects::predictions(fit0, newdata = nd))
  f0 <- frm_linpred(fit0, newdata = nd, se.fit = TRUE)
  cat("mean SE:", me0$std.error, "against frm_linpred", f0$se.fit, "\n")
  cat("sigma(fit0)", sigma(fit0), "; predictions(dpar = sigma)",
      q(marginaleffects::predictions(fit0, newdata = nd[1, ],
                                     dpar = "sigma"))$estimate, "\n")
  old <- fit0
  old$frame$autocor[[1L]]$sigma_is <- NULL
  cat("saved 0.76.0 fit, mean:",
      msg(q(marginaleffects::predictions(old, newdata = nd))), "\n")
  cat("saved 0.76.0 fit, dpar = sigma:",
      substr(msg(q(marginaleffects::predictions(old, newdata = nd,
                                                dpar = "sigma"))), 1, 160),
      "\n")
  cat("saved 0.76.0 fit, avg_slopes:",
      substr(msg(q(marginaleffects::avg_slopes(old, variables = "x"))), 1,
             160), "\n")
}

if (sec %in% c("mr", "all")) {
  cat("\n===== mointerop x reslot: one- and two-value mo() priors =====\n")
  set.seed(3)
  dm <- data.frame(m2 = rep(1:2, 30), m1 = 1L, x = rnorm(60))
  dm$y <- 0.5 * dm$m2 + rnorm(60)
  f2 <- frm(bf(y ~ mo(m2)), data = dm, dry_run = "objective")
  cat("two values, prior = \"flat\":", ent(rp(f2, "flat")$ri$entries), "\n")
  cat("two values, default:", ent(rp(f2, NULL)$ri$entries), "\n")
  cat("two values, dirichlet(2) on simo:",
      msg(r <- rp(f2, set_prior("dirichlet(2)", class = "simo",
                                coef = "mom21"))), ent(r$ri$entries), "\n")
  cat("two values, (flat) on simo:",
      msg(r <- rp(f2, set_prior("", class = "simo", coef = "mom21"))),
      ent(r$ri$entries), "\n")
  fm2 <- q(frm(bf(y ~ mo(m2)), data = dm,
               prior = set_prior("dirichlet(2)", class = "simo",
                                 coef = "mom21")))
  fm0 <- q(frm(bf(y ~ mo(m2)), data = dm))
  cat("frm() with dirichlet(2) on the one weight: logLik", logLik(fm2),
      "against none", logLik(fm0), "\n")
  ds <- q(frm_sample(bf(y ~ mo(m2)), data = dm, family = gaussian(),
                     prior = "flat", chains = 1, iter = 300, seed = 2,
                     refresh = 0))
  cat("frm_sample(prior = \"flat\") two values: variables",
      variables(ds), "\n")
  cat("one value, frm_sample(prior = \"flat\"):",
      substr(msg(q(frm_sample(bf(y ~ x + mo(m1)), data = dm,
                              family = gaussian(), prior = "flat",
                              chains = 1, iter = 100, refresh = 0))), 1, 200),
      "\n")
  set.seed(5)
  d4 <- data.frame(k = sample(1:4, 80, TRUE))
  d4$y <- d4$k + rnorm(80)
  f4 <- frm(bf(y ~ mo(k)), data = d4, dry_run = "objective")
  cat("four values, prior = \"flat\":", ent(rp(f4, "flat")$ri$entries), "\n")
  cat("four values, class theta on the simplex:",
      substr(msg(rp(f4, set_prior("normal(0, 1)", class = "theta"))), 1,
             160), "\n")
}
cat("\nOVERLAP DONE\n")

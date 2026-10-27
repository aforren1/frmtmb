# Reviewer: every route a prior takes to an autocorrelation parameter,
# before (rellib-r13) and after (the lane library). No sampling: the
# resolved entries and bounds of frm()'s and frm_sample()'s priors.
# Usage: Rscript dev/review-acprior-probe.R before|after
# Log: dev/review-acprior-log/probe-<arm>.txt
arm <- commandArgs(TRUE)[1]
lib <- c(if (identical(arm, "after")) "C:/Users/adf44/source/r/wt-acprior-lib",
         "C:/Users/adf44/source/r/rellib-r13",
         "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(lib)
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "; sample from",
    dirname(find.package("frmtmb.sample")), "\n")
set.seed(5)
ng <- 10; k <- 5
d <- data.frame(t = rep(seq_len(k), ng), g = rep(seq_len(ng), each = k),
                x = rnorm(ng * k))
d$y1 <- 1 + 0.5 * d$x + rnorm(ng)[d$g] * 0.5 + rnorm(ng * k)
d$y2 <- 1 - 0.3 * d$x + rnorm(ng)[d$g] * 0.5 + rnorm(ng * k)
d$y <- d$y1
msg <- function(e) cat("  ERROR:", conditionMessage(e), "\n")
show_ri <- function(lab, fit, pr) {
  cat("\n---", lab, "\n")
  tryCatch({
    ri <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, pr))$ri
    for (en in ri$entries) {
      if (!identical(en$comp, "thetaac")) next
      cat("  sample entry", en$comp, paste(en$idx, collapse = ","),
          en$dist$kind,
          if (identical(en$dist$kind, "trans")) en$dist$inner$kind, "\n")
    }
    lo <- ri$lower[grepl("thetaac", names(ri$lower))]
    up <- ri$upper[grepl("thetaac", names(ri$upper))]
    if (length(lo)) cat("  lower:", paste(names(lo), signif(lo, 6),
                                          collapse = " "), "\n")
    if (length(up)) cat("  upper:", paste(names(up), signif(up, 6),
                                          collapse = " "), "\n")
  }, error = msg)
}
show_ml <- function(lab, fo, pr, data = d, data2 = NULL) {
  cat("\n--- ML", lab, "\n")
  tryCatch({
    fit <- frm(fo, data = data, data2 = data2, prior = pr, dry_run = "objective")
    ri <- frmtmb:::resolve_prior_input(fit, pr)
    for (en in ri$entries) {
      cat("  entry", en$comp, paste(en$idx, collapse = ","), en$dist$kind,
          if (identical(en$dist$kind, "trans")) en$dist$inner$kind, "\n")
    }
    if (length(ri$lower)) cat("  lower:", paste(names(ri$lower),
                                                signif(ri$lower, 6)), "\n")
    if (length(ri$upper)) cat("  upper:", paste(names(ri$upper),
                                                signif(ri$upper, 6)), "\n")
  }, error = msg)
}
tab <- function(lab, fo, route = "fit", data = d) {
  cat("\n=== table", lab, "route", route, "\n")
  tryCatch({
    tb <- default_prior(fo, data = data, route = route)
    tb <- as.data.frame(tb)
    keep <- intersect(c("prior", "class", "coef", "group", "resp", "dpar",
                        "lb", "ub"), names(tb))
    print(tb[!tb$class %in% c("b", "Intercept"), keep], row.names = FALSE)
  }, error = msg)
}

fmv2 <- bf(y1 ~ x + ar(t, g, cov = TRUE)) + bf(y2 ~ x + ar(t, g, cov = TRUE)) +
  set_rescor(FALSE)
fmvc <- bf(y1 ~ x + ar(t, g, cov = TRUE)) + bf(y2 ~ x + cosy(t, g)) +
  set_rescor(FALSE)
fmv3 <- bf(y1 ~ x + cosy(t, g)) + bf(y2 ~ x + cosy(t, g)) + set_rescor(FALSE)
tab("mv ar+cosy", fmvc)
tab("mv ar+cosy", fmvc, "sample")
fit2 <- frm(fmv2, data = d, dry_run = "objective")
fitc <- frm(fmvc, data = d, dry_run = "objective")
fit3 <- frm(fmv3, data = d, dry_run = "objective")
show_ri("mv ar both: default", fit2, NULL)
show_ri("mv ar both: ar no resp, lb = 0.2 (brms refuses)", fit2,
        set_prior("", class = "ar", lb = 0.2))
show_ri("mv ar both: ar resp y1, lb = 0.2", fit2,
        set_prior("", class = "ar", resp = "y1", lb = 0.2))
show_ri("mv cosy both: default", fit3, NULL)
show_ri("mv cosy both: cosy no resp normal(0, 0.1) (brms refuses)", fit3,
        set_prior("normal(0, 0.1)", class = "cosy"))
show_ri("mv cosy both: cosy resp y2 normal(0, 0.1)", fit3,
        set_prior("normal(0, 0.1)", class = "cosy", resp = "y2"))
show_ml("mv ar both: ar no resp normal (brms refuses)", fmv2,
        set_prior("normal(0.5, 0.2)", class = "ar"))
for (j in 1:2) {
  show_ml(paste0("mv ar+cosy: theta thetaac_", j), fmvc,
          set_prior("normal(0, 1)", class = "theta",
                    coef = paste0("thetaac_", j)))
}
show_ml("ar p=2 cov=TRUE: theta thetaac_1 lb", y ~ x + ar(t, g, p = 2, cov = TRUE),
        set_prior("", class = "theta", coef = "thetaac_1", lb = 0))
show_ml("ar p=2 cov=TRUE: ar lb", y ~ x + ar(t, g, p = 2, cov = TRUE),
        set_prior("", class = "ar", lb = 0))
show_ml("ar p=2 cov=TRUE: ar normal", y ~ x + ar(t, g, p = 2, cov = TRUE),
        set_prior("normal(0, 0.5)", class = "ar"))
show_ml("ar p=2 cov=FALSE: ar lb", y ~ x + ar(t, g, p = 2),
        set_prior("", class = "ar", lb = 0))
show_ml("ar p=2 cov=FALSE: theta thetaac_2", y ~ x + ar(t, g, p = 2),
        set_prior("normal(0, 1)", class = "theta", coef = "thetaac_2"))
show_ml("arma(p=2,q=1): ma lb", y ~ x + arma(t, g, p = 2, q = 1),
        set_prior("", class = "ma", lb = -0.5))
show_ml("ar() no gr: theta", y ~ x + ar(tt),
        set_prior("normal(0, 1)", class = "theta", coef = "thetaac_1"),
        data = transform(d, tt = seq_len(nrow(d))))
show_ml("cosy(gr=): theta", y ~ x + cosy(t, g),
        set_prior("normal(0, 1)", class = "theta", coef = "thetaac_1"))
show_ml("unstr: theta thetaac_3", y ~ x + unstr(t, g),
        set_prior("normal(0, 1)", class = "theta", coef = "thetaac_3"))
show_ml("unstr: cortime lkj(2) + theta thetaac_1",
        y ~ x + unstr(t, g),
        c(set_prior("lkj(2)", class = "cortime"),
          set_prior("normal(0, 1)", class = "theta", coef = "thetaac_1")))
show_ml("unstr + (1|g): class theta, no coef", y ~ x + unstr(t, g) + (1 | g),
        set_prior("normal(0, 1)", class = "theta"))
show_ml("ar cov only: class theta, no coef", y ~ x + ar(t, g, cov = TRUE),
        set_prior("normal(0, 1)", class = "theta"))
show_ml("ar cov + (1|g): class theta, no coef", y ~ x + ar(t, g, cov = TRUE) +
          (1 | g), set_prior("normal(0, 1)", class = "theta"))
show_ml("fcor: theta", y ~ x + fcor(M),
        set_prior("normal(0, 1)", class = "theta", coef = "thetaac_1"),
        data2 = list(M = diag(nrow(d))))
show_ml("fcor + ar cov in mv? fcor alone + (1|g): class theta", y ~ x + fcor(M) +
          (1 | g), set_prior("normal(0, 1)", class = "theta"),
        data2 = list(M = diag(nrow(d))))
show_ml("prior_string ar", y ~ x + ar(t, g, cov = TRUE),
        prior_string("normal(0.5, 0.2)", class = "ar"))
show_ml("ar cov: lb = -1 (table end)", y ~ x + ar(t, g, cov = TRUE),
        set_prior("", class = "ar", lb = -1))
show_ml("ar cov: ub = -1 (wrong end)", y ~ x + ar(t, g, cov = TRUE),
        set_prior("", class = "ar", ub = -1))
show_ml("ar cov: lb = 1 (wrong end)", y ~ x + ar(t, g, cov = TRUE),
        set_prior("", class = "ar", lb = 1))
show_ml("ar cov: lb = -1 + 1e-16", y ~ x + ar(t, g, cov = TRUE),
        set_prior("", class = "ar", lb = -1 + 1e-16))
show_ml("cosy: lb = -0.25 (window end, d = 5)", y ~ x + cosy(t, g),
        set_prior("", class = "cosy", lb = -0.25))
show_ml("cosy: lb = -1/4 written as -0.2500000001", y ~ x + cosy(t, g),
        set_prior("", class = "cosy", lb = -0.2500000001))
show_ml("cosy: ub = 1.5 (brms takes)", y ~ x + cosy(t, g),
        set_prior("", class = "cosy", ub = 1.5))

# the sampling table passed back as prior =
for (fo in list(y ~ x + cosy(t, g), y ~ x + ar(t, g, cov = TRUE),
                y ~ x + ar(t, g, p = 2, cov = TRUE))) {
  cat("\n--- sample table round trip:", deparse(fo), "\n")
  tryCatch({
    tb <- default_prior(fo, data = d, route = "sample")
    fit <- frm(fo, data = d, dry_run = "objective")
    r0 <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, NULL))$ri
    r1 <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, tb))$ri
    p <- fit$frame$par_template
    f0 <- frmtmb:::neg_log_prior_fn(r0$entries)
    f1 <- frmtmb:::neg_log_prior_fn(r1$entries)
    vals <- vapply(list(-1, 0.3, 2), function(th) {
      p$thetaac[] <- th
      c(f0(p), f1(p))
    }, c(0, 0))
    print(vals)
    cat("  lower default:", signif(r0$lower, 6), " table:", signif(r1$lower, 6),
        "\n")
  }, error = msg)
}

# the same double slot for random effects and rescor
cat("\n=== RE / rescor double slots\n")
show_ml("RE: cor lkj(2) + theta theta_3", y ~ x + (1 + x | g),
        c(set_prior("lkj(2)", class = "cor"),
          set_prior("normal(0, 1)", class = "theta", coef = "theta_3")))
show_ml("RE: theta theta_3 + cor lkj(2)", y ~ x + (1 + x | g),
        c(set_prior("normal(0, 1)", class = "theta", coef = "theta_3"),
          set_prior("lkj(2)", class = "cor")))
show_ml("RE: sd normal + theta theta_1", y ~ x + (1 | g),
        c(set_prior("normal(0, 1)", class = "sd"),
          set_prior("normal(0, 1)", class = "theta", coef = "theta_1")))
show_ml("RE: theta theta_1 + sd normal", y ~ x + (1 | g),
        c(set_prior("normal(0, 1)", class = "theta", coef = "theta_1"),
          set_prior("normal(0, 1)", class = "sd")))
fmr <- bf(y1 ~ x) + bf(y2 ~ x) + set_rescor(TRUE)
show_ml("rescor lkj(2) + theta thetar_1", fmr,
        c(set_prior("lkj(2)", class = "rescor"),
          set_prior("normal(0, 1)", class = "theta", coef = "thetar_1")))
show_ml("theta thetar_1 + rescor lkj(2)", fmr,
        c(set_prior("normal(0, 1)", class = "theta", coef = "thetar_1"),
          set_prior("lkj(2)", class = "rescor")))

# frm_sample's bounded natural defaults under a user prior
cat("\n=== bounded defaults under a user prior (frm_sample)\n")
set.seed(9)
dd <- data.frame(x = rnorm(60))
dd$y <- 1 + dd$x + rnorm(60)
dd$yp <- 0.3 + exp(rnorm(60, 0, 0.5))
showb <- function(lab, fit, pr) {
  cat("\n---", lab, "\n")
  tryCatch({
    ri <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, pr))$ri
    cat("  lower:", paste(names(ri$lower), signif(ri$lower, 6)), "\n")
    cat("  upper:", paste(names(ri$upper), signif(ri$upper, 6)), "\n")
  }, error = msg)
}
fs <- frm(y ~ x, data = dd, family = brmsfamily("gaussian", link_sigma = "identity"),
          dry_run = "objective")
showb("sigma identity: default", fs, NULL)
showb("sigma identity: normal(1, 1)", fs,
      set_prior("normal(1, 1)", class = "sigma"))
ft <- frm(y ~ x, data = dd, family = brmsfamily("student", link_nu = "identity"),
          dry_run = "objective")
showb("nu identity: default", ft, NULL)
showb("nu identity: normal(5, 2)", ft, set_prior("normal(5, 2)", class = "nu"))
fn <- frm(yp ~ x, data = dd, family = brmsfamily("shifted_lognormal"),
          dry_run = "objective")
showb("ndt identity: default", fn, NULL)
showb("ndt identity: normal(0.2, 0.1)", fn,
      set_prior("normal(0.2, 0.1)", class = "ndt"))

# Reviewer probe 4: (i) the named-list prior as a second slot; (ii) round
# trips of a run's printed prior table passed back as prior =, compared
# on the prior-augmented objective's value and gradient (one tape per
# arm, no draws compared across tmbstan runs); (iii) frm() under the
# sampling table on ar(cov = TRUE).
# Usage: Rscript dev/review-reslot-p4.R after. Seeds: data 1, points 7,
# sampler 1.
source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
suppressPackageStartupMessages(library(frmtmb.sample))
dd <- re_data(1L)
q <- function(expr) suppressWarnings(suppressMessages(expr))
rp <- function(fit, prior) q(frmtmb.sample:::sample_resolve_priors(fit, prior))
show_e <- function(ri) {
  for (e in ri$entries) {
    cat(sprintf("  entry %-8s %-10s %-8s %-9s flat=%s\n", e$comp,
                paste(e$idx, collapse = ","), e$dist$kind, e$scale,
                isTRUE(e$dist$flat)))
  }
}
cat("== (i) named list beside the class stack, (1 + x + z | g)\n")
fit3 <- q(frm(bf(y ~ x + (1 + x + z | g)), data = dd, dry_run = "objective"))
r <- tryCatch(rp(fit3, list(theta_4 = prior_normal(0, 1))),
              error = function(e) e)
if (inherits(r, "error")) cat("  ERROR:", conditionMessage(r), "\n") else
  show_e(r$ri)
r <- tryCatch(rp(fit3, list(theta_2 = prior_normal(0, 1))),
              error = function(e) e)
if (inherits(r, "error")) cat("  ERROR:", conditionMessage(r), "\n") else
  show_e(r$ri)
r <- tryCatch(q(frm(bf(y ~ x + (1 + x + z | g)), data = dd,
                    prior = list(theta_5 = prior_normal(0, 1)))),
              error = function(e) e)
cat("  frm(prior = list(theta_5 = )):",
    if (inherits(r, "error")) conditionMessage(r) else "fits", "\n")

cat("\n== (ii) round trips of prior_summary() of a run\n")
cmp_obj <- function(fit, e0, e1, label) {
  o0 <- frmtmb.sample:::prior_augmented_obj(fit, e0)
  o1 <- frmtmb.sample:::prior_augmented_obj(fit, e1)
  set.seed(7)
  df <- dg <- 0
  for (j in 1:50) {
    p <- o0$par + stats::rnorm(length(o0$par), 0, 0.5)
    f0 <- o0$fn(p)
    df <- max(df, abs(o1$fn(p) - f0) / max(1, abs(f0)))
    g0 <- o0$gr(p)
    dg <- max(dg, max(abs(o1$gr(p) - g0)) / max(1, abs(g0)))
  }
  cat(sprintf("  %-40s rel fn diff %.3g, rel grad diff %.3g\n", label, df, dg))
}
fit2 <- q(frm(bf(y ~ x + (1 + x + z | g)), data = dd))
for (pr in list("flat", NULL)) {
  lab <- if (is.null(pr)) "default" else pr
  ds <- q(frm_sample(fit2, prior = pr, chains = 1, iter = 300, warmup = 150,
                     seed = 1, refresh = 0))
  ps <- prior_summary(ds)
  cat("-- prior_summary() of the", lab, "run:\n")
  print(as.data.frame(ps)[, intersect(c("prior", "class", "coef", "group",
                                        "lb", "ub", "source"),
                                      names(as.data.frame(ps)))],
        row.names = FALSE)
  e0 <- rp(fit2, pr)$ri$entries
  for (form in c("object", "data.frame")) {
    back <- if (form == "object") ps else as.data.frame(ps)
    r1 <- tryCatch(rp(fit2, back), error = function(e) e)
    if (inherits(r1, "error")) {
      cat("  ", lab, form, "passed back: ERROR", conditionMessage(r1), "\n")
      next
    }
    cmp_obj(fit2, e0, r1$ri$entries, paste(lab, "run,", form, "passed back"))
  }
}

cat("\n== (iii) frm() under the sampling table, ar(cov = TRUE)\n")
set.seed(13)
G <- 6L; n <- 8L
da <- data.frame(g = factor(rep(seq_len(G), each = n)), t = rep(seq_len(n), G))
da$y <- 1 + unlist(lapply(seq_len(G), function(j)
  as.numeric(stats::arima.sim(list(ar = 0.6), n))))
fa <- bf(y ~ 1 + ar(t, g, cov = TRUE))
tab <- default_prior(fa, data = da, route = "sample")
print(as.data.frame(tab)[, c("prior", "class", "lb", "ub")], row.names = FALSE)
f_ml <- q(frm(fa, data = da))
f_tab <- q(frm(fa, data = da, prior = tab))
fitd <- q(frm(fa, data = da))
ps <- prior_summary(q(frm_sample(fitd, chains = 1, iter = 200, warmup = 100,
                                 seed = 1, refresh = 0)))
f_ps <- q(frm(fa, data = da, prior = ps))
f_uni <- q(frm(fa, data = da, prior = set_prior("uniform(-1, 1)", class = "ar") +
                 set_prior("student_t(3, 0, 2.5)", class = "sigma") +
                 set_prior(tab$prior[tab$class == "Intercept"],
                           class = "Intercept")))
est <- function(f) unlist(f$estimates[c("beta", "betad", "thetaac")])
cat("  ML               :", format(est(f_ml), digits = 7), "\n")
cat("  prior = tab      :", format(est(f_tab), digits = 7), "\n")
cat("  prior = summary  :", format(est(f_ps), digits = 7), "\n")
cat("  uniform(-1,1) ar :", format(est(f_uni), digits = 7), "\n")

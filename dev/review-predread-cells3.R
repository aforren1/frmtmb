# Item 1: the cells whose Fixed_param route has no translation (mixture x cs(), mv x
# cs(), categorical with per-category predictors): layout from a short brms fit, values
# from the model's own algebra read off each package's parameters. Data seed 7311.
# Also the readers of a layered posterior_linpred().
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7311)
n <- 400
d <- data.frame(x = rnorm(n), z = rnorm(n))
eta <- 0.7 * d$x
d$ord <- as.integer(cut(eta + ifelse(d$z > 0, 1.2, -1.2) * rlogis(n), c(-Inf, -0.5, 0.8, 2, Inf), labels = FALSE))
d$cat <- as.integer(cut(eta + 0.5 * d$z + rlogis(n), c(-Inf, 0, 1.5, Inf), labels = FALSE))
d$y <- eta + rnorm(n)
comp <- rbinom(n, 1, 0.5)
d$mord <- as.integer(cut(ifelse(comp == 1, eta + 1, -eta) + rlogis(n), c(-Inf, -0.5, 0.8, 2, Inf), labels = FALSE))
nd <- data.frame(x = c(-1, 0.3, 1.2, 0.5), z = c(0.7, -1.1, 1.4, 0.2))
fmt <- function(v) paste(signif(as.numeric(v), 5), collapse = " ")
br <- function(ftxt, famtxt) {
  fam <- if (is.null(famtxt)) NULL else rev_ev(famtxt, "brms")
  quiet(brms::brm(rev_ev(ftxt, "brms"), family = fam, data = d, chains = 1, iter = 500, refresh = 0,
                  seed = 11, silent = 2, init = 0))
}
msg <- function(e) {
  r <- tryCatch(quiet(e), error = function(err) paste("ERROR:", substr(conditionMessage(err), 1, 170)))
  if (is.character(r)) r else fmt(r)
}
cat("\n=== categorical with per-category predictors\n")
ft <- "bf(cat ~ 1, mu2 ~ x + z, mu3 ~ z)"
fa <- "categorical()"
fit <- quiet(frm(rev_ev(ft, "frmtmb"), family = rev_ev(fa, "frmtmb"), data = d))
bfit <- tryCatch(br(ft, fa), error = function(e) {cat("brms ERR:", conditionMessage(e), "
"); NULL})
b <- if (!is.null(bfit)) as.matrix(bfit)
pl <- if (!is.null(bfit)) posterior_linpred(bfit, newdata = nd, re_formula = NA)
cat("brms dim", dim(pl), "\n")
lf <- frm_linpred(fit, newdata = nd)
cat("frm dim", dim(lf), "\n")
est <- fit$estimates$beta
print(round(est, 3))
if (!is.null(pl)) cat("brms algebra rel:", rev_E$eq_rel(pl[1, , ], cbind(b[1, "b_mu2_Intercept"] + b[1, "b_mu2_x"] * nd$x + b[1, "b_mu2_z"] * nd$z, b[1, "b_mu3_Intercept"] + b[1, "b_mu3_z"] * nd$z)), "\n")
cat("frm algebra rel:", rev_E$eq_rel(lf, cbind(est[1] + est[2] * nd$x + est[3] * nd$z, est[4] + est[5] * nd$z)), "\n")
cat("frm_curve est rel:", rev_E$eq_rel(frm_curve(fit, newdata = nd, seed = 1)$.estimate, lf), "\n")
li <- frm_linpred(fit, newdata = nd, se.fit = TRUE)
cat("se.fit dims", dim(li$se.fit), " vals ", fmt(li$se.fit), "\n")
V <- vcov(fit, full = FALSE)
if (nrow(V) == length(est)) {
  cat("se by hand mu3:", fmt(sqrt(diag(cbind(1, nd$z) %*% V[4:5, 4:5] %*% t(cbind(1, nd$z))))), "\n")
}
cat("\n=== readers of a layered posterior_linpred()\n")
sf <- quiet(frm_sample(bf(ord ~ x + cs(z)), family = sratio(), data = d, chains = 1, iter = 400, seed = 3, refresh = 0))
bf2 <- br("bf(ord ~ x + cs(z))", "sratio()")
for (lab in c("loo_linpred", "posterior_linpred", "CE dpar mu", "CE linpred")) {
  fb <- switch(lab,
    loo_linpred = function(o) dim(loo_linpred(o)),
    posterior_linpred = function(o) dim(posterior_linpred(o, ndraws = 5)),
    "CE dpar mu" = function(o) brms::conditional_effects(o, "x", dpar = "mu", resolution = 3)[[1]]$estimate__,
    "CE linpred" = function(o) brms::conditional_effects(o, "x", method = "posterior_linpred", resolution = 3)[[1]]$estimate__)
  ff <- switch(lab,
    loo_linpred = function(o) dim(frmtmb.sample::loo_linpred(o)),
    posterior_linpred = function(o) dim(frmtmb.sample::posterior_linpred(o, ndraws = 5)),
    "CE dpar mu" = function(o) conditional_effects(o, "x", dpar = "mu", resolution = 3)[[1]]$estimate__,
    "CE linpred" = function(o) conditional_effects(o, "x", method = "posterior_linpred", resolution = 3)[[1]]$estimate__)
  cat(lab, ": brms", msg(fb(bf2)), " | draws", msg(ff(sf)), "\n")
}
fx <- rev_ev("bf(ord ~ x + cs(z))", "frmtmb")
fm <- quiet(frm(fx, family = sratio(), data = d))
cat("core CE dpar mu:", msg(conditional_effects(fm, "x", dpar = "mu", resolution = 3)[[1]]$estimate__), "\n")
cat("core CE method linpred:", msg(conditional_effects(fm, "x", method = "posterior_linpred", resolution = 3)[[1]]$estimate__), "\n")
cat("predictive_error:", msg(dim(frmtmb.sample::predictive_error(sf, ndraws = 5))), "\n")
cat("hypothesis-free summaries, residuals:", msg(dim(frmtmb.sample:::residuals.frmtmb_draws(sf, ndraws = 5, summary = FALSE))), "\n")
cat("as.numeric in curve_predict_eta is length ", length(frmtmb:::curve_predict_eta(fm, nd, NULL, NULL, NA, FALSE)), "\n")
cat("insight link:", msg(as.numeric(insight::get_predicted(fm, data = nd, predict = "link"))), "\n")
cat("insight default:", msg(as.numeric(insight::get_predicted(fm, data = nd))), "\n")
cat("get_predict default:", msg(marginaleffects::get_predict(fm, newdata = nd)$estimate), "\n")
cat("predictions link se:", msg(marginaleffects::predictions(fm, newdata = nd, type = "link")$std.error), "\n")
cat("frm_linpred se:", msg(frm_linpred(fm, newdata = nd, se.fit = TRUE)$se.fit), "\n")
cat("avg_slopes link:", msg(marginaleffects::avg_slopes(fm, variables = "z", newdata = nd, type = "link")$estimate), " est bcs ", msg(fm$estimates[[grep("cs", names(fm$estimates), value = TRUE)[1L]]]), "\n")
cat("\n=== feature and curve with resp on a univariate fit\n")
fu <- quiet(frm(bf(y ~ x), family = gaussian(), data = d))
gx <- data.frame(x = seq(-3, 3, by = 0.5))
for (r in list(character(0), "y", c("y", "y"), "nosuch", NA_character_, 1L)) {
  cat(paste(deparse(r), collapse = ""), ": feature", msg(frm_curve_feature(fu, "x", type = "crossing", at = 0, newdata = gx, resp = r)$.estimate), "\n")
}
cat("\n=== plot of mv conditional effects, points\n")
d$y2 <- -d$x + rnorm(n)
mv2 <- quiet(frm(mvbf(bf(y ~ x), bf(y2 ~ x + z), rescor = FALSE), family = gaussian(), data = d))
ce <- conditional_effects(mv2, "x", resolution = 4)
cat("attr points present:", !is.null(attr(ce[[1]], "points")), " mean resp__ by panel:", round(mean(attr(ce[[1]], "points")$resp__), 4), round(mean(attr(ce[[2]], "points")$resp__), 4), " data means ", round(mean(d$y), 4), round(mean(d$y2), 4), "\n")
grDevices::pdf(NULL)
pp <- plot(ce, points = TRUE, ask = FALSE, plot = FALSE)
cat("class:", class(pp), " class 1:", class(pp[[1]]), "\n")
cat("print each panel:", msg({lapply(pp, print); 1}), "\n")
cat("plot(ce) direct:", msg({plot(ce, ask = FALSE, points = TRUE); 1}), "\n")
md <- rev_draws_at(mv2, 6L)
ced <- conditional_effects(md, "x", resolution = 4)
cat("draws plot:", msg({plot(ced, ask = FALSE, points = TRUE); 1}), "\n")
cat("draws points mean:", round(mean(attr(ced[[2]], "points")$resp__), 4), "\n")
cat("predict method mv points:", msg({plot(conditional_effects(mv2, "x", resolution = 4, method = "predict"), ask = FALSE, points = TRUE); 1}), "\n")
cat("select_points on mv:", msg(nrow(attr(conditional_effects(mv2, "x", resolution = 4, select_points = 0.5)[[2]], "points"))), "\n")
grDevices::dev.off()
cat("\n=== pp_mixture multivariate against the univariate fit of that response\n")
d$ym <- ifelse(comp == 1, 3 + d$x, -3 - d$x) + rnorm(n, 0, 0.7)
d$ym2 <- ifelse(comp == 1, 1 + d$z, -2 - d$z) + rnorm(n, 0, 0.7)
fm2 <- quiet(frm(frmtmb::mvbf(bf(ym ~ x, family = mixture(gaussian(), gaussian())),
                              bf(ym2 ~ z, family = mixture(gaussian(), gaussian())), rescor = FALSE), data = d))
fu2 <- quiet(frm(bf(ym2 ~ z), family = mixture(gaussian(), gaussian()), data = d))
a <- pp_mixture(fm2, resp = "ym2")
b1 <- pp_mixture(fu2)
cat("fit: first rows, mv", round(a[1:3, "Estimate", 1], 4), " uni ", round(b1[1:3, "Estimate", 1], 4), " max abs ", signif(max(abs(a[, "Estimate", ] - b1[, "Estimate", ])), 3), "\n")
dmm <- rev_draws_at(fm2, 3L)
a2 <- frmtmb.sample::pp_mixture(dmm, resp = "ym2")
cat("draws: max abs vs fit", signif(max(abs(a2[, "Estimate", ] - a[, "Estimate", ])), 3), "\n")
cat("brms stop text:", msg(pp_mixture(br("bf(y ~ x)", "gaussian()"))), "\n")

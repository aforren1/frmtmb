# Items 1, 3, 4, 6: readers the other scripts do not reach. Data seed 7901. Each
# block prints before and after when run with REV_LIB=before and without.
#   Rscript dev/review-predread-misc.R
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
library(ggplot2)
before <- identical(Sys.getenv("REV_LIB"), "before")
set.seed(7901)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n), f = factor(sample(c("a", "b", "c"), n, TRUE)),
                g = factor(rep(sprintf("g%02d", 1:30), each = 10)),
                u = factor(rep(sprintf("u%02d", 1:30), 10)))
d$t <- rep(1:10, 30)
d$sx <- 0.4
d$lat <- rnorm(n)
d$xo <- d$lat + rnorm(n, 0, 0.4)
gl <- rnorm(30)[d$g]
d$sg <- 0.4
d$xg <- gl + rnorm(30, 0, 0.4)[d$g]
d$y <- 1 + 1.2 * d$lat - 0.5 * d$z + rnorm(n, 0, 0.5)
d$yg <- 1 + 1.2 * gl - 0.5 * d$z + rnorm(n, 0, 0.5)
d$y2 <- -1 + 0.6 * d$x + 0.5 * d$z + rnorm(n, 0, 0.7)
d$y3 <- 0.5 + 0.8 * d$x - 0.5 * d$z + rnorm(n, 0, 0.7)
d$cnt <- rpois(n, exp(0.5 + 0.8 * d$lat))
d$bin <- rbinom(n, 1, plogis(0.3 + 1.2 * d$lat))
d$pos <- exp(0.3 + 0.8 * d$lat + rnorm(n, 0, 0.5))
d$ord <- as.integer(cut(0.8 * d$x + ifelse(d$z > 0, 1.2, -1.2) * rlogis(n),
                        c(-Inf, -0.5, 0.8, 2, Inf), labels = FALSE))
comp <- rbinom(n, 1, 0.5)
d$mixc <- rpois(n, ifelse(comp == 1, exp(1.5 + 0.4 * d$x), exp(0.2 - 0.3 * d$x)))
d$mixp <- ifelse(comp == 1, exp(1 + 0.5 * d$x + rnorm(n, 0, 0.3)), rgamma(n, 4, 4 / exp(-0.5 - 0.5 * d$x)))
d$cat <- as.integer(cut(0.7 * d$x + rlogis(n), c(-Inf, 0, 1.5, Inf), labels = FALSE))
d$mord <- as.integer(cut(ifelse(comp == 1, 0.8 * d$x + 1, -0.8 * d$x) + rlogis(n),
                         c(-Inf, -0.5, 0.8, 2, Inf), labels = FALSE))
msg <- function(e) {
  w <- character()
  v <- tryCatch(withCallingHandlers(e, warning = function(c) {
    w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")
  }, message = function(c) invokeRestart("muffleMessage")),
  error = function(err) structure("ERR", m = conditionMessage(err)))
  if (identical(c(v), "ERR")) return(paste("ERROR:", substr(attr(v, "m"), 1, 170)))
  s <- if (is.numeric(v)) paste(signif(as.numeric(v)[seq_len(min(8, length(v)))], 5), collapse = " ") else "OK"
  if (length(w)) s <- paste(s, "WARN:", substr(w[1L], 1, 110))
  s
}
fit_plain <- quiet(frm(bf(y ~ me(xo, sx) + z), family = gaussian(), data = d))
fit_gr <- quiet(frm(bf(yg ~ me(xg, sg, gr = g) + z), family = gaussian(), data = d))
cat("=== A. conditional_effects on me(gr =), noisy variable and default\n")
cat("CE xg     :", msg(conditional_effects(fit_gr, "xg", resolution = 4)[[1]]$estimate__), "\n")
cat("CE default:", msg(names(conditional_effects(fit_gr, resolution = 4))), "\n")
cat("CE z      :", msg(conditional_effects(fit_gr, "z", resolution = 4)[[1]]$estimate__), "\n")
cat("CE xg, conditions g:", msg(conditional_effects(fit_gr, "xg", resolution = 4, conditions = data.frame(g = "g01"))[[1]]$estimate__), "\n")
cat("predict CE xg:", msg(conditional_effects(fit_gr, "xg", resolution = 4, method = "predict")[[1]]$estimate__), "\n")
cat("frm_curve xg in one level:", msg(frm_curve(fit_gr, newdata = data.frame(xg = c(0, 1), sg = 0.4, z = 0, g = "g01"), simultaneous = FALSE)$.se), "\n")
if (!before) {
  ds_gr <- rev_draws_at(fit_gr, 5L)
  cat("draws CE xg:", msg(conditional_effects(ds_gr, "xg", resolution = 4)[[1]]$estimate__), "\n")
}
cat("\n=== B. conditional_effects(method = 'predict', seed =) under two global seeds\n")
fit_g <- quiet(frm(bf(y ~ z), family = gaussian(), data = d))
f1 <- {set.seed(1); conditional_effects(fit_g, "z", method = "predict", resolution = 3, seed = 3)[[1]]$estimate__}
f2 <- {set.seed(2); conditional_effects(fit_g, "z", method = "predict", resolution = 3, seed = 3)[[1]]$estimate__}
cat("plain gaussian: identical", identical(f1, f2), "\n")
f1 <- {set.seed(1); conditional_effects(fit_plain, "xo", method = "predict", resolution = 3, seed = 3)[[1]]$estimate__}
f2 <- {set.seed(2); conditional_effects(fit_plain, "xo", method = "predict", resolution = 3, seed = 3)[[1]]$estimate__}
cat("me: identical", identical(f1, f2), "\n")
cat("formals CE seed doc:", "seed" %in% names(formals(frmtmb:::conditional_effects.frmtmb_fit)), "\n")
cat("\n=== C. resp on a univariate fit: readers the lane did not name\n")
nd <- data.frame(xo = c(-1, 1), sx = 0.4, z = 0)
for (r in list(character(0), "y", c("y", "y"), "nosuch", NA_character_, 1L)) {
  lab <- paste(deparse(r), collapse = "")
  cat(lab, ": curve", substr(msg(frm_curve(fit_plain, newdata = nd, resp = r, simultaneous = FALSE)$.se), 1, 60), "\n")
  cat(lab, ": deriv", substr(msg(frm_curve_deriv(fit_plain, "xo", newdata = nd, resp = r)$.se), 1, 80), "\n")
  cat(lab, ": feature", substr(msg(frm_curve_feature(fit_g, "z", type = "level", at = 1, newdata = data.frame(z = seq(-2, 2, by = 0.5)), resp = r)$.estimate), 1, 80), "\n")
  cat(lab, ": get_predict", substr(msg(marginaleffects::get_predict(fit_plain, newdata = nd, resp = r)$estimate), 1, 80), "\n")
  cat(lab, ": emmeans", substr(msg(as.data.frame(emmeans::emmeans(fit_g, ~ z, at = list(z = 1), resp = r))$SE), 1, 90), "\n")
}
if (!before) {
  bg <- rev_pair("bf(y ~ z)", "gaussian()", d, ndraws = 20)
  for (r in list(character(0), "y", c("y", "y"), "nosuch", NA_character_, 1L)) {
    cat(paste(deparse(r), collapse = ""), ": brms emmeans", substr(msg(as.data.frame(emmeans::emmeans(bg$fx, ~ z, at = list(z = 1), resp = r))$emmean), 1, 100), "\n")
  }
}
cat("\n=== D. panel order and panels: both responses read x and z\n")
if (!before) {
  mvp <- rev_pair("mvbf(bf(y2 ~ x + z), bf(y3 ~ x + z + f), rescor = FALSE)", "gaussian()", d)
  md <- rev_draws_at(mvp$fit, 6L)
  cat("brms :", names(brms::conditional_effects(mvp$fx)), "\n")
  cat("fit  :", names(conditional_effects(mvp$fit, resolution = 3)), "\n")
  cat("draws:", names(conditional_effects(md, resolution = 3)), "\n")
  cat("brms effects x,z:", names(brms::conditional_effects(mvp$fx, c("x", "z"))), "\n")
  cat("fit  effects x,z:", names(conditional_effects(mvp$fit, c("x", "z"), resolution = 3)), "\n")
  cat("draws effects x,z:", names(conditional_effects(md, c("x", "z"), resolution = 3)), "\n")
  cat("brms one resp, x z:", names(brms::conditional_effects(mvp$fx, c("x", "z"), resp = "y3")), "\n")
  cat("fit one resp, x z:", names(conditional_effects(mvp$fit, c("x", "z"), resp = "y3", resolution = 3)), "\n")
  cat("fit print ok:", msg({capture.output(print(conditional_effects(mvp$fit, c("x"), resolution = 3))); 1}), "\n")
}
cat("\n=== E. layouts of fitted(scale = 'linear') and the linpred arrays\n")
lay <- function(a) paste(dim(a), collapse = "x")
layn <- function(a) paste(sapply(dimnames(a), function(v) if (is.null(v)) "NULL" else paste(head(v, 6), collapse = ",")), collapse = " | ")
if (!before) {
  for (sp in list(
    list("mix_pois", "bf(mixc ~ x)", "mixture(poisson(), poisson())"),
    list("categorical", "bf(cat ~ x)", "categorical()"))) {
    pr <- rev_pair(sp[[2]], sp[[3]], d)
    nd2 <- data.frame(x = c(-1, 1), z = 0)
    ds2 <- rev_draws_at(pr$fit, 3L)
    a <- fitted(pr$fx, newdata = nd2, scale = "linear")
    b <- fitted(pr$fit, newdata = nd2, scale = "linear")
    c1 <- frmtmb.sample:::fitted.frmtmb_draws(ds2, newdata = nd2, scale = "linear")
    cat(sp[[1]], "fitted linear brms", lay(a), layn(a), "\n")
    cat(sp[[1]], "fitted linear fit ", lay(b), layn(b), "\n")
    cat(sp[[1]], "fitted linear draw", lay(c1), layn(c1), "\n")
    a <- fitted(pr$fx, newdata = nd2, scale = "linear", summary = FALSE)
    c1 <- frmtmb.sample:::fitted.frmtmb_draws(ds2, newdata = nd2, scale = "linear", summary = FALSE)
    cat(sp[[1]], "unsummarized brms", lay(a), layn(a), "\n")
    cat(sp[[1]], "unsummarized draw", lay(c1), layn(c1), "\n")
    l1 <- frm_linpred(pr$fit, newdata = nd2)
    cat(sp[[1]], "frm_linpred", lay(l1), layn(l1), "\n")
    s1 <- frm_linpred(pr$fit, newdata = nd2, se.fit = TRUE)
    cat(sp[[1]], "se.fit dim", lay(s1$se.fit), " values ", signif(s1$se.fit, 4), "\n")
    pp <- marginaleffects::predictions(pr$fit, newdata = nd2, type = "link")
    cat(sp[[1]], "predictions link groups", paste(unique(pp$group), collapse = ","), " est ", signif(pp$estimate, 4), " se ", signif(pp$std.error, 4), "\n")
  }
}
cat("\n=== F. curve objects: transform under layers, a curve as first argument, print, plot\n")
if (!before) {
  fit_mp <- quiet(frm(bf(mixp ~ x), family = mixture(lognormal(), Gamma("log")), data = d))
  fit_mc <- quiet(frm(bf(mixc ~ x), family = mixture(poisson(), poisson()), data = d))
  fit_cat <- quiet(frm(bf(cat ~ x), family = categorical(), data = d))
  fit_cs <- quiet(frm(bf(ord ~ x + cs(z)), family = sratio(), data = d))
  gx <- data.frame(x = c(-1, 0, 1), z = 0.5)
  cat("mixture identity+log transform:", msg(frm_curve(fit_mp, newdata = gx, transform = TRUE, simultaneous = FALSE)$.estimate), "\n")
  cat("mixture log+log transform:", msg(frm_curve(fit_mc, newdata = gx, transform = TRUE, simultaneous = FALSE)$.estimate), "\n")
  cat("mixture link:", msg(frm_curve(fit_mc, newdata = gx, simultaneous = FALSE)$.estimate), "\n")
  cat("mixture mu2 transform:", msg(frm_curve(fit_mc, newdata = gx, dpar = "mu2", transform = TRUE, simultaneous = FALSE)$.estimate), "\n")
  cat("categorical transform:", msg(frm_curve(fit_cat, newdata = gx, transform = TRUE, simultaneous = FALSE)$.estimate), "\n")
  cat("cs transform:", msg(frm_curve(fit_cs, newdata = gx, transform = TRUE, simultaneous = FALSE)$.estimate), "\n")
  cs1 <- frm_curve(fit_cs, newdata = gx, seed = 1)
  cat("cs curve names:", names(cs1), "\n")
  cat("cs curve rows:", nrow(cs1), " .eta levels ", levels(cs1$.eta), "\n")
  print(cs1)
  cat("deriv of a curve object:", msg(frm_curve_deriv(cs1, "x")$.estimate), "\n")
  cat("deriv in z of a curve object:", msg(frm_curve_deriv(cs1, "z")$.estimate), "\n")
  cat("feature of a curve object:", msg(frm_curve_feature(cs1, "x", type = "level", at = 0)$.estimate), "\n")
  grDevices::pdf(NULL)
  cat("plot:", msg({p <- plot(cs1); 1}), "\n")
  cat("autoplot:", msg({p <- ggplot2::autoplot(cs1); 1}), "\n")
  grDevices::dev.off()
  cat("as.data.frame cols:", names(as.data.frame(cs1)), "\n")
  cat("frm_curve contrast cs:", msg(frm_curve(fit_cs, newdata = gx, contrast = transform(gx, z = -0.5), simultaneous = FALSE)$.se), "\n")
  fe <- tryCatch(quiet(frm_curve_feature(fit_cs, "x", type = "level", at = 0, newdata = data.frame(x = seq(-3, 3, by = 0.5), z = 0.5))), error = function(e) e)
  if (inherits(fe, "error")) cat("feature ERR", conditionMessage(fe), "\n") else print(as.data.frame(fe)[, c(".eta", ".estimate", ".se", ".value")])
  cat("mixture feature:", msg(frm_curve_feature(fit_mc, "x", type = "level", at = 1, newdata = data.frame(x = seq(-3, 3, by = 0.5)))$.estimate), "\n")
  cat("feature of curve with se reading .value_se:", msg(frm_curve_feature(fit_cs, "x", type = "level", at = 0, newdata = data.frame(x = seq(-3, 3, by = 0.5), z = 0.5))$.value_se), "\n")
}
cat("\n=== G. me() beside a cov = FALSE and a cov = TRUE term on simulate(newdata)\n")
d$yar <- d$y
for (g0 in levels(d$g)) {
  w <- which(d$g == g0)
  d$yar[w] <- 1 + 1.2 * d$lat[w] - 0.5 * d$z[w] + as.numeric(arima.sim(list(ar = 0.7), length(w), sd = 0.4))
}
ndar <- data.frame(xo = c(0.5, 0.5, 0.5, 0.5), sx = 0.3, z = 0, t = 1:4, g = "g01", yar = c(1, 2, 3, 4))
Par <- tryCatch(rev_pair("bf(yar ~ me(xo, sx) + z + ar(time = t, gr = g, cov = FALSE))", "gaussian()", d, ndraws = 3000), error = function(e) e)
if (inherits(Par, "error")) {
  cat("pair ERR:", conditionMessage(Par), "\n")
  far <- quiet(frm(bf(yar ~ me(xo, sx) + z + ar(time = t, gr = g, cov = FALSE)), family = gaussian(), data = d))
} else {
  far <- Par$fit
  bpr <- posterior_predict(Par$fx, newdata = ndar)
  cat("brms predict mean", signif(colMeans(bpr), 4), " sd ", signif(apply(bpr, 2, sd), 4), "\n")
  ber <- posterior_epred(Par$fx, newdata = ndar)
  cat("brms epred mean", signif(colMeans(ber), 4), " sd ", signif(apply(ber, 2, sd), 4), "\n")
}
set.seed(5)
s1 <- tryCatch(as.matrix(t(vapply(simulate(far, nsim = 3000, seed = 4, newdata = ndar), as.numeric, numeric(4)))), error = function(e) e)
if (!inherits(s1, "error")) cat("simulate mean", signif(colMeans(s1), 4), " sd ", signif(apply(s1, 2, sd), 4), "\n") else cat("simulate ERR", conditionMessage(s1), "\n")
p1 <- tryCatch(predict(far, newdata = ndar, summary = FALSE, ndraws = 3000, propagate_error = FALSE), error = function(e) e)
if (!inherits(p1, "error")) cat("predict mean", signif(colMeans(p1), 4), " sd ", signif(apply(p1, 2, sd), 4), "\n") else cat("predict ERR", conditionMessage(p1), "\n")
cat("frm_linpred ", msg(frm_linpred(far, newdata = ndar)), "\n")
cat("frm_linpred se", msg(frm_linpred(far, newdata = ndar, se.fit = TRUE)$se.fit), "\n")
cat("fitted Est.Error", msg(fitted(far, newdata = ndar)[, "Est.Error"]), "\n")
# the same grid without me(): the sum over the one-step mean, for reference
far0 <- quiet(frm(bf(yar ~ xo + z + ar(time = t, gr = g, cov = FALSE)), family = gaussian(), data = d))
s0 <- as.matrix(t(vapply(simulate(far0, nsim = 3000, seed = 4, newdata = ndar), as.numeric, numeric(4))))
cat("no me() simulate mean", signif(colMeans(s0), 4), " sd ", signif(apply(s0, 2, sd), 4), "\n")
cat("no me() frm_linpred", msg(frm_linpred(far0, newdata = ndar)), "\n")
cat("\n=== H. me() on a dpar's natural scale, and the mean through it\n")
Ps <- rev_pair("bf(y ~ z, sigma ~ me(xo, sx))", "gaussian()", d, ndraws = 4000)
ndh <- data.frame(xo = c(-1, 1, 1), sx = c(0.3, 0.3, 1), z = 0)
cat("brms epred sigma mean", signif(colMeans(posterior_epred(Ps$fx, newdata = ndh, dpar = "sigma")), 5), "\n")
cat("fitted sigma         ", signif(fitted(Ps$fit, newdata = ndh, dpar = "sigma")[, "Estimate"], 5), "\n")
cat("frm_linpred sigma se ", signif(frm_linpred(Ps$fit, newdata = ndh, dpar = "sigma", se.fit = TRUE)$se.fit, 5), "\n")
Pl <- rev_pair("bf(pos ~ me(xo, sx))", "lognormal()", d, ndraws = 4000)
cat("lognormal brms epred mean", signif(colMeans(posterior_epred(Pl$fx, newdata = ndh)), 5), "\n")
cat("lognormal fitted         ", signif(fitted(Pl$fit, newdata = ndh)[, "Estimate"], 5), "\n")
cat("\n=== N. fitted() interval and support, me() on the response scale\n")
for (sp in list(list("bf(cnt ~ me(xo, sx))", "poisson()"), list("bf(bin ~ me(xo, sx))", "bernoulli()"))) {
  ff <- quiet(frm(rev_ev(sp[[1]], "frmtmb"), family = rev_ev(sp[[2]], "frmtmb"), data = d))
  ndn <- data.frame(xo = c(0, 1, 1, 1.5), sx = c(0.3, 0.8, 1.2, 1.2))
  cat(sp[[2]], "fitted\n")
  print(round(fitted(ff, newdata = ndn), 3))
  cat(" frm_curve transform band\n")
  print(round(as.data.frame(frm_curve(ff, newdata = ndn, transform = TRUE, simultaneous = FALSE))[, c(".estimate", ".se", ".lower_ci", ".upper_ci")], 3))
  ce <- conditional_effects(ff, "xo", resolution = 3, conditions = data.frame(sx = 1.2))[[1]]
  print(round(as.data.frame(ce)[, c("xo", "estimate__", "se__", "lower__", "upper__")], 3))
}

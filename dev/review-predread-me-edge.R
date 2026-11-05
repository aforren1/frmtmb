# Item 4: me() edge cases and readers. Data seed 7401.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7401)
n <- 200
d <- data.frame(z = rnorm(n), g = factor(rep(sprintf("g%02d", 1:20), each = 10)))
d$lat <- rnorm(n)
d$sx <- 0.4
d$xo <- d$lat + rnorm(n, 0, d$sx)
gl <- rnorm(20)[d$g]
d$sg <- 0.4
d$xg <- gl + rnorm(20, 0, 0.4)[d$g]
d$y <- 1 + 1.2 * d$lat - 0.5 * d$z + rnorm(n, 0, 0.5)
d$yg <- 1 + 1.2 * gl - 0.5 * d$z + rnorm(n, 0, 0.5)
d$u <- factor(rep(sprintf("u%02d", 1:20), 10))
d$yu <- d$y + rnorm(20, 0, 0.7)[d$u]
d$cnt <- rpois(n, exp(0.5 + 0.8 * d$lat))
P <- rev_pair("bf(y ~ me(xo, sx) + z)", "gaussian()", d, ndraws = 2000)
G <- rev_pair("bf(yg ~ me(xg, sg, gr = g) + z)", "gaussian()", d, ndraws = 2000)
U <- rev_pair("bf(yu ~ me(xo, sx) + z + (1 | u))", "gaussian()", d, ndraws = 2000)
Pds <- rev_draws_at(P$fit, 2000L)
Gds <- rev_draws_at(G$fit, 2000L)
Uds <- rev_draws_at(U$fit, 2000L)
msg <- function(e) {
  w <- character()
  v <- tryCatch(withCallingHandlers(e, warning = function(c) {
    w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")
  }, message = function(c) invokeRestart("muffleMessage")),
  error = function(err) structure("ERR", m = conditionMessage(err)))
  if (identical(c(v), "ERR")) {
    return(paste("ERROR:", substr(attr(v, "m"), 1, 150)))
  }
  s <- if (is.numeric(v)) {
    paste("OK", if (is.null(dim(v))) length(v) else paste(dim(v), collapse = "x"),
          "NA:", sum(is.na(v)))
  } else {
    "OK"
  }
  if (length(w)) s <- paste(s, "WARN:", substr(w[1L], 1, 110))
  s
}
readers <- function(pr, ds, nd, name, ...) {
  cat("\n--- ", name, "\n")
  cat("brms posterior_epred : ", msg(posterior_epred(pr$fx, newdata = nd, ...)), "\n")
  cat("brms posterior_predict: ", msg(posterior_predict(pr$fx, newdata = nd, ...)), "\n")
  cat("brms fitted          : ", msg(fitted(pr$fx, newdata = nd, ...)), "\n")
  cat("draws posterior_epred: ", msg(frmtmb.sample::posterior_epred(ds, newdata = nd, ...)), "\n")
  cat("draws posterior_predict:", msg(frmtmb.sample::posterior_predict(ds, newdata = nd, ...)), "\n")
  cat("fit fitted           : ", msg(fitted(pr$fit, newdata = nd, ...)), "\n")
  cat("fit frm_linpred      : ", msg(frm_linpred(pr$fit, newdata = nd, ...)), "\n")
  cat("fit frm_linpred se   : ", msg(frm_linpred(pr$fit, newdata = nd, se.fit = TRUE, ...)$se.fit), "\n")
  cat("fit predict          : ", msg(predict(pr$fit, newdata = nd, ndraws = 20, ...)), "\n")
  cat("fit simulate         : ", msg(as.matrix(simulate(pr$fit, nsim = 3, seed = 1, newdata = nd, ...))), "\n")
  cat("fit frm_curve        : ", msg(frm_curve(pr$fit, newdata = nd, simultaneous = FALSE, ...)$.se), "\n")
}
nd0 <- data.frame(xo = c(-1, 0.5, 1), sx = c(0.3, 0.6, 0.9), z = 0)
readers(P, Pds, nd0, "complete")
nd <- nd0
nd$sx[2] <- 0
readers(P, Pds, nd, "sdx = 0 in one row")
nd <- nd0
nd$sx[2] <- -1
readers(P, Pds, nd, "sdx < 0 in one row")
nd <- nd0
nd$sx[2] <- NA
readers(P, Pds, nd, "NA sdx in one row")
nd <- nd0
nd$xo[2] <- NA
readers(P, Pds, nd, "NA x in one row")
nd <- nd0
nd$sx <- NULL
readers(P, Pds, nd, "no sdx column")
nd <- data.frame(xg = c(0.2, 0.2, 1), sg = c(0.4, 0.4, 0.8), z = 0,
                 g = c("g01", "g01", "g02"))
readers(G, Gds, nd, "gr: complete")
nd2 <- nd
nd2$g <- c("g01", "g01", "new")
readers(G, Gds, nd2, "gr: new level")
cat("brms allow_new_levels: ", msg(posterior_epred(G$fx, newdata = nd2, allow_new_levels = TRUE)), "\n")
cat("fit allow_new_levels : ", msg(frm_linpred(G$fit, newdata = nd2, allow_new_levels = TRUE, se.fit = TRUE)$se.fit), "\n")
cat("fit predict allow    : ", msg(predict(G$fit, newdata = nd2, ndraws = 20, allow_new_levels = TRUE)), "\n")
cat("draws allow          : ", msg(frmtmb.sample::posterior_epred(Gds, newdata = nd2, allow_new_levels = TRUE)), "\n")
nd2 <- nd
nd2$xg[2] <- 0.3
readers(G, Gds, nd2, "gr: x varies in a level")
nd2 <- nd
nd2$g <- NULL
readers(G, Gds, nd2, "gr: no g column")
nd2 <- nd
nd2$g <- c("g01", NA, "g02")
nd2$xg <- c(.2, .3, 1)
readers(G, Gds, nd2, "gr: NA level (character)")
nd2$g <- factor(nd2$g)
readers(G, Gds, nd2, "gr: NA level (factor)")

cat("\n--- gr: sharing, per row, on draws and predict\n")
ndg <- data.frame(xg = c(0.2, 0.2, 1, 0.2), sg = c(0.4, 0.4, 0.8, 0.4), z = 0,
                  g = c("g01", "g01", "g02", "g01"))
set.seed(7402)
e <- frmtmb.sample::posterior_epred(Gds, newdata = ndg)
cat("draws: identical rows 1,2,4:", identical(e[, 1], e[, 2]), identical(e[, 1], e[, 4]),
    " cor(1,3)", round(cor(e[, 1], e[, 3]), 4), "\n")
b <- posterior_epred(G$fx, newdata = ndg)
cat("brms : identical rows 1,2,4:", identical(b[, 1], b[, 2]), identical(b[, 1], b[, 4]),
    " cor(1,3)", round(cor(b[, 1], b[, 3]), 4), "\n")
L <- frm_lp_basis(G$fit, newdata = ndg, extra_cov = TRUE)
cat("basis extra_cov:\n")
print(round(as.matrix(L$extra_cov), 4))
cat("cov of brms draws:\n")
print(round(cov(b), 4))
cat("cov of draws:\n")
print(round(cov(e), 4))
set.seed(1)
o1 <- frmtmb.sample::posterior_epred(Gds, newdata = ndg)
set.seed(1)
o2 <- frmtmb.sample::posterior_epred(Gds, newdata = ndg[c(3, 1, 2, 4), ])
cat("order: sd of each row under the two orders",
    round(apply(o1, 2, sd), 4), "/", round(apply(o2, 2, sd), 4), "\n")
cat("\n--- reproducibility\n")
same <- function(a, b) identical(a, b)
cat("predict(seed) twice:", same({set.seed(5); predict(P$fit, newdata = nd0, ndraws = 50)},
                                 {set.seed(5); predict(P$fit, newdata = nd0, ndraws = 50)}), "\n")
cat("predict(summary=FALSE, seed) twice:",
    same({set.seed(5); predict(P$fit, newdata = nd0, ndraws = 50, summary = FALSE)},
         {set.seed(5); predict(P$fit, newdata = nd0, ndraws = 50, summary = FALSE)}), "\n")
cat("simulate(seed) twice:", same(simulate(P$fit, nsim = 3, seed = 5, newdata = nd0),
                                  simulate(P$fit, nsim = 3, seed = 5, newdata = nd0)), "\n")
f1 <- {set.seed(9); frmtmb.sample::posterior_epred(Pds, newdata = nd0)}
f2 <- {set.seed(9); frmtmb.sample::posterior_epred(Pds, newdata = nd0)}
cat("posterior_epred set.seed twice:", same(f1, f2), "\n")
f1 <- {set.seed(9); frmtmb.sample::posterior_predict(Pds, newdata = nd0)}
f2 <- {set.seed(9); frmtmb.sample::posterior_predict(Pds, newdata = nd0)}
cat("posterior_predict set.seed twice:", same(f1, f2), "\n")
cat("two different seeds differ:",
    !same({set.seed(5); predict(P$fit, newdata = nd0, ndraws = 50)},
          {set.seed(6); predict(P$fit, newdata = nd0, ndraws = 50)}), "\n")
f1 <- {set.seed(1); conditional_effects(P$fit, "xo", method = "predict", resolution = 4, seed = 3)[[1]]$estimate__}
f2 <- {set.seed(2); conditional_effects(P$fit, "xo", method = "predict", resolution = 4, seed = 3)[[1]]$estimate__}
cat("CE predict same seed arg, global seed differs:", same(f1, f2), "\n")
set.seed(4)
a <- runif(1)
set.seed(4)
invisible(fitted(P$fit, newdata = nd0))
b <- runif(1)
cat("fitted(newdata) leaves the stream alone:", a == b, "\n")
set.seed(4)
a <- runif(1)
set.seed(4)
invisible(predict(P$fit, newdata = nd0, ndraws = 3))
b <- runif(1)
cat("predict(seed=) leaves the global stream alone:", a == b, "\n")

cat("\n--- training rows: as newdata and in sample, with brms\n")
tr <- d[1:4, ]
b1 <- quiet(posterior_epred(P$fx))
b2 <- posterior_epred(P$fx, newdata = tr)
cat("brms in-sample epred sd rows 1:4:", round(apply(b1, 2, sd)[1:4], 4), " (no save_pars)\n")
cat("brms newdata=train epred sd rows 1:4:", round(apply(b2, 2, sd), 4), "\n")
cat("brms in-sample mean", round(colMeans(b1)[1:4], 4), "\n")
cat("frmtmb fitted in-sample", round(fitted(P$fit)[1:4, "Estimate"], 4), " se ",
    round(fitted(P$fit)[1:4, "Est.Error"], 4), "\n")
cat("frmtmb fitted newdata=train", round(fitted(P$fit, newdata = tr)[, "Estimate"], 4),
    " se ", round(fitted(P$fit, newdata = tr)[, "Est.Error"], 4), "\n")
set.seed(1); ip <- predict(P$fit, summary = FALSE, ndraws = 500, propagate_error = FALSE)
cat("frmtmb predict in-sample sd", round(apply(ip, 2, sd)[1:4], 4), " sigma",
    round(exp(P$fit$estimates$betad[1]), 4), "\n")
set.seed(1); ip2 <- predict(P$fit, newdata = tr, summary = FALSE, ndraws = 500, propagate_error = FALSE)
cat("frmtmb predict newdata=train sd", round(apply(ip2, 2, sd)[1:4], 4), "\n")

cat("\n--- re_formula\n")
ndu <- data.frame(xo = c(-1, 0.5, 1), sx = c(0.3, 0.6, 0.9), z = 0, u = c("u01", "u02", "u03"))
for (rf in list(NULL, NA, ~ (1 | u))) {
  lab <- if (is.null(rf)) "NULL" else if (is.na(rf)[1]) "NA" else "~(1|u)"
  be <- posterior_epred(U$fx, newdata = ndu, re_formula = rf)
  fe <- frmtmb.sample::posterior_epred(Uds, newdata = ndu, re_formula = rf)
  pl <- frm_linpred(U$fit, newdata = ndu, re_formula = rf, se.fit = TRUE)
  fp <- predict(U$fit, newdata = ndu, re_formula = rf, summary = FALSE, ndraws = 2000,
                propagate_error = FALSE)
  bp <- posterior_predict(U$fx, newdata = ndu, re_formula = rf)
  cat(lab, ": brms mean", round(colMeans(be), 3), " sd", round(apply(be, 2, sd), 3), "\n")
  cat(lab, ": draws mean", round(colMeans(fe), 3), " sd", round(apply(fe, 2, sd), 3), "\n")
  cat(lab, ": frm_linpred", round(pl$fit, 3), " se", round(pl$se.fit, 3), "\n")
  cat(lab, ": pp brms mean", round(colMeans(bp), 3), " sd", round(apply(bp, 2, sd), 3), "\n")
  cat(lab, ": predict mean", round(colMeans(fp), 3), " sd", round(apply(fp, 2, sd), 3), "\n")
}

cat("\n--- other readers on a me() fit: double count or missing\n")
fit <- P$fit
ndn <- data.frame(xo = c(-1, 0.5, 1, 1.5), sx = c(0.3, 0.6, 0.9, 0.9), z = c(0, 1, 0, 0))
pl <- frm_linpred(fit, newdata = ndn, se.fit = TRUE)
lb <- frm_lp_basis(fit, newdata = ndn)
V <- as.matrix(lb$V)
rowvar <- rowSums((lb$A %*% V) * lb$A)
cf <- summary(fit)$fixed
print(round(cf, 4))
bsp <- cf[grep("me", rownames(cf))[1], 1]
cat("(b sdx):", round(abs(bsp) * ndn$sx, 4), "\n")
cat("frm_linpred se      ", round(pl$se.fit, 4), "\n")
cat("basis sqrt(extra)   ", round(sqrt(lb$extra_var), 4), "\n")
cat("basis sqrt(AVA')    ", round(sqrt(pmax(rowvar, 0)), 4), "\n")
cat("sqrt(AVA'+extra)    ", round(sqrt(rowvar + lb$extra_var), 4), "\n")
fv <- fitted(fit, newdata = ndn)
cat("fitted Est.Error    ", round(fv[, "Est.Error"], 4), "\n")
cu <- frm_curve(fit, newdata = ndn, simultaneous = TRUE, seed = 1)
cat("frm_curve .se       ", round(cu$.se, 4), " crit ", round(cu$.crit[1], 3), " crit_sim ",
    round(cu$.crit_sim[1], 3), "\n")
S <- attr(cu, "Sigma")
cat("Sigma off-diag max  ", round(max(abs(S - diag(diag(S)))), 4), " AVA' off-diag max ",
    round(max(abs((lb$A %*% V %*% t(lb$A)) - diag(rowvar))), 4), "\n")
pr <- marginaleffects::predictions(fit, newdata = ndn, re_formula = NA)
cat("marginaleffects est ", round(pr$estimate, 4), "\n")
cat("marginaleffects se  ", round(pr$std.error, 4), "\n")
ame <- tryCatch(marginaleffects::avg_slopes(fit, variables = "xo", newdata = ndn, re_formula = NA),
                error = function(e) e)
if (inherits(ame, "error")) {
  cat("avg_slopes ERROR", conditionMessage(ame), "\n")
} else {
  cat("avg_slopes xo", round(ame$estimate, 4), " se ", round(ame$std.error, 4), "\n")
}
em <- as.data.frame(emmeans::emmeans(fit, ~ xo, at = list(xo = c(-1, 1), z = 0, sx = 0.6)))
cat("emmeans xo emmean", round(em$emmean, 4), " SE ", round(em$SE, 4), "\n")
emb <- as.data.frame(emmeans::emmeans(P$fx, ~ xo, at = list(xo = c(-1, 1), z = 0, sx = 0.6)))
cat("brms emmeans ", round(emb$emmean, 4), " L/U ", round(emb$lower.HPD, 3), round(emb$upper.HPD, 3), "\n")
cn <- as.data.frame(emmeans::contrast(emmeans::emmeans(fit, ~ xo, at = list(xo = c(-1, 1), z = 0, sx = 0.6)), "revpairwise"))
cnb <- as.data.frame(emmeans::contrast(emmeans::emmeans(P$fx, ~ xo, at = list(xo = c(-1, 1), z = 0, sx = 0.6)), "revpairwise"))
cat("emmeans contrast", round(cn$estimate, 4), " SE ", round(cn$SE, 4), "\n")
cat("brms contrast", round(cnb$estimate, 4), " L/U ", round(cnb$lower.HPD, 3), round(cnb$upper.HPD, 3), "\n")
ema <- as.data.frame(emmeans::emmeans(fit, ~ z, at = list(z = c(0, 1))))
emab <- as.data.frame(emmeans::emmeans(P$fx, ~ z, at = list(z = c(0, 1))))
cat("emmeans z avg emmean", round(ema$emmean, 4), " SE ", round(ema$SE, 4), "\n")
cat("brms z avg", round(emab$emmean, 4), " L/U ", round(emab$lower.HPD, 3), round(emab$upper.HPD, 3), "\n")
emc <- as.data.frame(emmeans::contrast(emmeans::emmeans(fit, ~ z, at = list(z = c(0, 1))), "revpairwise"))
emcb <- as.data.frame(emmeans::contrast(emmeans::emmeans(P$fx, ~ z, at = list(z = c(0, 1))), "revpairwise"))
cat("emm z contrast", round(emc$estimate, 4), " SE ", round(emc$SE, 4),
    " (fixed z, se in summary above); brms contrast HPD", round(emcb$lower.HPD, 3), round(emcb$upper.HPD, 3), "\n")
ctr <- frm_curve(fit, newdata = data.frame(xo = c(1, 1), sx = 0.6, z = c(0, 1)),
                 contrast = data.frame(xo = c(1, 1), sx = 0.6, z = c(1, 0)), simultaneous = FALSE)
cat("frm_curve contrast same xo: est", round(ctr$.estimate, 4), " se ", round(ctr$.se, 4),
    " (coef z se ", round(cf["z", 2], 4), ")\n")
der <- try_num(frm_curve_deriv(fit, "xo", newdata = data.frame(xo = c(-1, 0.5), sx = 0.6, z = 0)))
if (is.data.frame(der)) {
  print(as.data.frame(der)[, c("xo", ".estimate", ".se")])
} else {
  cat("deriv ERR", attr(der, "err"), "\n")
}
cat("frm_curve_deriv under gr xg: ", msg(frm_curve_deriv(G$fit, "xg", newdata = data.frame(xg = c(0.2, 1), sg = 0.4, z = 0, g = c("g01", "g02")))$.se), "\n")
cat("frm_curve_deriv under gr z: ", msg(frm_curve_deriv(G$fit, "z", newdata = data.frame(xg = c(0.2, 1), sg = 0.4, z = 0, g = c("g01", "g02")))$.se), "\n")
cat("CE fit xo: ", msg(conditional_effects(fit, "xo", resolution = 4)[[1]]$se__), "\n")
cat("brms CE xo: ", msg(brms::conditional_effects(P$fx, "xo", resolution = 4)[[1]]$se__), "\n")
cat("CE fit gr xg: ", msg(conditional_effects(G$fit, "xg", resolution = 4)[[1]]$se__), "\n")
cat("brms CE gr xg: ", msg(brms::conditional_effects(G$fx, "xg", resolution = 4)[[1]]$se__), "\n")
cat("emmeans gr z: ", msg(as.data.frame(emmeans::emmeans(G$fit, ~ z, at = list(z = 0:1)))$SE), "\n")
cat("brms emmeans gr z: ", msg(as.data.frame(emmeans::emmeans(G$fx, ~ z, at = list(z = 0:1)))$emmean), "\n")
cat("emmeans gr xg: ", msg(as.data.frame(emmeans::emmeans(G$fit, ~ xg, at = list(xg = c(-1, 1))))$SE), "\n")
cat("brms emmeans gr xg: ", msg(as.data.frame(emmeans::emmeans(G$fx, ~ xg, at = list(xg = c(-1, 1))))$emmean), "\n")
ce <- conditional_effects(fit, "xo", resolution = 4, method = "predict", seed = 1)[[1]]
ceb <- brms::conditional_effects(P$fx, "xo", resolution = 4, method = "posterior_predict")[[1]]
cat("CE predict width", round((ce$upper__ - ce$lower__), 3), " brms ", round((ceb$upper__ - ceb$lower__), 3), "\n")
ce <- conditional_effects(fit, "xo", resolution = 4)[[1]]
ceb <- brms::conditional_effects(P$fx, "xo", resolution = 4)[[1]]
cat("CE fitted width", round((ce$upper__ - ce$lower__), 3), " brms ", round((ceb$upper__ - ceb$lower__), 3),
    " (brms: no coefficient uncertainty here)\n")
cat("sx in CE grid: frmtmb", ce$sx[1], " brms ", ceb$sx[1], "\n")
ins <- tryCatch(insight::get_predicted(fit, data = ndn), error = function(e) e)
cat("insight get_predicted:", if (inherits(ins, "error")) conditionMessage(ins) else round(as.numeric(ins), 4), "\n")
if (!inherits(ins, "error")) cat(" SE ", round(attr(ins, "SE") %||% NA, 4), "\n")

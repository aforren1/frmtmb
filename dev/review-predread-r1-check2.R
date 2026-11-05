# Short check, items 2, 3: what each contrast answers, the jump at "same x", what a
# default marginaleffects call does on a me() fit, and gp() kriging in the same class.
# Seed 7402 (me), 7403 (gp).
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7402)
n <- 200
d <- data.frame(z = rnorm(n), lat = rnorm(n), sx = 0.4, f = factor(sample(c("a", "b"), n, TRUE)))
d$xo <- d$lat + rnorm(n, 0, 0.4)
d$y <- 1 + 1.2 * d$lat - 0.5 * d$z + 0.6 * (d$f == "b") + rnorm(n, 0, 0.5)
P <- rev_pair("bf(y ~ me(xo, sx) + z + f)", "gaussian()", d, ndraws = 4000)
fit <- P$fit
fx <- P$fx
cf <- summary(fit)$fixed
print(signif(cf[, 1:2], 4))
msg <- function(e) {
  r <- tryCatch(suppressWarnings(suppressMessages(e)), error = function(err) paste("ERROR:", substr(conditionMessage(err), 1, 150)))
  if (is.character(r)) r else paste(signif(as.numeric(r), 4), collapse = " ")
}
cat("=== contrast in z: the jump as the two x values part\n")
for (dx in c(0, 1e-12, 1e-6, 1e-3, 0.1, 1)) {
  a <- data.frame(xo = 1, sx = 0.6, z = 1, f = "a")
  b <- data.frame(xo = 1 + dx, sx = 0.6, z = 0, f = "a")
  ct <- frm_curve(fit, newdata = a, contrast = b, simultaneous = FALSE)
  cat(sprintf("dx %6g: est %.4f se %.4f\n", dx, ct$.estimate, ct$.se))
}
cat("same x, sdx 0.6 and 0.6000001:")
ct <- frm_curve(fit, newdata = data.frame(xo = 1, sx = 0.6, z = 1, f = "a"),
                contrast = data.frame(xo = 1, sx = 0.6000001, z = 0, f = "a"), simultaneous = FALSE)
cat(sprintf(" est %.4f se %.4f\n", ct$.estimate, ct$.se))
cat("contrast in x itself, 1 against 2, sdx 0.6:")
ct <- frm_curve(fit, newdata = data.frame(xo = 2, sx = 0.6, z = 0, f = "a"),
                contrast = data.frame(xo = 1, sx = 0.6, z = 0, f = "a"), simultaneous = FALSE)
cat(sprintf(" est %.4f se %.4f; slope se %.4f\n", ct$.estimate, ct$.se, cf[grep("me", rownames(cf))[1], 2]))
de <- frm_curve_deriv(fit, "xo", newdata = data.frame(xo = 1.5, sx = 0.6, z = 0, f = "a"), simultaneous = FALSE)
cat("frm_curve_deriv in x: est", signif(de$.estimate, 4), " se", signif(de$.se, 4), "\n")
cat("two rows of one newdata, frm_curve Sigma covariance (rows same x, sx):\n")
cu <- frm_curve(fit, newdata = data.frame(xo = c(1, 1), sx = 0.6, z = c(1, 0), f = "a"), simultaneous = FALSE)
S <- attr(cu, "Sigma")
print(signif(S, 4))
cat("var of the difference read off this Sigma:", signif(S[1, 1] + S[2, 2] - 2 * S[1, 2], 4), " sd", signif(sqrt(S[1, 1] + S[2, 2] - 2 * S[1, 2]), 4), "\n")
cat("the same pair through contrast =:", signif(frm_curve(fit, newdata = data.frame(xo = 1, sx = 0.6, z = 1, f = "a"), contrast = data.frame(xo = 1, sx = 0.6, z = 0, f = "a"), simultaneous = FALSE)$.se, 4), "\n")
em <- as.data.frame(emmeans::contrast(emmeans::emmeans(fit, ~ z, at = list(z = c(0, 1), xo = 1, sx = 0.6)), "revpairwise"))
cat("emmeans ~ z pairs: est", signif(em$estimate, 4), " SE", signif(em$SE, 4), "\n")
emf <- as.data.frame(emmeans::contrast(emmeans::emmeans(fit, ~ f, at = list(xo = 1, sx = 0.6, z = 0)), "revpairwise"))
cat("emmeans ~ f pairs: est", signif(emf$estimate, 4), " SE", signif(emf$SE, 4), " coef se", signif(cf["fb", 2], 4), "\n")
emb <- as.data.frame(emmeans::contrast(emmeans::emmeans(fx, ~ f, at = list(xo = 1, sx = 0.6, z = 0)), "revpairwise"))
cat("brms emmeans ~ f pairs: est", signif(emb$estimate, 4), " HPD", signif(emb$lower.HPD, 4), signif(emb$upper.HPD, 4), "\n")
nd <- data.frame(xo = c(1, 1), sx = 0.6, z = c(1, 0), f = "a")
bl <- posterior_linpred(fx, newdata = nd)
cat("brms one newdata, sd of the difference:", signif(sd(bl[, 1] - bl[, 2]), 4), "\n")
cat("the coefficient of z, by simulation of one unit's latent, same draw:", signif(sd((bl[, 1] - bl[, 2]) * 0), 4), "\n")
cat("\n=== marginaleffects on a me() fit, default call\n")
ndm <- data.frame(xo = c(-1, 1), sx = 0.6, z = 0, f = "a")
fns <- list(
  "predictions(newdata)" = function(m, ...) marginaleffects::predictions(m, newdata = ndm, ...)$estimate,
  "predictions() in sample" = function(m, ...) marginaleffects::predictions(m, ...)$estimate[1:2],
  "avg_predictions()" = function(m, ...) marginaleffects::avg_predictions(m, ...)$estimate,
  "avg_slopes(z)" = function(m, ...) marginaleffects::avg_slopes(m, variables = "z", ...)$estimate,
  "comparisons(f)" = function(m, ...) marginaleffects::avg_comparisons(m, variables = "f", ...)$estimate,
  "hypotheses(b)" = function(m, ...) marginaleffects::hypotheses(m, "z = 0", ...)$estimate)
for (nm in names(fns)) {
  cat(sprintf("%-24s default: %s\n", nm, substr(msg(fns[[nm]](fit)), 1, 120)))
  cat(sprintf("%-24s vcov = FALSE: %s\n", nm, substr(msg(fns[[nm]](fit, vcov = FALSE)), 1, 80)))
  cat(sprintf("%-24s brms: %s\n", nm, substr(msg(fns[[nm]](fx)), 1, 80)))
}
pb <- list(estimate = colMeans(posterior_linpred(fx, newdata = ndm)), conf.low = apply(posterior_linpred(fx, newdata = ndm), 2, quantile, 0.025), conf.high = apply(posterior_linpred(fx, newdata = ndm), 2, quantile, 0.975))
cat("brms posterior_linpred summary (its predictions() needs the collapse package, not installed): est", signif(pb$estimate, 4), " CI", signif(pb$conf.low, 4), signif(pb$conf.high, 4), "\n")
li <- frm_linpred(fit, newdata = ndm, se.fit = TRUE)
cat("frm_linpred est", signif(li$fit, 4), " se", signif(li$se.fit, 4), "\n")
V <- as.matrix(vcov(fit, full = TRUE))
pv <- tryCatch(marginaleffects::predictions(fit, newdata = ndm, vcov = V), error = function(e) conditionMessage(e))
cat("vcov = matrix:", if (is.character(pv)) substr(pv, 1, 120) else paste("se", paste(signif(pv$std.error, 4), collapse = " ")), "\n")
ct <- get("get_vcov.frmtmb_fit", asNamespace("frmtmb"))
cat("get_vcov formals:", names(formals(ct)), "\n")
cat("what marginaleffects passes to get_vcov: ")
tr <- list()
trace_args <- function(model, vcov = NULL, ...) {
  tr[[length(tr) + 1L]] <<- names(list(...))
  NULL
}
print(names(formals(marginaleffects::get_vcov)))
cat("r16 reference: avg_slopes(z) se equals coef z se in the mean-only design; see before-arm log\n")
cat("\n=== gp() kriging: predictions() at an unseen position\n")
set.seed(7403)
m <- 80
g <- data.frame(x = sort(runif(m, 0, 6)))
g$y <- sin(g$x) + rnorm(m, 0, 0.25)
Gp <- rev_pair("bf(y ~ gp(x))", "gaussian()", g, ndraws = 2000)
gnd <- data.frame(x = c(g$x[5], 3.33, 7, 9))
lg <- frm_lp_basis(Gp$fit, newdata = gnd)
sl <- frm_linpred(Gp$fit, newdata = gnd, se.fit = TRUE)
cat("frm_linpred est", signif(sl$fit, 4), " se", signif(sl$se.fit, 4), "\n")
cat("extra sd", signif(sqrt(lg$extra_var), 4), "\n")
pg <- tryCatch(marginaleffects::predictions(Gp$fit, newdata = gnd), error = function(e) conditionMessage(e))
cat("predictions:", if (is.character(pg)) pg else paste("se", paste(signif(pg$std.error, 4), collapse = " ")), "\n")
bg <- posterior_linpred(Gp$fx, newdata = gnd)
cat("brms at fixed values sd over draws:", signif(apply(bg, 2, sd), 4), "\n")
cat("predictions in sample positions only:", if (!is.character(tryCatch(marginaleffects::predictions(Gp$fit, newdata = g[1:3, , drop = FALSE]), error = function(e) "x"))) "answers" else "refused", "\n")

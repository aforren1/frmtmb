# Item 4, second part: by, new data, re_formula, new levels, prior tables,
# level-order checks by hand. Seed 9820.
source("dev/review-refuse-lib.R")
set.seed(9820)
n <- 150
lev <- paste0("L", sprintf("%02d", 1:12))
g1 <- sample(lev[1:10], n, TRUE, prob = c(0.35, rep(0.65 / 9, 9)))
g2 <- sample(lev[2:12], n, TRUE, prob = c(rep(0.6 / 9, 9), 0.2, 0.2))
d <- data.frame(g1 = g1, g2 = g2, x = rnorm(n), w = rnorm(n),
                w1 = runif(n, 0.2, 1), w2 = runif(n, 0.2, 1),
                v1 = runif(n, 0.2, 1), v2 = runif(n, 0.2, 1),
                x1 = rnorm(n), x2 = rnorm(n),
                f1 = factor(sample(c("a", "b"), n, TRUE)))
d$f2 <- d$f1
re <- rnorm(12, 0, 0.8); names(re) <- lev; res <- rnorm(12, 0, 0.3); names(res) <- lev
d$y <- 1 + 0.5 * d$x + (d$w1 * re[d$g1] + d$w2 * re[d$g2]) / (d$w1 + d$w2) +
  rnorm(n, 0, exp(-0.2 + (res[d$g1] + res[d$g2]) / 2))
show <- function(lab, r) {
  if (is_err(r)) { pr("%-40s ERROR %s", lab, msg(r, 170)); return(invisible()) }
  pr("%-40s verdict=%s const=%s pts=[%s] grad=%s errs: frm[%s] brms[%s] l1[%s]", lab, r$verdict,
     format(r$l1_const %||% NA, digits = 3),
     paste(format(unlist(r$l1_extra_const) %||% NA, digits = 3), collapse = ","),
     format(r$l1_grad %||% NA, digits = 3),
     substr(r$frm_error %||% "", 1, 110), substr(r$brms_error %||% "", 1, 110),
     substr(r$l1_error %||% "", 1, 110))
}
sp <- function(bf, data = d, ...) c(list(bf = bf, data = data), list(...))
A0 <- exp(-abs(outer(1:12, 1:12, "-")) / 3); dimnames(A0) <- list(lev, lev)
d$f1 <- factor(ifelse(d$g1 %in% lev[1:6], "a", "b")); d$f2 <- factor(ifelse(d$g2 %in% lev[1:6], "a", "b"))
show("mm by matrix control", tr(meas(sp("bf(y ~ x + (1 | mm(g1, g2, by = cbind(f1, f2))))"))))
show("|ID| with by matrix", tr(meas(sp(
  "bf(y ~ x + (1 | p | mm(g1, g2, by = cbind(f1, f2))), sigma ~ 1 + (1 | p | mm(g1, g2, by = cbind(f1, f2))))"))))
show("|ID| with by, id =", tr(meas(sp(
  "bf(y ~ x + (1 | mm(g1, g2, by = cbind(f1, f2), id = 'p')), sigma ~ 1 + (1 | mm(g1, g2, by = cbind(f1, f2), id = 'p')))"))))

fit <- frm(bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w1, w2))),
              sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(v1, v2)))), data = d)
re1 <- ranef(fit)[[1]]
cf <- fixef(fit)
nd <- d[1:6, ]
ft <- fitted(fit, newdata = nd, dpar = "mu")
hmu <- function(nd, r = re1[, 1]) cf["Intercept", 1] + cf["x", 1] * nd$x +
  (nd$w1 * r[nd$g1] + nd$w2 * r[nd$g2]) / (nd$w1 + nd$w2)
pr("hand mu vs fitted(newdata): %.3e", max(abs(hmu(nd) - ft[, 1])))
hsg <- function(nd, r = re1[, 2]) exp(cf["sigma_Intercept", 1] + (nd$v1 * r[nd$g1] + nd$v2 * r[nd$g2]) / (nd$v1 + nd$v2))
fs <- fitted(fit, newdata = nd, dpar = "sigma")
pr("hand sigma vs fitted(newdata, dpar = sigma): %.3e", max(abs(hsg(nd) - fs[, 1])))
nd2 <- nd; nd2$g2[1] <- "NEW"; nd2$g1[2] <- "NEW"; nd2$g1[3] <- "NEW2"; nd2$g2[3] <- "NEW3"
a <- tr(fitted(fit, newdata = nd2, dpar = "mu"))
pr("new level without allow_new_levels: %s", msg(a, 110))
a <- tr(fitted(fit, newdata = nd2, dpar = "mu", allow_new_levels = TRUE))
if (is_err(a)) pr("allow_new_levels %s", msg(a, 200)) else {
  pr("allow_new_levels row 1: fitted %.4f; old member alone at its weight %.4f; old member weight rescaled 1 %.4f; population %.4f", a[1, 1],
     cf["Intercept", 1] + cf["x", 1] * nd2$x[1] + nd2$w1[1] * re1[nd2$g1[1], 1] / (nd2$w1[1] + nd2$w2[1]),
     cf["Intercept", 1] + cf["x", 1] * nd2$x[1] + re1[nd2$g1[1], 1],
     cf["Intercept", 1] + cf["x", 1] * nd2$x[1])
}
for (snl in c("uncertainty", "gaussian", "old_levels")) {
  a <- tr(fitted(fit, newdata = nd2, dpar = "mu", allow_new_levels = TRUE, sample_new_levels = snl))
  pr("  fitted sample_new_levels = %s: %s %s", snl, msg(a, 130), if (is_err(a)) "" else paste(round(a[1:3, 1], 3), collapse = " "))
}
for (snl in c("uncertainty", "gaussian", "old_levels")) {
  a <- tr(predict(fit, newdata = nd2, allow_new_levels = TRUE, sample_new_levels = snl))
  pr("  predict sample_new_levels = %s: %s %s", snl, msg(a, 130), if (is_err(a)) "" else paste(round(a[1:3, 1], 3), collapse = " "))
}
a <- tr(fitted(fit, newdata = nd, re_formula = NA, dpar = "mu"))
pr("re_formula = NA vs fixed: %s %s", msg(a, 100), if (is_err(a)) "" else format(max(abs(a[, 1] - (cf["Intercept", 1] + cf["x", 1] * nd$x))), digits = 3))
a <- tr(fitted(fit, newdata = nd, re_formula = ~ (1 | p | mm(g1, g2, weights = cbind(w1, w2))), dpar = "mu"))
pr("re_formula the term written back: %s %s", msg(a, 170), if (is_err(a)) "" else format(max(abs(a[, 1] - ft[, 1])), digits = 3))
a <- tr(fitted(fit, newdata = nd, re_formula = ~ (1 | mm(g1, g2, weights = cbind(w1, w2))), dpar = "mu"))
pr("re_formula without the key: %s %s", msg(a, 170), if (is_err(a)) "" else format(max(abs(a[, 1] - ft[, 1])), digits = 3))
nd3 <- nd; nd3$w1 <- 1; nd3$w2 <- 0; nd3$v1 <- 0; nd3$v2 <- 1
a <- fitted(fit, newdata = nd3, dpar = "mu")
pr("newdata weights: mu vs hand %.3e", max(abs(a[, 1] - hmu(nd3))))
a <- fitted(fit, newdata = nd3, dpar = "sigma")
pr("newdata weights: sigma vs hand %.3e", max(abs(a[, 1] - hsg(nd3))))
a <- tr(simulate(fit, nsim = 1, seed = 1)); pr("simulate: %s", msg(a, 100))
a <- tr(frm_linpred(fit, newdata = nd)); pr("frm_linpred newdata: %s", msg(a, 100))
a <- tr(conditional_effects(fit, "x")); pr("conditional_effects: %s", msg(a, 100))
a <- tr(coef(fit)); pr("coef: %s", msg(a, 100))
a <- tr(VarCorr(fit)); pr("VarCorr: %s", msg(a, 100))
a <- tr(predict(fit, newdata = nd, re_formula = NA)); pr("predict re_formula = NA: %s", msg(a, 100))
pr("names: %s", paste(grep("^(sd|cor)_", frmtmb::variables(fit), value = TRUE), collapse = " "))

# mm cov = A, level order: hand density with A reordered
Ax <- rbind(cbind(A0, 0.1), 0.1); Ax[13, 13] <- 1
dimnames(Ax) <- list(c(lev, "EXTRA"), c(lev, "EXTRA"))
Ax <- Ax[c(13, 12:1), c(13, 12:1)]
fa <- frm(bf(y ~ x + (1 | mm(g1, g2, weights = cbind(w1, w2), cov = A))), data = d, data2 = list(A = Ax))
fb <- frm(bf(y ~ x + (1 | mm(g1, g2, weights = cbind(w1, w2), cov = A))), data = d, data2 = list(A = A0))
pr("A with an extra level, reordered vs plain: logLik %.8f vs %.8f", as.numeric(logLik(fa)), as.numeric(logLik(fb)))
fc <- tr(frm(bf(y ~ x + (1 | gr(g1, cov = A))), data = d, data2 = list(A = Ax)))
ra <- ranef(fa)[[1]]
pr("ranef levels: %s", paste(rownames(ra), collapse = ","))
a <- tr(fitted(fa, newdata = d[1:4, ])); pr("cov: fitted(newdata) vs own %s", if (is_err(a)) msg(a, 90) else format(max(abs(a[, 1] - fitted(fa)[1:4, 1])), digits = 3))
nd4 <- d[1:4, ]; nd4$g2[1] <- "EXTRA"; nd4$g1[2] <- "NEWZ"
a <- tr(fitted(fa, newdata = nd4)); pr("cov, a level A names but the data do not, no allow_new_levels: %s", msg(a, 130))
a <- tr(fitted(fa, newdata = nd4, allow_new_levels = TRUE))
pr("cov, new levels: %s %s", msg(a, 130), if (is_err(a)) "" else paste(round(a[1:2, 1], 3), collapse = " "))
b1 <- tr(fitted(frm(bf(y ~ x + (1 | gr(g1, cov = A))), data = d, data2 = list(A = A0)), newdata = transform(d[1:2, ], g1 = "NEWZ"), allow_new_levels = TRUE))
pr("control gr(cov): new levels %s", msg(b1, 100))

# a non-symmetric A: which triangle is read
Aas <- A0; Aas[1, 2] <- 0.1
fu <- frm(bf(y ~ x + (1 | mm(g1, g2, cov = A))), data = d, data2 = list(A = Aas))
Aas2 <- A0; Aas2[2, 1] <- 0.1
fl <- frm(bf(y ~ x + (1 | mm(g1, g2, cov = A))), data = d, data2 = list(A = Aas2))
f0 <- frm(bf(y ~ x + (1 | mm(g1, g2, cov = A))), data = d, data2 = list(A = A0))
pr("non-symmetric A: logLik with [1,2] changed %.5f, [2,1] changed %.5f, symmetric %.5f", as.numeric(logLik(fu)), as.numeric(logLik(fl)), as.numeric(logLik(f0)))
gu <- frm(bf(y ~ x + (1 | gr(g1, cov = A))), data = d, data2 = list(A = Aas))
pr("control gr(cov): %.5f (silently fits too)", as.numeric(logLik(gu)))

# prior tables
forms <- list(
  id = "bf(y ~ x + (1 + x | p | mm(g1, g2)), sigma ~ 1 + (1 | p | mm(g1, g2)))",
  idw = "bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w1, w2))), sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(v1, v2))))",
  cov = "bf(y ~ x + (1 + x | mm(g1, g2, cov = A)))",
  mmc = "bf(y ~ x + (1 + mmc(x1, x2) | p | mm(g1, g2)), sigma ~ 1 + (1 | p | mm(g1, g2)))",
  by = "bf(y ~ x + (1 | p | mm(g1, g2, by = cbind(f1, f2))), sigma ~ 1 + (1 | p | mm(g1, g2, by = cbind(f1, f2))))",
  idmv = "mvbf(bf(y ~ x + (1 | p | mm(g1, g2))), bf(w ~ x + (1 | p | mm(g1, g2))), rescor = FALSE)")
key <- function(t) {
  t <- as.data.frame(t)
  paste(t$class, t$coef, t$group, t$resp, t$dpar, t$nlpar, sep = "|")
}
for (nm in names(forms)) {
  a <- tr(default_prior(eval(parse(text = forms[[nm]]), envir = asNamespace("frmtmb")), data = d, data2 = list(A = A0)))
  b <- tr(brms::default_prior(eval(parse(text = forms[[nm]]), envir = asNamespace("brms")), data = d, data2 = list(A = A0)))
  if (is_err(a) || is_err(b)) { pr("table %s: frm %s brms %s", nm, msg(a), msg(b)); next }
  ka <- key(a); kb <- key(b)
  pr("table %s: same set %s same order %s; only frm [%s] only brms [%s]", nm, setequal(ka, kb), identical(ka, kb),
     paste(setdiff(ka, kb), collapse = "; "), paste(setdiff(kb, ka), collapse = "; "))
}
# frm_sample on id with different weights: which reading
pr("frm_sample block: id with the two weights uses the same core Z (reads re_blocks)")

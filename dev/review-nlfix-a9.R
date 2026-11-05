# Loose ends, one block each. Seed 9901.
source("dev/review-nlfix-common.R")
`%||%` <- function(a, b) if (is.null(a)) b else a
ll <- function(r) if (iserr(r)) paste("ERR", substr(r$msg, 1, 260)) else paste0(format(as.numeric(logLik(r)), digits = 10), " conv=", r$opt$convergence, " beta=", paste(format(r$estimates$beta, digits = 4), collapse = ","))
set.seed(9901)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 20)), h = factor(rep(1:8, length.out = n)), t = rep(1:20, 10))
cat("\n--- 1. the guard on a large response and a small covariate\n")
d$xs <- d$x * 1e-10
d$ybig <- 1e6 + 5e4 * d$x + rnorm(n, 0, 1e3)
cat("linear  :", ll(tt(frm(ybig ~ xs, data = d))), "\n")
cat("a + b*xs:", ll(tt(frm(bf(ybig ~ a + b * xs, a + b ~ 1, nl = TRUE), data = d))), "\n")
cat("b*xs only:", ll(tt(frm(bf(ybig ~ b * xs, b ~ 1, nl = TRUE), data = d))), "\n")
cat("linear no intercept:", ll(tt(frm(ybig ~ 0 + xs, data = d))), "\n")
d$yr <- 1 + 5e-3 * d$x + rnorm(n, 0, 1e-3); d$xm <- d$x * 1e-17
cat("a*xm at y ~ 1:", ll(tt(frm(bf(yr ~ a * xm, a ~ 1, nl = TRUE), data = d))), " linear:", ll(tt(frm(yr ~ 0 + xm, data = d))), "\n")
cat("\n--- 2. the sampler route\n")
suppressPackageStartupMessages(library(frmtmb.sample))
d$ys <- ifelse(d$x > 0.5, 1, 0) + rnorm(n, 0, 0.3)
d$x3 <- d$x + 3
r <- tt(frm_sample(bf(ys ~ inv_logit(10 * (x3 - c)), c ~ 1, nl = TRUE), data = d, chains = 1, iter = 100, refresh = 0, seed = 2, prior = set_prior("normal(0, 5)", nlpar = "c"))); show("frm_sample smooth step", if (iserr(r)) r else 1)
r <- tt(frm_sample(bf(ys ~ step(x3 - c), c ~ 1, nl = TRUE), data = d, chains = 1, iter = 100, refresh = 0, seed = 2, prior = set_prior("normal(0, 5)", nlpar = "c"))); show("frm_sample step, default inits", if (iserr(r)) r else 1)
r <- tt(frm_sample(bf(y2 ~ a + 2 * b - c, a + b + c ~ 1, nl = TRUE), data = transform(d, y2 = x), chains = 1, iter = 100, refresh = 0, seed = 2, prior = set_prior("normal(0, 1)", nlpar = "a"))); show("frm_sample a + 2b - c", if (iserr(r)) r else 1)
r <- tt(frm_sample(bf(y2 ~ 2 * a - b, a + b ~ 1, nl = TRUE), data = transform(d, y2 = x), chains = 1, iter = 100, refresh = 0, seed = 2, prior = c(set_prior("normal(0, 1)", nlpar = "a"), set_prior("normal(0, 1)", nlpar = "b")))); show("frm_sample 2 a - b", if (iserr(r)) r else 1)
cat("brms builds 2a - b:", nchar(as.character(brms::stancode(brms::bf(x ~ 2 * a - b, a + b ~ 1, nl = TRUE), data = d))) > 0, "\n")
cat("\n--- 3. conditional_effects(band = 'boot') on an ordinal body: named by the refusal\n")
d$yo <- cut(0.8 * d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = 1:4, ordered_result = TRUE)
a <- frm(bf(yo ~ b1 * x, b1 ~ 1, nl = TRUE), family = cumulative(), data = d)
print(names(formals(frmtmb:::conditional_effects.frmtmb_fit)))
x <- tt(system.time(ce <- conditional_effects(a, effects = "x", band = "boot"))); show("band boot", if (iserr(x)) x else x[["elapsed"]])
x <- tt(conditional_effects(a, effects = "x", dpar = "mu")); show("dpar mu", if (iserr(x)) x else nrow(x[[1]]))
cat("\n--- 4. partial re_formula on an ordinal body\n")
f2 <- frm(bf(yo ~ b1 * x + b2, b1 ~ 1 + (1 | g), b2 ~ 0 + (1 | h), nl = TRUE), family = cumulative(), data = d)
l2 <- frm(bf(yo ~ 0 + x + (0 + x | g) + (1 | h)), family = cumulative(), data = d)
cat("logLik:", as.numeric(logLik(f2)), as.numeric(logLik(l2)), "\n")
nd <- d[1:6, ]
for (rf in list(NULL, NA, ~ (0 + x | g), ~ (1 | h))) {
  lab <- if (is.null(rf)) "NULL" else if (is.na(rf)[1]) "NA" else deparse1(rf)
  rfn <- if (is.null(rf)) NULL else if (is.na(rf)[1]) NA else if (grepl("x", lab)) ~ (1 | g) else ~ (1 | h)
  a1 <- tt(fitted(f2, newdata = nd, re_formula = rfn)); b1 <- tt(fitted(l2, newdata = nd, re_formula = rf))
  cat(sprintf("%-16s fitted: %s\n", lab, if (iserr(a1)) paste("nl ERR", substr(a1$msg, 1, 200)) else if (iserr(b1)) paste("lin ERR", substr(b1$msg, 1, 200)) else format(max(abs(a1[, "Estimate", ] - b1[, "Estimate", ])), digits = 3)))
  a2 <- tt(frm_linpred(f2, newdata = nd, re_formula = rfn)); b2 <- tt(frm_linpred(l2, newdata = nd, re_formula = rf))
  cat(sprintf("%-16s linpred: %s\n", lab, if (iserr(a2)) paste("nl ERR", substr(a2$msg, 1, 200)) else if (iserr(b2)) paste("lin ERR", substr(b2$msg, 1, 200)) else format(max(abs(as.numeric(a2) - as.numeric(b2))), digits = 3)))
}
cat("\n--- 5. na.exclude, mi() and re_formula under a body with a cov = FALSE term\n")
de <- d; de$y <- 1 + 0.5 * de$x + ave(rnorm(n, 0, 0.5), de$g, FUN = function(z) as.numeric(stats::filter(z, 0.5, "recursive")))
st <- list(beta = c("a_(Intercept)" = 1, "b_(Intercept)" = 0.5))
dn <- de; dn$x[c(5, 50)] <- NA
fn <- tt(frm(bf(y ~ a + b * x, a + b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = FALSE)), data = dn, start = st, na.action = na.exclude))
fl <- tt(frm(bf(y ~ x, autocor = ~ ar(time = t, gr = g, cov = FALSE)), data = dn, na.action = na.exclude))
cat("logLik", ll(fn), "|", ll(fl), "\n")
if (!iserr(fn) && !iserr(fl)) {
  v <- tt(frm_linpred(fn)); w <- tt(frm_linpred(fl)); cat("linpred na.exclude: ", if (iserr(v)) v$msg else format(max(abs(v - w), na.rm = TRUE), digits = 3), " NA", sum(is.na(v)), sum(is.na(w)), "\n")
  v <- tt(fitted(fn)); w <- tt(fitted(fl)); cat("fitted na.exclude: ", if (iserr(v)) v$msg else format(max(abs(v[, 1] - w[, 1]), na.rm = TRUE), digits = 3), "\n")
  v <- tt(residuals(fn)); w <- tt(residuals(fl)); cat("residuals na.exclude: ", if (iserr(v)) v$msg else format(max(abs(v[, 1] - w[, 1]), na.rm = TRUE), digits = 3), "\n")
}
dm <- de; dm$y[c(7, 8, 60)] <- NA
fn <- tt(frm(bf(y | mi() ~ a + b * x, a + b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = FALSE)), data = dm, start = st))
fl <- tt(frm(bf(y | mi() ~ x, autocor = ~ ar(time = t, gr = g, cov = FALSE)), data = dm))
cat("mi() logLik", ll(fn), "|", ll(fl), "\n")
if (!iserr(fn) && !iserr(fl)) { v <- tt(fitted(fn)); w <- tt(fitted(fl)); cat("mi() fitted: ", if (iserr(v)) v$msg else format(max(abs(v[, 1] - w[, 1]), na.rm = TRUE), digits = 3), "\n") }
frn <- tt(frm(bf(y ~ a + b * x, a ~ 1 + (1 | h), b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = FALSE)), data = de, start = st))
frl <- tt(frm(bf(y ~ x + (1 | h), autocor = ~ ar(time = t, gr = g, cov = FALSE)), data = de))
cat("RE + cov=FALSE:", ll(frn), "|", ll(frl), "\n")
for (rf in list(NULL, NA)) {
  v <- tt(fitted(frn, re_formula = rf)); w <- tt(fitted(frl, re_formula = rf))
  cat(" re_formula", format(rf), ":", if (iserr(v)) v$msg else format(max(abs(v[, 1] - w[, 1])), digits = 3), "\n")
}
cat("\n--- 6. readers on per-location bodies\n")
set.seed(9902)
dc <- data.frame(x = rnorm(n)); dc$cat <- factor(sapply(seq_len(n), function(i) sample(c("a", "b", "c"), 1, prob = c(1, exp(0.8 * dc$x[i]), exp(-0.5 * dc$x[i])))))
pn <- tt(frm(bf(cat ~ 1) + nlf(mub ~ ab + bb * x) + nlf(muc ~ ac + bc * x) + lf(ab ~ 1) + lf(bb ~ 1) + lf(ac ~ 1) + lf(bc ~ 1), family = categorical(), data = dc))
pl <- tt(frm(cat ~ x, family = categorical(), data = dc))
cat("categorical per-location:", ll(pn), "|", ll(pl), "\n")
if (!iserr(pn)) {
  cat(" fitted diff:", max(abs(fitted(pn)[, "Estimate", ] - fitted(pl)[, "Estimate", ])), " linpred muc diff:", max(abs(frm_linpred(pn, dpar = "muc") - frm_linpred(pl, dpar = "muc"))), "\n")
  x <- tt(predict(pn, newdata = dc[1:4, ])); show(" predict", if (iserr(x)) x else dim(x))
  x <- tt(simulate(pn, nsim = 2)); show(" simulate", if (iserr(x)) x else length(x))
  print(as.data.frame(get_prior(pn))[, c("prior", "class", "coef", "dpar", "nlpar")][1:4, ])
}
cat("\n--- 7. lbeta and other base names on the tape\n")
e <- frmtmb:::ad_overload_env(baseenv(), quote(lbeta(a, b)))
print(ls(e, all.names = TRUE)[1:8])
f <- function(p) lbeta(p[1], p[2])
environment(f) <- e
tp <- tt(RTMB::MakeTape(f, c(2, 3)))
if (!iserr(tp)) { cat("lbeta(2,3) tape", tp(c(2, 3)), "base", base::lbeta(2, 3), " grad", tp$jacobian(c(2, 3)), " ref", digamma(2) - digamma(5), digamma(3) - digamma(5), "\n") } else cat(tp$msg, "\n")
cat("\n--- 8. compat table for the term x nl, and the extension rows\n")
for (ft in c("ar()", "ma()", "arma()", "cosy()", "unstr()", "fcor()")) print(as.data.frame(frm_compat(ft, "nl"))[, c("feature_a", "feature_b", "status")])
print(as.data.frame(frm_compat("nl", "poisson"))$status); print(as.data.frame(frm_compat("ar()", "poisson"))$status)
cat("\n--- 9. a - exp(b * x) under Gamma against other starts\n")
set.seed(7501)
n2 <- 300
dg <- data.frame(x = runif(n2, 0.1, 3))
dg$y <- rgamma(n2, shape = 4, rate = 4 / (2 * exp(0.3 * dg$x)))
best <- -Inf
for (i in 1:15) {
  s <- c("a_(Intercept)" = runif(1, 1, 8), "b_(Intercept)" = rnorm(1, 0, 0.5))
  v <- tryCatch(as.numeric(logLik(suppressWarnings(suppressMessages(frm(bf(y ~ a - exp(b * x), a + b ~ 1, nl = TRUE), family = Gamma(), data = dg, start = list(beta = s)))))), error = function(e) -Inf)
  if (is.finite(v) && v > best) best <- v
}
cat("best over 15 starts", best, " default start", as.numeric(logLik(suppressMessages(frm(bf(y ~ a - exp(b * x), a + b ~ 1, nl = TRUE), family = Gamma(), data = dg)))), "\n")
cat("\n--- 10. the sratio spec's better mode\n")
E <- new.env(parent = asNamespace("frmtmb")); E$`%||%` <- `%||%`
for (f in c("dev/equiv-specs.R")) sys.source(f, envir = E)
ds <- E$eq_dat("ordinal", seed = 4932)
set.seed(11)
bb <- -Inf; bs <- NULL
for (i in 1:30) {
  s <- c("b1_(Intercept)" = rnorm(1, 0, 3), "b2_(Intercept)" = rnorm(1, 0, 2), "b2_z" = rnorm(1, 0, 1))
  r <- tryCatch(suppressWarnings(suppressMessages(frm(bf(y ~ b1 * exp(b2 * x), b1 ~ 1, b2 ~ 1 + z, nl = TRUE), family = sratio(), data = ds, start = list(beta = s)))), error = function(e) NULL)
  if (!is.null(r) && as.numeric(logLik(r)) > bb) { bb <- as.numeric(logLik(r)); bs <- r$estimates$beta }
}
cat("best over 30 frm() starts:", bb, " at", format(bs, digits = 4), "\n")
cat("\n--- 11. frm_curve on a linear cov = FALSE, before and after arms read the same\n")
print(tt(as.data.frame(frm_curve(frl, newdata = de[1:3, ]))[1, ".estimate"])[1])

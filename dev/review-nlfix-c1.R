# Fix round 1, item 1 and 2: try to fool nl_fit_off_saddle() on bodies off the
# lane's table. Each row: the default-start fit, the best of random starts of frm().
# Seeds in the script.
source("dev/review-nlfix-common.R")
`%||%` <- function(a, b) if (is.null(a)) b else a
ns <- asNamespace("frmtmb")
msh <- function(z, k = 190) if (length(z)) substr(z[1], 1, k) else ""
best_ms <- function(f, fam, d, nms, k = 10, scale = 1, seed = 1) {
  set.seed(seed)
  b <- -Inf
  for (i in seq_len(k)) {
    st <- setNames(rnorm(length(nms), 0, scale), nms)
    r <- tryCatch(suppressWarnings(suppressMessages(frm(f, family = fam, data = d, start = list(beta = st)))), error = function(e) NULL)
    if (!is.null(r)) {
      v <- as.numeric(logLik(r))
      if (is.finite(v) && v > b && v < 1e6) b <- v
    }
  }
  b
}
one <- function(lab, f, fam, d, nms = NULL, ref = NULL, k = 10, ...) {
  t0 <- proc.time()[["elapsed"]]
  r <- tt(frm(f, family = fam, data = d, ...))
  el <- proc.time()[["elapsed"]] - t0
  if (is.null(ref) && !is.null(nms)) ref <- best_ms(f, fam, d, nms, k = k)
  if (iserr(r)) {
    cat(sprintf("%-34s STOP %s | ref=%s\n", lab, substr(r$msg, 1, 230), if (is.null(ref)) "NA" else format(ref, digits = 9)))
    return(invisible(r))
  }
  ll <- as.numeric(logLik(r))
  cat(sprintf("%-34s ll=%.5f conv=%s ref=%s gap=%s t=%.2f beta=%s\n    MSG[%s]\n    WARN[%s]\n", lab, ll, r$opt$convergence,
              if (is.null(ref)) "NA" else format(ref, digits = 9), if (is.null(ref)) "NA" else format(ref - ll, digits = 3), el,
              paste(format(r$estimates$beta, digits = 3), collapse = ","), msh(attr(r, "m")), msh(attr(r, "w"))))
  invisible(r)
}
set.seed(11001)
n <- 300
d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3), g = factor(rep(1:15, each = 20)), h = factor(rep(1:12, length.out = n)))
u <- rnorm(15, 0, 0.5)
lat <- 1.5 * exp(0.6 * d$x) + 0.7 * d$z
d$yo <- ordered(cut(lat + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
d$yn <- 1 + 0.7 * d$x + 0.6 * d$z + rnorm(n, 0, exp(-0.5 + 0.4 * d$z))
d$yo2 <- ordered(cut(1.2 * exp(0.5 * d$z) + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
d$yore <- ordered(cut(1.5 * exp(0.6 * d$x) + u[as.integer(d$g)] + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
d$ynre <- 1 + 0.7 * d$x + u[as.integer(d$g)] + rnorm(n, 0, 0.5)
d$ymix <- ifelse(runif(n) < 0.4, -2 + 0.8 * d$x, 2 - 0.6 * d$x) + rnorm(n, 0, 0.5)
d$w <- ifelse(seq_len(n) %% 20 == 0, d$x, 0)
d$yw <- 1 + 2 * d$w + rnorm(n, 0, 0.5)
d$yow <- ordered(cut(2 * d$w + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = 1:4))
cat("=== A. saddle with a nonconstant body\n")
one("gauss c0 + a x + b1 b2 z", bf(yn ~ c0 + a * x + b1 * b2 * z, c0 + a + b1 + b2 ~ 1, nl = TRUE), gaussian(), d, c("c0_(Intercept)", "a_(Intercept)", "b1_(Intercept)", "b2_(Intercept)"))
one("cumul a x + b1 b2 z", bf(yo ~ a * x + b1 * b2 * z, a + b1 + b2 ~ 1, nl = TRUE), cumulative(), d, c("a_(Intercept)", "b1_(Intercept)", "b2_(Intercept)"))
one("cumul a z + b1 exp(b2 x)", bf(yo ~ a * z + b1 * exp(b2 * x), a + b1 + b2 ~ 1, nl = TRUE), cumulative(), d, c("a_(Intercept)", "b1_(Intercept)", "b2_(Intercept)"))
one("cumul b1 exp(b2 x) + 0.5 z", bf(yo ~ b1 * exp(b2 * x) + 0.5 * z, b1 + b2 ~ 1, nl = TRUE), cumulative(), d, c("b1_(Intercept)", "b2_(Intercept)"))
one("cumul b1 exp(b2 x) + 1e-6 z", bf(yo ~ b1 * exp(b2 * x) + 1e-6 * z, b1 + b2 ~ 1, nl = TRUE), cumulative(), d, c("b1_(Intercept)", "b2_(Intercept)"))
one("cumul same, se = FALSE", bf(yo ~ b1 * exp(b2 * x) + 0.5 * z, b1 + b2 ~ 1, nl = TRUE), cumulative(), d, se = FALSE)
cat("\n=== B. body constant on most rows\n")
one("gauss a + b1 b2 w, 5% rows", bf(yw ~ a + b1 * b2 * w, a + b1 + b2 ~ 1, nl = TRUE), gaussian(), d, c("a_(Intercept)", "b1_(Intercept)", "b2_(Intercept)"))
one("cumul b1 exp(b2 w), 5% rows", bf(yow ~ b1 * exp(b2 * w), b1 + b2 ~ 1, nl = TRUE), cumulative(), d, c("b1_(Intercept)", "b2_(Intercept)"))
cat("\n=== C. several nlpars with random effects\n")
one("cumul b1 exp(b2 x), RE b1", bf(yore ~ b1 * exp(b2 * x), b1 ~ 1 + (1 | g), b2 ~ 1, nl = TRUE), cumulative(), d, c("b1_(Intercept)", "b2_(Intercept)"), k = 4)
one("cumul RE b1, RE b2", bf(yore ~ b1 * exp(b2 * x), b1 ~ 1 + (1 | g), b2 ~ 1 + (1 | h), nl = TRUE), cumulative(), d, c("b1_(Intercept)", "b2_(Intercept)"), k = 4)
one("gauss a + b1 b2 x, RE a, b1", bf(ynre ~ a + b1 * b2 * x, a ~ 1 + (1 | g), b1 ~ 1 + (1 | h), b2 ~ 1, nl = TRUE), gaussian(), d, c("a_(Intercept)", "b1_(Intercept)", "b2_(Intercept)"), k = 4)
cat("\n=== D. mixtures\n")
one("mix(cum, cum) shared saddle", bf(yo ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), mixture(cumulative(), cumulative()), d, c("b1_(Intercept)", "b2_(Intercept)"), k = 4)
fm <- bf(ymix ~ 1) + nlf(mu1 ~ a1 + b1 * c1 * x) + nlf(mu2 ~ a2 + b2 * c2 * x) + lf(a1 ~ 1) + lf(b1 ~ 1) + lf(c1 ~ 1) + lf(a2 ~ 1) + lf(b2 ~ 1) + lf(c2 ~ 1)
r <- one("mix per-location two saddles", fm, mixture(gaussian(), gaussian()), d, c("a1_(Intercept)", "b1_(Intercept)", "c1_(Intercept)", "a2_(Intercept)", "b2_(Intercept)", "c2_(Intercept)"), k = 8)
cat("\n=== E. two bodies, each on its saddle\n")
fm2 <- mvbf(bf(yo ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE, family = cumulative()), bf(yo2 ~ c1 * exp(c2 * z), c1 + c2 ~ 1, nl = TRUE, family = cumulative()), rescor = FALSE)
r1 <- tt(frm(bf(yo ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), family = cumulative(), data = d))
r2 <- tt(frm(bf(yo2 ~ c1 * exp(c2 * z), c1 + c2 ~ 1, nl = TRUE), family = cumulative(), data = d))
cat("separately:", as.numeric(logLik(r1)), as.numeric(logLik(r2)), " sum", as.numeric(logLik(r1)) + as.numeric(logLik(r2)), "\n")
one("mv cumulative, two saddles", fm2, NULL, d, ref = as.numeric(logLik(r1)) + as.numeric(logLik(r2)))
fm3 <- bf(yn ~ a + b1 * c1 * x) + nlf(sigma ~ exp(s0 + s1 * s2 * z)) + lf(a ~ 1) + lf(b1 ~ 1) + lf(c1 ~ 1) + lf(s0 ~ 1) + lf(s1 ~ 1) + lf(s2 ~ 1)
fm3b <- bf(yn ~ 1) + nlf(mu ~ a + b1 * c1 * x) + nlf(sigma ~ exp(s0 + s1 * s2 * z)) + lf(a ~ 1) + lf(b1 ~ 1) + lf(c1 ~ 1) + lf(s0 ~ 1) + lf(s1 ~ 1) + lf(s2 ~ 1)
one("gauss mu and sigma, two saddles", fm3b, gaussian(), d, c("a_(Intercept)", "b1_(Intercept)", "c1_(Intercept)", "s0_(Intercept)", "s1_(Intercept)", "s2_(Intercept)"), k = 8)
cat("linear-spelled limit, mu ~ x, sigma ~ z:", as.numeric(logLik(frm(bf(yn ~ x, sigma ~ z), data = d))), "\n")
fm4 <- bf(yo ~ 1) + nlf(mu ~ b1 * exp(b2 * x)) + nlf(disc ~ exp(d1 * d2 * z)) + lf(b1 ~ 1) + lf(b2 ~ 1) + lf(d1 ~ 1) + lf(d2 ~ 1)
one("cumul mu and disc, two saddles", fm4, cumulative(), d, c("b1_(Intercept)", "b2_(Intercept)", "d1_(Intercept)", "d2_(Intercept)"), k = 8)
cat("\n=== F. a slope of exactly zero in the data: no saddle to escape\n")
dz <- d
dz$y0 <- residuals(lm(yn ~ x, data = d)) + 1
one("gauss a + b1 b2 x, slope 0", bf(y0 ~ a + b1 * b2 * x, a + b1 + b2 ~ 1, nl = TRUE), gaussian(), dz)
dz$y1 <- dz$y0 + 2e-5 * dz$x
one("gauss, slope 2e-5", bf(y1 ~ a + b1 * b2 * x, a + b1 + b2 ~ 1, nl = TRUE), gaussian(), dz)
cat(" linear fit, slope 0:", as.numeric(logLik(frm(y0 ~ x, data = dz))), "\n")
cat(" share of null data sets that stop, gaussian, n = 40, y independent of x, 400 sets:\n")
set.seed(11002)
stops <- 0L; tot <- 0L; firstmsg <- NULL
for (i in 1:400) {
  dn <- data.frame(x = rnorm(40), y = rnorm(40))
  r <- tryCatch(suppressWarnings(suppressMessages(frm(bf(y ~ a + b1 * b2 * x, a + b1 + b2 ~ 1, nl = TRUE), data = dn))), error = function(e) e)
  tot <- tot + 1L
  if (inherits(r, "error")) { stops <- stops + 1L; firstmsg <- conditionMessage(r) }
}
cat(" stops", stops, "of", tot, "\n"); if (!is.null(firstmsg)) cat(" ", substr(firstmsg, 1, 200), "\n")
cat("\n=== G. item 2(a): a strict local maximum that is not the maximum\n")
set.seed(11003)
dq <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
dq$yo <- ordered(cut(-0.15 * dq$x + 0.8 * dq$z + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = 1:4))
dq$yn <- -0.15 * dq$x + 0.8 * dq$z + rnorm(n)
one("cumul b^2 x + b^4 z", bf(yo ~ b^2 * x + b^4 * z, b ~ 1, nl = TRUE), cumulative(), dq, "b_(Intercept)")
one("gauss a + b^2 x + b^4 z", bf(yn ~ a + b^2 * x + b^4 * z, a + b ~ 1, nl = TRUE), gaussian(), dq, c("a_(Intercept)", "b_(Intercept)"))
one("cumul b^2 x (slope negative)", bf(yo ~ b^2 * x, b ~ 1, nl = TRUE), cumulative(), dq, "b_(Intercept)")
one("cumul b1^2 x + b2^2 z (mixed)", bf(yo ~ b1^2 * x + b2^2 * z, b1 + b2 ~ 1, nl = TRUE), cumulative(), dq, c("b1_(Intercept)", "b2_(Intercept)"))
cat("\n=== H. item 2(c): a true maximum on nearly noiseless data\n")
for (sg in c(1e-3, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9)) {
  set.seed(11004)
  dn <- data.frame(x = rnorm(200)); dn$y <- 1 + 2 * dn$x + rnorm(200, 0, sg)
  lin <- as.numeric(logLik(frm(y ~ x, data = dn)))
  r <- tt(frm(bf(y ~ a + b1 * b2 * x, a + b1 + b2 ~ 1, nl = TRUE), data = dn))
  r2 <- tt(frm(bf(y ~ a + b1 * b2 * x, a + b1 + b2 ~ 1, nl = TRUE), data = dn, start = list(beta = c("a_(Intercept)" = 1, "b1_(Intercept)" = 1, "b2_(Intercept)" = 1))))
  cat(sprintf("sigma %.0e: linear %.3f | body default %s | body started %s\n", sg, lin, if (iserr(r)) paste("STOP", substr(r$msg, 1, 100)) else format(as.numeric(logLik(r)), digits = 8), if (iserr(r2)) "ERR" else format(as.numeric(logLik(r2)), digits = 8)))
  # nl_restart_credible() was removed in fix round 2
}
cat("\n=== I. control and the flag\n")
set.seed(9202)
dk <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
dk$y <- ordered(cut(1.5 * exp(0.6 * dk$x) + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
fk <- bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)
one("acat default", fk, acat(), dk)
one("acat check_se = ignore", fk, acat(), dk, control = frmtmb_control(check_se = "ignore"))
one("acat se = FALSE", fk, acat(), dk, se = FALSE)
dg <- transform(dk, yg = 1 + 0.5 * x + rnorm(n))
fg <- bf(yg ~ a * exp(b * x), a + b ~ 1, nl = TRUE)
st0 <- as.numeric(logLik(tt(frm(fg, data = dg, start = list(beta = c("a_(Intercept)" = 1, "b_(Intercept)" = 0.5))))))
cat("ordinary fit:", st0, "\n")
one("iter.max = 0 (the objective at the start)", bf(yg ~ a + b * x, a + b ~ 1, nl = TRUE), gaussian(), dg, control = frmtmb_control(optCtrl = list(iter.max = 0L, eval.max = 1L)))
one("iter.max = 2", bf(yg ~ a * exp(b * x), a + b ~ 1, nl = TRUE), gaussian(), dg, control = frmtmb_control(optCtrl = list(iter.max = 2L)))
one("iter.max = 2 linear spelling", yg ~ x, gaussian(), dg, control = frmtmb_control(optCtrl = list(iter.max = 2L)))
cat("\n=== J. RNG, order, and what stays with the fit\n")
set.seed(1); s0 <- .Random.seed
a1 <- suppressWarnings(suppressMessages(frm(fk, family = cumulative(), data = dk)))
cat("RNG state unchanged by a restarted fit:", identical(s0, .Random.seed), "\n")
set.seed(999); a2 <- suppressWarnings(suppressMessages(frm(fk, family = cumulative(), data = dk)))
cat("another seed: identical estimates:", identical(a1$estimates, a2$estimates), "\n")
perm <- rev(seq_len(n))
a3 <- suppressWarnings(suppressMessages(frm(fk, family = cumulative(), data = dk[perm, ])))
cat("rows reversed: ll", format(c(as.numeric(logLik(a1)), as.numeric(logLik(a3))), digits = 12), " beta", format(a3$estimates$beta, digits = 6), "\n")
u1 <- tt(update(a1)); cat("update():", if (iserr(u1)) u1$msg else as.numeric(logLik(u1)), " msg:", msh(attr(u1, "m"), 90), "\n")
u2 <- tt(update(a1, data = dk[1:250, ])); cat("update(data):", if (iserr(u2)) u2$msg else as.numeric(logLik(u2)), " msg:", msh(attr(u2, "m"), 90), "\n")
u3 <- tt(as.numeric(logLik(frm(fk, family = cumulative(), data = dk[1:250, ], start = list(beta = c("b1_(Intercept)" = 1)))))); cat(" that subset from a start:", u3[1], "\n")
b <- tt(frm_bootstrap(a1, nsim = 12, seed = 1))
if (iserr(b)) cat("frm_bootstrap:", b$msg, "\n") else { cat("frm_bootstrap: class", class(b)[1], "\n"); bt <- b$t %||% b; print(round(apply(as.matrix(bt), 2, function(v) c(mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE), na = sum(is.na(v)), zero = sum(v == 0, na.rm = TRUE))), 3)[, 1:min(2, NCOL(bt)), drop = FALSE]) }
ce <- tt(conditional_effects(a1, effects = "x", band = "boot", boot = 6, seed = 1)); cat("conditional_effects boot:", if (iserr(ce)) ce$msg else nrow(ce[[1]]), "\n")
v <- tt(confint(a1, parm = "b2_Intercept", method = "profile")); cat("profile b2:", if (iserr(v)) substr(v$msg, 1, 200) else format(as.numeric(v), digits = 5), "\n")
print(a1$call)
cat("\n=== K. fits that are already at a maximum, with a flat direction: does a restart change them?\n")
for (lab in c("a + b", "b0 + b1 x, cumul", "shared mixture")) {
  t1 <- proc.time()[["elapsed"]]
  r <- switch(lab,
    "a + b" = tt(frm(bf(yg ~ a + b, a + b ~ 1, nl = TRUE), data = dg)),
    "b0 + b1 x, cumul" = tt(frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family = cumulative(), data = dk)),
    "shared mixture" = tt(frm(bf(yg ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family = mixture(gaussian(), gaussian()), data = dg)))
  el <- proc.time()[["elapsed"]] - t1
  cat(sprintf("%-20s ll=%s t=%.2f restart message: %s\n", lab, if (iserr(r)) r$msg else format(as.numeric(logLik(r)), digits = 10), el, msh(attr(r, "m"), 120)))
}

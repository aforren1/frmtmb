# Review of lane readers, short check of the fix round. Lane build.
# A: the divergence rule's reach: families where a draw enters a scale or
#    shape (seeds 7201..7207); B: the 11-node value and its 9-node warning
#    against an independent Gauss-Hermite rule (data seed 7101, inf.R case
#    1); C: the closed form against Monte Carlo (seed 71, 1e7 draws) and
#    where it must not apply (count of me_gh_moments() calls); D: the
#    dropped-draw shortcut on a periodic body (data seed 7208); E: the rule
#    check on the review's synthetic cases; F: m7, a body's Est.Error under
#    allow_new_levels against Monte Carlo (data seed 8201, MC seed 73, 2e5).
source("dev/review-readers-setup.R")
dm <- function(f, nd, dp) as.numeric(frm_linpred(f, newdata = nd, dpar = dp))
slope <- function(f, nd, dp, v) {
  nd1 <- nd
  nd1[[v]] <- nd1[[v]] + 1
  dm(f, nd1, dp) - dm(f, nd, dp)
}
fw <- function(f, nd, ...) {
  r <- with_warn(fitted(f, newdata = nd, ...))
  list(f = r$value, w = wclass(r$warnings),
       msg = vapply(r$warnings, conditionMessage, ""))
}
tryfit <- function(...) tryCatch(frm(...), error = function(e) {
  cat("  fit error:", conditionMessage(e), "\n")
  NULL
})
show <- function(lbl, r) {
  cat(sprintf("%-38s Est %s Err %s warn %s\n", lbl,
              paste(format(r$f[, "Estimate"], digits = 6), collapse = " "),
              paste(format(r$f[, "Est.Error"], digits = 4), collapse = " "),
              paste(unique(r$w), collapse = ",")))
}
n <- 800
mk <- function(seed) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), l = rnorm(n), sz = 0.3)
  d$zo <- d$l + rnorm(n, 0, 0.3)
  d
}

cat("== A. which means diverge, and which are flagged\n")
d <- mk(7201)
comp <- runif(n) < 0.4
d$y <- ifelse(comp, exp(1.5 + rnorm(n, 0, exp(-0.7 + 0.5 * d$l))),
              exp(-0.5 + rnorm(n, 0, 0.4)))
fm <- tryfit(bf(y ~ 1, sigma1 ~ me(zo, sz)),
             family = mixture(lognormal(), lognormal()), data = d)
nd <- data.frame(zo = 1, sz = c(0.3, 0.8))
if (!is.null(fm)) {
  show("mixture(lognormal x2), sigma1 ~ me()", fw(fm, nd))
  es <- dm(fm, nd, "sigma1"); bs <- slope(fm, nd[1, ], "sigma1", "zo")
  em <- dm(fm, nd, "mu1"); th <- as.numeric(frm_linpred(fm, newdata = nd,
                                                        dpar = "theta1",
                                                        type = "response"))
  cat("  sigma1 link sd", bs * nd$sz, " brms-style 20 x 4000 at sz 0.8:",
      format(range(vapply(1:20, function(s) {
        set.seed(s)
        mean(exp(em[2] + exp(2 * (es[2] + bs * 0.8 * rnorm(4000))) / 2))
      }, 0)), digits = 4), "(component 1 mean only)\n")
}
d <- mk(7202)
d$y <- exp(0.5 + 0.4 * d$x + rnorm(n, 0, exp(-0.5 + 0.5 * d$l)))
fs <- tryfit(bf(y ~ x, nlf(sigma ~ exp(s1)), s1 ~ me(zo, sz)),
             family = lognormal(), data = d)
if (!is.null(fs)) {
  show("lognormal, nlf(sigma ~ exp(s1)), s1 ~ me", fw(fs, data.frame(
    x = 0, zo = 1, sz = c(0.3, 0.8))))
}
fh <- NULL
d <- mk(7203)
yy <- exp(0.5 + 0.4 * d$x + rnorm(n, 0, exp(-0.5 + 0.5 * d$l)))
d$yh <- ifelse(runif(n) < 0.3, 0, yy)
fh <- tryfit(bf(yh ~ x, sigma ~ me(zo, sz)), family = hurdle_lognormal(),
             data = d)
if (!is.null(fh)) {
  show("hurdle_lognormal, sigma ~ me()", fw(fh, data.frame(
    x = 0, zo = 1, sz = c(0.3, 0.8))))
}
d <- mk(7204)
d$ys <- 1 + 0.5 * d$x + rt(n, exp(1 + 0.5 * d$l))
ft <- tryfit(bf(ys ~ x, nu ~ me(zo, sz)), family = student(), data = d)
nd <- data.frame(x = 0.5, zo = 1, sz = c(0.3, 3))
if (!is.null(ft)) {
  r <- fw(ft, nd)
  show("student, nu ~ me()", r)
  cat("  Est / mu - 1:", format(r$f[, "Estimate"] / dm(ft, nd, "mu") - 1,
                                digits = 3), "\n")
}
d <- mk(7205)
d$yw <- rweibull(n, exp(0.7 + 0.3 * d$l), exp(0.5 + 0.3 * d$x))
fwb <- tryfit(bf(yw ~ x, shape ~ me(zo, sz)), family = weibull(), data = d)
if (!is.null(fwb)) {
  r <- fw(fwb, nd)
  show("weibull, shape ~ me()", r)
  cat("  Est / exp(mu) - 1:",
      format(r$f[, "Estimate"] / exp(dm(fwb, nd, "mu")) - 1, digits = 3),
      "\n")
}
d <- mk(7206)
d$ysn <- 1 + 0.5 * d$x + rnorm(n) + 0.3 * abs(rnorm(n)) * d$l
fsn <- tryfit(bf(ysn ~ x, alpha ~ me(zo, sz)), family = skew_normal(),
              data = d)
if (!is.null(fsn)) {
  r <- fw(fsn, nd)
  show("skew_normal, alpha ~ me()", r)
  cat("  Est / mu - 1:", format(r$f[, "Estimate"] / dm(fsn, nd, "mu") - 1,
                                digits = 3), "\n")
}
d <- mk(7207)
d$yg <- 1 + 0.5 * d$x + 0.5 * (-log(-log(runif(n))))
fg <- tryfit(bf(yg ~ x, xi ~ me(zo, sz)), family = gen_extreme_value(),
             data = d)
if (!is.null(fg)) {
  r <- fw(fg, data.frame(x = 0.5, zo = 1, sz = c(0.3, 3)))
  show("gen_extreme_value, xi ~ me()", r)
  xi <- dm(fg, data.frame(x = 0.5, zo = 1, sz = 0.3), "xi")
  bx <- slope(fg, data.frame(x = 0.5, zo = 1, sz = 0.3), "xi", "zo")
  cat("  xi at the row", xi, " slope", bx, " P(xi >= 1) at sz 0.3 / 3:",
      pnorm((xi - 1) / abs(bx * c(0.3, 3))), "\n")
}

cat("\n== B. the 11-node value and the 9-node figure, inf.R case 1\n")
d <- mk(7101)
set.seed(7101)
d <- data.frame(x = rnorm(n), l = rnorm(n), sz = 0.3)
d$zo <- d$l + rnorm(n, 0, 0.3)
d$y <- exp(0.5 + 0.4 * d$x + rnorm(n, 0, exp(-0.5 + 0.5 * d$l)))
f1 <- frm(bf(y ~ x, sigma ~ me(zo, sz)), family = lognormal(), data = d)
gh <- function(K) {
  # Golub-Welsch for the probabilists' Hermite rule, independent of the lane
  J <- matrix(0, K, K)
  for (i in seq_len(K - 1)) J[i, i + 1] <- J[i + 1, i] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = e$vectors[1, ]^2)
}
for (sz in c(0.05, 0.3, 0.5, 0.8)) {
  nd <- data.frame(x = 0, zo = 1, sz = sz)
  em <- dm(f1, nd, "mu"); es <- dm(f1, nd, "sigma")
  bs <- slope(f1, nd, "sigma", "zo")
  m <- function(K) {
    g <- gh(K)
    sum(g$w * exp(em + exp(2 * (es + bs * sz * g$x)) / 2))
  }
  r <- fw(f1, nd)
  num <- regmatches(r$msg, regexpr("move that mean by [^ ]+", r$msg))
  cat(sprintf(paste("sz %.2f: fitted %.6g | own 11-node %.6g, 9-node %.6g,",
                    "|9/11 - 1| %.2g | warning: '%s' | 41-node %.4g\n"),
              sz, r$f[, "Estimate"], m(11), m(9), abs(m(9) / m(11) - 1),
              paste(num, collapse = ""), m(41)))
}

cat("\n== C. closed form against Monte Carlo, seed 71, 1e7 draws\n")
cnt <- new.env()
cnt$n <- 0L
suppressMessages(trace("me_gh_moments", where = ns, print = FALSE,
                       tracer = quote(cnt$n <- cnt$n + 1L)))
d <- mk(7209)
d$yp <- rpois(n, exp(0.2 + 1 * d$l))
d$yl <- exp(0.5 + 1 * d$l + rnorm(n, 0, 0.5))
lam <- exp(0.4 + 1 * d$l)
cc <- rpois(n, lam)
while (any(cc == 0)) cc[cc == 0] <- rpois(sum(cc == 0), lam[cc == 0])
d$yh <- ifelse(runif(n) < 0.3, 0L, cc)
d$yb <- rbinom(n, 1, plogis(0.2 + d$l))
fits <- list(
  poisson = frm(bf(yp ~ me(zo, sz)), family = poisson(), data = d),
  lognormal_mu = frm(bf(yl ~ me(zo, sz)), family = lognormal(), data = d),
  hurdle_mu = frm(bf(yh ~ me(zo, sz)), family = hurdle_poisson(), data = d),
  bernoulli = frm(bf(yb ~ me(zo, sz)), family = bernoulli(), data = d))
N <- 1e7
for (nm in names(fits)) {
  f <- fits[[nm]]
  b <- slope(f, data.frame(zo = 0, sz = 1), "mu", "zo")
  for (s in c(0.5, 1, 2, 3)) {
    nd <- data.frame(zo = 0.3, sz = s / b)
    cnt$n <- 0L
    mm <- ns$fitted_me_moments(f, nd, NULL, NULL, NULL, FALSE)
    calls <- cnt$n
    e <- dm(f, nd, "mu")
    set.seed(71)
    u <- e + s * rnorm(N)
    v <- switch(nm,
      poisson = exp(u),
      lognormal_mu = exp(u + as.numeric(frm_linpred(
        f, newdata = nd, dpar = "sigma", type = "response"))^2 / 2),
      hurdle_mu = {
        hu <- as.numeric(frm_linpred(f, newdata = nd, dpar = "hu",
                                     type = "response"))
        l <- exp(u)
        (1 - hu) * l / (-expm1(-l))
      },
      bernoulli = plogis(u))
    cat(sprintf(paste("%-13s sd %.1f: GH calls %2d | mean %.6g MC %.6g",
                      "z %6.2f | sd %.6g MC sd %.6g ratio %.4f\n"),
                nm, s, calls, mm$mean, mean(v),
                (mm$mean - mean(v)) / (sd(v) / sqrt(N)), sqrt(mm$var), sd(v),
                sqrt(mm$var) / sd(v)))
  }
}
untrace("me_gh_moments", where = ns)

cat("\n== D. a draw whose unit step moves nothing: y ~ sin(a)^2 + c\n")
set.seed(7208)
d <- data.frame(l1 = rnorm(n), l2 = rnorm(n), s1 = 0.2, s2 = 0.2)
d$x1 <- d$l1 + rnorm(n, 0, 0.2)
d$x2 <- d$l2 + rnorm(n, 0, 0.2)
d$y <- sin(0.3 + 0.8 * d$l1)^2 + (0.5 + 0.5 * d$l2) + rnorm(n, 0, 0.1)
fp <- tryfit(bf(y ~ sin(a)^2 + c, a ~ me(x1, s1), c ~ me(x2, s2),
                nl = TRUE), family = gaussian(), data = d)
if (!is.null(fp)) {
  ba <- slope(fp, data.frame(x1 = 0, x2 = 0, s1 = 1, s2 = 1), "a", "x1")
  a0 <- dm(fp, data.frame(x1 = 0, x2 = 0, s1 = 1, s2 = 1), "a")
  # a row where a sits at 0 and the draw's sd is pi: sin(+-pi)^2 = sin(0)^2
  x1 <- -a0 / ba
  nd <- data.frame(x1 = x1, x2 = 0, s1 = pi / abs(ba), s2 = 0.5)
  r <- fw(fp, nd)
  c0 <- dm(fp, nd, "c")
  bc <- slope(fp, nd, "c", "x2")
  set.seed(74)
  z <- rnorm(2e6)
  ref <- mean(sin(dm(fp, nd, "a") + pi * z)^2) + c0
  cat("a at the row", dm(fp, nd, "a"), " fitted", r$f[, "Estimate"],
      " MC", ref, " exact", (1 - exp(-2 * pi^2)) / 2 + c0,
      " value with a's draw at 0", sin(dm(fp, nd, "a"))^2 + c0, "\n")
}

cat("\n== E. me_rule_check() on the review's synthetic cases\n")
for (r in 3:4) {
  for (s in c(1, 1.5, 2, 2.5, 3)) {
    K <- ns$me_gh_size[r]
    val <- function(z) matrix(exp(s * z[[1L]] + 0.5 *
                                    Reduce(`+`, z[-1L], 0)), 1)
    mo <- ns$me_gh_moments(val, r, K)
    v <- s^2 + 0.25 * (r - 1)
    err <- sqrt(mo$m2 - mo$m1^2) / sqrt(exp(2 * v) - exp(v)) - 1
    w <- with_warn(ns$me_rule_check(val, r, K, mo, "y"))
    cat(sprintf("r %d K %2d s %.1f: sd error %9.2e | warns %s\n", r, K, s,
                err, if (length(w$warnings)) "yes" else "no"))
  }
}

cat("\n== F. m7: a body under allow_new_levels = TRUE\n")
set.seed(8201)
n2 <- 300; J <- 12
d <- data.frame(x = runif(n2, 0, 3), g = factor(sample(J, n2, TRUE)))
d$y <- (3 + rnorm(J, 0, 0.4)[d$g]) * exp(-0.8 * d$x) + rnorm(n2, 0, 0.1)
fA0 <- frm(bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE),
           family = gaussian(), data = d)
fA <- frm(bf(y ~ a * exp(-b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE),
          family = gaussian(), data = d)
nd <- data.frame(x = c(0.5, 2), g = d$g[1])
cat("no group term, flag / no flag Est.Error:",
    fitted(fA0, newdata = nd, allow_new_levels = TRUE)[, "Est.Error"], "/",
    fitted(fA0, newdata = nd)[, "Est.Error"], "\n")
cat("seen level, flag / no flag:",
    fitted(fA, newdata = nd, allow_new_levels = TRUE)[, "Est.Error"], "/",
    fitted(fA, newdata = nd)[, "Est.Error"], "\n")
ndn <- transform(nd, g = "new")
fn <- with_warn(fitted(fA, newdata = ndn, allow_new_levels = TRUE))
fe <- fixef(fA)[, "Estimate"]
V <- vcov(fA)[1:2, 1:2]
sg <- VarCorr(fA)$g$sd[1, 1]
set.seed(73)
B <- MASS::mvrnorm(2e5, fe[1:2], V)
u <- rnorm(2e5, 0, sg)
sims <- vapply(ndn$x, function(x) (B[, 1] + u) * exp(-B[, 2] * x),
               numeric(2e5))
cat("new level: Est", fn$value[, "Estimate"], " MC mean", colMeans(sims),
    "\n  Est.Error", fn$value[, "Est.Error"], " MC sd", apply(sims, 2, sd),
    " ratio", fn$value[, "Est.Error"] / apply(sims, 2, sd),
    "\n  Q", fn$value[, "Q2.5"], fn$value[, "Q97.5"], " MC q",
    apply(sims, 2, quantile, 0.025), apply(sims, 2, quantile, 0.975),
    "\n  warnings:", wclass(fn$warnings), "\n")
cat("frm_linpred(se.fit) at the new level:\n")
print(tryCatch(frm_linpred(fA, newdata = ndn, allow_new_levels = TRUE,
                           se.fit = TRUE),
               error = function(e) paste(class(e)[1], conditionMessage(e))))

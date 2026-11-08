# Review of lane readers, confirmation of the S1 fix. Lane build.
# A: mixture(lognormal, gaussian), sigma1 ~ me() on the lognormal part
#    (data seed 7301; the gaussian part at 8, so every response is
#    positive); B: hurdle_lognormal with a draw on hu only (seed
#    7302); C: Gamma with a draw on shape (seed 7303); D: the warning on two
#    separate calls; E: the sin(a)^2 + c row of dev/review-readers-sc.R D
#    (data seed 7208, MC seed 74, 2e6 draws).
source("dev/review-readers-setup.R")
dm <- function(f, nd, dp) as.numeric(frm_linpred(f, newdata = nd, dpar = dp))
slope <- function(f, nd, dp, v) {
  nd1 <- nd
  nd1[[v]] <- nd1[[v]] + 1
  dm(f, nd1, dp) - dm(f, nd, dp)
}
fw <- function(f, nd, ...) {
  r <- with_warn(fitted(f, newdata = nd, ...))
  list(f = r$value, w = wclass(r$warnings))
}
show <- function(lbl, r) {
  cat(sprintf("%-44s Est %s Err %s warn %s\n", lbl,
              paste(format(r$f[, "Estimate"], digits = 6), collapse = " "),
              paste(format(r$f[, "Est.Error"], digits = 4), collapse = " "),
              paste(r$w, collapse = ",")))
}
n <- 800
mk <- function(seed) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), l = rnorm(n), sz = 0.3)
  d$zo <- d$l + rnorm(n, 0, 0.3)
  d
}

cat("== A. mixture(lognormal, gaussian), sigma1 ~ me()\n")
d <- mk(7301)
comp <- runif(n) < 0.5
d$y <- ifelse(comp, exp(1.5 + rnorm(n, 0, exp(-0.7 + 0.5 * d$l))),
              8 + rnorm(n, 0, 0.5))
fa <- tryCatch(frm(bf(y ~ 1, sigma1 ~ me(zo, sz)),
                   family = mixture(lognormal(), gaussian()), data = d),
               error = function(e) { cat("fit error:", conditionMessage(e),
                                         "\n"); NULL })
if (!is.null(fa)) {
  cat("family:", fa$spec$responses[[1]]$family$family, " dpars:",
      names(fa$spec$responses[[1]]$dpars), "\n")
  bs <- slope(fa, data.frame(zo = 0, sz = 1), "sigma1", "zo")
  for (sz in c(0.3, 3, 10)) {
    show(sprintf("sz %.1f (sigma1 link sd %.3f)", sz, bs * sz),
         fw(fa, data.frame(zo = 1, sz = sz)))
  }
}

cat("\n== B. hurdle_lognormal, hu ~ me() only\n")
d <- mk(7302)
yy <- exp(0.5 + 0.4 * d$x + rnorm(n, 0, 0.6))
d$yh <- ifelse(runif(n) < plogis(-0.5 + 0.8 * d$l), 0, yy)
fb <- frm(bf(yh ~ x, hu ~ me(zo, sz)), family = hurdle_lognormal(), data = d)
nd <- data.frame(x = 0.5, zo = 1, sz = c(0.3, 3))
r <- fw(fb, nd)
show("hurdle_lognormal, hu draw", r)
sig <- as.numeric(frm_linpred(fb, newdata = nd, dpar = "sigma",
                              type = "response"))
eh <- dm(fb, nd, "hu"); bh <- slope(fb, nd[1, ], "hu", "zo")
ref <- vapply(1:2, function(i) {
  stats::integrate(function(z) (1 - plogis(eh[i] + bh * nd$sz[i] * z)) *
                     dnorm(z), -Inf, Inf, rel.tol = 1e-12)$value
}, 0) * exp(dm(fb, nd, "mu") + sig^2 / 2)
cat("  reference", ref, " rel", r$f[, "Estimate"] / ref - 1, "\n")

cat("\n== C. Gamma(log), shape ~ me()\n")
d <- mk(7303)
mu <- exp(0.5 + 0.3 * d$x)
k <- exp(1 + 0.5 * d$l)
d$yg <- rgamma(n, shape = k, rate = k / mu)
fc <- frm(bf(yg ~ x, shape ~ me(zo, sz)), family = Gamma(link = "log"),
          data = d)
r <- fw(fc, nd)
show("gamma, shape draw", r)
cat("  Est / mu - 1:", r$f[, "Estimate"] / exp(dm(fc, nd, "mu")) - 1, "\n")

cat("\n== D. the warning on two separate calls, and inside one\n")
d <- mk(7101)
set.seed(7101)
d <- data.frame(x = rnorm(n), l = rnorm(n), sz = 0.3)
d$zo <- d$l + rnorm(n, 0, 0.3)
d$y <- exp(0.5 + 0.4 * d$x + rnorm(n, 0, exp(-0.5 + 0.5 * d$l)))
f1 <- frm(bf(y ~ x, sigma ~ me(zo, sz)), family = lognormal(), data = d)
nd <- data.frame(x = 0, zo = 1, sz = 0.5)
for (i in 1:3) {
  r <- fw(f1, nd)
  cat("fitted() call", i, ": warnings", length(r$w), paste(r$w, collapse = ","),
      "\n")
}
r <- with_warn(lapply(1:2, function(i) fitted(f1, newdata = nd)))
cat("two fitted() calls inside one lapply():", length(r$warnings), "\n")
r <- with_warn(conditional_effects(f1, effects = "x",
                                   int_conditions = list(x = c(-1, 1)),
                                   conditions = data.frame(zo = 1, sz = 0.5)))
cat("conditional_effects() call:", length(r$warnings), "\n")
r <- with_warn(conditional_effects(f1, effects = "x",
                                   int_conditions = list(x = c(-1, 1)),
                                   conditions = data.frame(zo = 1, sz = 0.5)))
cat("second conditional_effects() call:", length(r$warnings), "\n")
r <- with_warn(predict(f1, newdata = nd))
cat("predict() after them:", length(r$warnings), "\n")
r <- fw(f1, nd)
cat("fitted() after them:", length(r$w), "\n")

cat("\n== E. y ~ sin(a)^2 + c at a = 0, draw sd pi\n")
set.seed(7208)
d <- data.frame(l1 = rnorm(n), l2 = rnorm(n), s1 = 0.2, s2 = 0.2)
d$x1 <- d$l1 + rnorm(n, 0, 0.2)
d$x2 <- d$l2 + rnorm(n, 0, 0.2)
d$y <- sin(0.3 + 0.8 * d$l1)^2 + (0.5 + 0.5 * d$l2) + rnorm(n, 0, 0.1)
fp <- suppressMessages(frm(bf(y ~ sin(a)^2 + c, a ~ me(x1, s1),
                              c ~ me(x2, s2), nl = TRUE),
                           family = gaussian(), data = d))
ba <- slope(fp, data.frame(x1 = 0, x2 = 0, s1 = 1, s2 = 1), "a", "x1")
a0 <- dm(fp, data.frame(x1 = 0, x2 = 0, s1 = 1, s2 = 1), "a")
for (sdv in c(pi, 2 * pi, pi / 1.7, 10 * pi, 1)) {
  nd <- data.frame(x1 = -a0 / ba, x2 = 0, s1 = sdv / abs(ba), s2 = 0.5)
  r <- fw(fp, nd)
  c0 <- dm(fp, nd, "c")
  ex <- (1 - exp(-2 * sdv^2)) / 2 + c0
  cat(sprintf("draw sd %.4f: fitted %.6f exact %.6f rel %.2e warn %s\n",
              sdv, r$f[, "Estimate"], ex, r$f[, "Estimate"] / ex - 1,
              paste(r$w, collapse = ",")))
}

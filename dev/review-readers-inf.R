# Review of lane readers: the infinite-mean rule of me_gh_moments().
# Data seeds per case below. brms's estimator at fixed values is simulated
# (seeds 1..20, 4000 draws each): the mean of exp(mu + sigma^2 / 2) over
# draws of sigma's me() value, which is what posterior_epred() averages.
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

cat("== 1. lognormal, sigma ~ me(zo, sz): seed 7101\n")
set.seed(7101)
n <- 800
d <- data.frame(x = rnorm(n), l = rnorm(n), sz = 0.3)
d$zo <- d$l + rnorm(n, 0, 0.3)
d$y <- exp(0.5 + 0.4 * d$x + rnorm(n, 0, exp(-0.5 + 0.5 * d$l)))
f1 <- frm(bf(y ~ x, sigma ~ me(zo, sz)), family = lognormal(), data = d)
for (sz in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8)) {
  nd <- data.frame(x = 0, zo = 1, sz = sz)
  r <- fw(f1, nd)
  es <- dm(f1, nd, "sigma"); bs <- slope(f1, nd, "sigma", "zo")
  cat(sprintf("sz %.2f (link sd %.3f): Estimate %s Est.Error %s Q %s %s warn %s\n",
              sz, bs * sz, format(r$f[, "Estimate"], digits = 6),
              format(r$f[, "Est.Error"], digits = 6),
              format(r$f[, "Q2.5"], digits = 4),
              format(r$f[, "Q97.5"], digits = 4),
              paste(r$w, collapse = ",")))
}
nd <- data.frame(x = 0, zo = 1, sz = 0.05)
es <- dm(f1, nd, "sigma"); bs <- slope(f1, nd, "sigma", "zo")
em <- dm(f1, nd, "mu")
cat("brms-style mean of 4000 draws at sz 0.05, seeds 1..20:\n")
print(signif(vapply(1:20, function(s) {
  set.seed(s)
  mean(exp(em + exp(2 * (es + bs * 0.05 * rnorm(4000))) / 2))
}, 0), 6))
cat("same at sz 0.8:\n")
print(signif(vapply(1:20, function(s) {
  set.seed(s)
  mean(exp(em + exp(2 * (es + bs * 0.8 * rnorm(4000))) / 2))
}, 0), 4))

cat("\n== 1b. conditional_effects() on that fit\n")
r <- with_warn(conditional_effects(f1, effects = "x",
                                   int_conditions = list(x = c(-1, 1)),
                                   conditions = data.frame(zo = 1, sz = 0.8)))
print(r$value[[1]][, c("x", "estimate__", "se__", "lower__", "upper__")])
cat("warnings:", paste(wclass(r$warnings), collapse = ","), "\n")

cat("\n== 2. lognormal, sigma ~ (1 | g), new level: seed 7102\n")
set.seed(7102)
J <- 30
d <- data.frame(x = rnorm(n), g = factor(sample(J, n, TRUE)))
d$y <- exp(0.5 + 0.4 * d$x + rnorm(n, 0, exp(-0.5 + rnorm(J, 0, 0.3)[d$g])))
f2 <- frm(bf(y ~ x, sigma ~ 1 + (1 | g)), family = lognormal(), data = d)
cat("sd(g) on sigma:", VarCorr(f2)$g$sd[1, 1], "\n")
r <- fw(f2, data.frame(x = 0, g = "new"), allow_new_levels = TRUE)
print(r$f)
cat("warnings:", paste(r$w, collapse = ","), "\n")

cat("\n== 3. poisson, y ~ me(x, sx), a log-linear mean, large sd: seed 7103\n")
set.seed(7103)
d <- data.frame(l = rnorm(n), sx = 0.3)
d$xo <- d$l + rnorm(n, 0, 0.3)
d$y <- rpois(n, exp(0.2 + 1 * d$l))
f3 <- frm(bf(y ~ me(xo, sx)), family = poisson(), data = d)
for (s in c(3, 4, 5, 6, 6.4, 6.6, 7, 8)) {
  nd0 <- data.frame(xo = 0, sx = 1)
  b <- slope(f3, nd0, "mu", "xo")
  nd <- data.frame(xo = 0, sx = s / b)
  r <- fw(f3, nd)
  e <- dm(f3, nd, "mu")
  cat(sprintf("link sd %.1f: Estimate / exp(eta + s^2/2) %s Est.Error / exact sd %s warn %s\n",
              s, format(r$f[, "Estimate"] / exp(e + s^2 / 2), digits = 6),
              format(r$f[, "Est.Error"] / sqrt((exp(s^2) - 1) *
                                                 exp(2 * e + s^2)),
                     digits = 4),
              paste(r$w, collapse = ",")))
}

cat("\n== 4. gamma and negbinomial with a draw on shape: seed 7104\n")
set.seed(7104)
d <- data.frame(x = rnorm(n), l = rnorm(n), sz = 0.3)
d$zo <- d$l + rnorm(n, 0, 0.3)
mu <- exp(0.5 + 0.3 * d$x)
k <- exp(1 + 0.5 * d$l)
d$yg <- rgamma(n, shape = k, rate = k / mu)
d$yn <- rnbinom(n, mu = mu, size = k)
nd <- data.frame(x = 0.5, zo = 1, sz = c(0.3, 3))
fg <- frm(bf(yg ~ x, shape ~ me(zo, sz)), family = Gamma(link = "log"),
          data = d)
r <- fw(fg, nd)
cat("gamma Estimate / mu - 1:", format(r$f[, "Estimate"] / exp(dm(fg, nd, "mu")) - 1,
                                   digits = 3), " Est.Error",
    format(r$f[, "Est.Error"], digits = 4), " warn", r$w, "\n")
fn <- frm(bf(yn ~ x, shape ~ me(zo, sz)), family = negbinomial(), data = d)
r <- fw(fn, nd)
cat("negbinomial Estimate / mu - 1:",
    format(r$f[, "Estimate"] / exp(dm(fn, nd, "mu")) - 1, digits = 3),
    " warn", r$w, "\n")

cat("\n== 5. lognormal, y ~ me(zo, sz) on mu only, large sd: seed 7105\n")
set.seed(7105)
d <- data.frame(l = rnorm(n), sz = 0.3)
d$zo <- d$l + rnorm(n, 0, 0.3)
d$y <- exp(0.5 + 1 * d$l + rnorm(n, 0, 0.5))
f5 <- frm(bf(y ~ me(zo, sz)), family = lognormal(), data = d)
for (s in c(1, 3, 5)) {
  nd0 <- data.frame(zo = 0, sz = 1)
  b <- slope(f5, nd0, "mu", "zo")
  nd <- data.frame(zo = 0, sz = s / b)
  sig <- as.numeric(frm_linpred(f5, newdata = nd, dpar = "sigma",
                                type = "response"))
  r <- fw(f5, nd)
  cat(sprintf("link sd %.1f: Estimate / closed form %s warn %s\n", s,
              format(r$f[, "Estimate"] / exp(dm(f5, nd, "mu") + s^2 / 2 +
                                               sig^2 / 2), digits = 8),
              paste(r$w, collapse = ",")))
}

cat("\n== 6. shifted_lognormal, sigma ~ me(zo, sz): seed 7106\n")
set.seed(7106)
d <- data.frame(l = rnorm(n), sz = 0.3)
d$zo <- d$l + rnorm(n, 0, 0.3)
d$y <- 0.3 + exp(0.2 + rnorm(n, 0, exp(-0.7 + 0.4 * d$l)))
f6 <- tryCatch(frm(bf(y ~ 1, sigma ~ me(zo, sz)),
                   family = shifted_lognormal(), data = d),
               error = function(e) { cat("fit error:", conditionMessage(e),
                                         "\n"); NULL })
if (!is.null(f6)) {
  r <- fw(f6, data.frame(zo = 1, sz = c(0.05, 0.8)))
  print(r$f)
  cat("warn:", paste(r$w, collapse = ","), "\n")
}

cat("\n== 2b. brms-style mean at case 2's estimates, seeds 1..20, 4000 draws\n")
s2 <- VarCorr(f2)$g$sd[1, 1]
nd2 <- data.frame(x = 0, g = "new")
em2 <- as.numeric(frm_linpred(f2, newdata = nd2, re_formula = NA))
es2 <- as.numeric(frm_linpred(f2, newdata = nd2, re_formula = NA,
                              dpar = "sigma"))
print(signif(vapply(1:20, function(s) {
  set.seed(s)
  mean(exp(em2 + exp(2 * (es2 + s2 * rnorm(4000))) / 2))
}, 0), 6))
cat("and with 4e6 draws, seed 1:")
set.seed(1)
cat(mean(exp(em2 + exp(2 * (es2 + s2 * rnorm(4e6))) / 2)), "\n")
cat("0.79.1's value, the median of sigma:", exp(em2 + exp(2 * es2) / 2), "\n")

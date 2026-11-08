# Review of lane readers, short check: the divergent cases the family-name
# rule misses. A: mixture(lognormal, lognormal), sigma1 ~ me() (data seed
# 7201, as dev/review-readers-sc.R), at larger sz. B: a lognormal whose
# sigma is a nonlinear body reading an nlpar with me() (data seed 7202).
# Independent reference: Golub-Welsch Gauss-Hermite at 11 and 41 nodes of
# component 1's mean, and brms-style means of 4000 draws, seeds 1..20.
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
n <- 800
set.seed(7201)
d <- data.frame(x = rnorm(n), l = rnorm(n), sz = 0.3)
d$zo <- d$l + rnorm(n, 0, 0.3)
comp <- runif(n) < 0.4
d$y <- ifelse(comp, exp(1.5 + rnorm(n, 0, exp(-0.7 + 0.5 * d$l))),
              exp(-0.5 + rnorm(n, 0, 0.4)))
fm <- frm(bf(y ~ 1, sigma1 ~ me(zo, sz)),
          family = mixture(lognormal(), lognormal()), data = d)
cat("dpars:", names(fm$spec$responses[[1]]$dpars), " family:",
    fm$spec$responses[[1]]$family$family, "\n")
for (sz in c(3, 6, 10)) {
  nd <- data.frame(zo = 1, sz = sz)
  r <- fw(fm, nd)
  es <- dm(fm, nd, "sigma1"); bs <- slope(fm, nd, "sigma1", "zo")
  em <- dm(fm, nd, "mu1")
  th1 <- as.numeric(frm_linpred(fm, newdata = nd, dpar = "theta1",
                                type = "response"))
  m2 <- as.numeric(frm_linpred(fm, newdata = nd, dpar = "mu2"))
  s2 <- as.numeric(frm_linpred(fm, newdata = nd, dpar = "sigma2",
                               type = "response"))
  rest <- (1 - th1) * exp(m2 + s2^2 / 2)
  bm <- range(vapply(1:20, function(s) {
    set.seed(s)
    th1 * mean(exp(em + exp(2 * (es + bs * sz * rnorm(4000))) / 2)) + rest
  }, 0))
  cat(sprintf(paste("sz %4.1f (sigma1 link sd %.3f): Estimate %.6g",
                    "Est.Error %.4g warn '%s' | brms-style 20 x 4000: %.4g to",
                    "%.4g\n"),
              sz, bs * sz, r$f[, "Estimate"], r$f[, "Est.Error"],
              paste(unique(r$w), collapse = ","), bm[1], bm[2]))
}

cat("\n== B. lognormal, sigma a body reading an nlpar with me()\n")
set.seed(7202)
d <- data.frame(x = rnorm(n), l = rnorm(n), sz = 0.3)
d$zo <- d$l + rnorm(n, 0, 0.3)
d$y <- exp(0.5 + 0.4 * d$x + rnorm(n, 0, exp(-0.5 + 0.5 * d$l)))
fs <- NULL
for (spec in list(
  quote(bf(y ~ x) + nlf(sigma ~ exp(s1)) + lf(s1 ~ me(zo, sz))),
  quote(bf(y ~ x, nlf(sigma ~ exp(s1)), lf(s1 ~ me(zo, sz)))))) {
  fs <- tryCatch(frm(eval(spec), family = lognormal(), data = d),
                 error = function(e) {
                   cat("  fit error for", deparse(spec), ":",
                       conditionMessage(e), "\n")
                   NULL
                 })
  if (!is.null(fs)) break
}
if (!is.null(fs)) {
  for (sz in c(0.3, 0.8)) {
    r <- fw(fs, data.frame(x = 0, zo = 1, sz = sz))
    cat(sprintf("sz %.1f: Estimate %.6g Est.Error %.4g warn '%s'\n", sz,
                r$f[, "Estimate"], r$f[, "Est.Error"],
                paste(unique(r$w), collapse = ",")))
  }
}

cat("\n== C. the class of frm_linpred(se.fit)'s refusal at a body's new level\n")
set.seed(8201)
n2 <- 300; J <- 12
d <- data.frame(x = runif(n2, 0, 3), g = factor(sample(J, n2, TRUE)))
d$y <- (3 + rnorm(J, 0, 0.4)[d$g]) * exp(-0.8 * d$x) + rnorm(n2, 0, 0.1)
fA <- frm(bf(y ~ a * exp(-b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE),
          family = gaussian(), data = d)
e <- tryCatch(frm_linpred(fA, newdata = data.frame(x = 1, g = "new"),
                          allow_new_levels = TRUE, se.fit = TRUE),
              error = function(e) e)
cat("classes:", class(e), "\n")

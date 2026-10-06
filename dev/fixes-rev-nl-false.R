# Reviewer of lane fixes, claim 2: identified nonlinear models the
# guard must let through, and unidentified ones it misses.
#   Rscript dev/fixes-rev-nl-false.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
run <- function(lab, expr) {
  r <- tryCatch({
    fit <- suppressMessages(suppressWarnings(expr))
    fe <- fixef(fit)[, "Estimate"]
    se <- fixef(fit)[, "Est.Error"]
    sprintf("FIT conv=%d  est %s  se %s", fit$opt$convergence,
            paste(sprintf("%s=%.4g", names(fe), fe), collapse = " "),
            paste(sprintf("%.3g", se), collapse = " "))
  }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 200)))
  cat(sprintf("%-52s %s\n", lab, r))
}
# 1. the test-nl-rtmb-scope.R response-preparation model, without its
# priors: an identified ML model (three distinct intercepts)
set.seed(7)
n <- 3000
pt <- stats::runif(n, 0, 0.6)
m1 <- 0.20; m2 <- 0.32; sdev <- 0.045
f1 <- stats::pnorm(pt, m1, sdev); f2 <- stats::pnorm(pt, m2, sdev)
p <- f2 * 0.95 + f1 * (1 - f2) * 0.2 + (1 - f1) * (1 - f2) * 0.5
d1 <- data.frame(pt = pt, y = stats::rbinom(n, 1, p))
form <- bf(y ~ pnorm(pt, m2, exp(ls)) * 0.95 +
             pnorm(pt, m1, exp(ls)) * (1 - pnorm(pt, m2, exp(ls))) * 0.2 +
             (1 - pnorm(pt, m1, exp(ls))) * (1 - pnorm(pt, m2, exp(ls))) * 0.5,
           m1 ~ 1, m2 ~ 1, ls ~ 1, nl = TRUE)
run("prep model, no prior, start given",
    frm(form, family = bernoulli(link = "identity"), data = d1,
        start = list(beta = c(0.15, 0.4, -3))))
run("prep model, the test's priors",
    frm(form, family = bernoulli(link = "identity"), data = d1,
        prior = prior(normal(0.15, 1), nlpar = "m1") +
          prior(normal(0.40, 1), nlpar = "m2") +
          prior(normal(-3, 1), nlpar = "ls")))
# 2. a psychometric function, pnorm(x, mu, sigma), mu ~ 1, ls ~ 1
set.seed(11)
x <- runif(800, -2, 2)
d2 <- data.frame(x = x, y = rbinom(800, 1, pnorm(x, 0.3, 0.7)))
run("psychometric pnorm(x, mu, exp(ls)), no prior",
    frm(bf(y ~ pnorm(x, m0, exp(ls)), m0 ~ 1, ls ~ 1, nl = TRUE),
        family = bernoulli(link = "identity"), data = d2,
        start = list(beta = c(0, 0))))
run("same with (x - mu) / exp(ls) in the first argument",
    frm(bf(y ~ pnorm((x - m0) / exp(ls)), m0 ~ 1, ls ~ 1, nl = TRUE),
        family = bernoulli(link = "identity"), data = d2,
        start = list(beta = c(0, 0))))
# 3. a Gaussian bump, dnorm(x, a, s) scaled
set.seed(12)
x <- runif(300, -3, 3)
d3 <- data.frame(x = x, y = 2 * dnorm(x, 0.5, 1.2) + rnorm(300, 0, 0.05))
run("bump h * dnorm(x, c0, exp(ls)), no prior",
    frm(bf(y ~ h * dnorm(x, c0, exp(ls)), h ~ 1, c0 ~ 1, ls ~ 1, nl = TRUE),
        data = d3, start = list(beta = c(1, 0, 0))))
# 4. the identified shapes the brief names
set.seed(13)
n <- 200
d4 <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, 20)))
d4$y <- 1 + 0.5 * d4$x + 0.4 * d4$z + rnorm(10)[d4$g] * 0.5 + rnorm(n, 0, 0.3)
d4$yp <- 2 * exp(0.3 * d4$x) + rnorm(n, 0, 0.2)
run("a + b, a ~ 1 + x, b ~ 0 + z", frm(bf(y ~ a + b, a ~ 1 + x, b ~ 0 + z,
                                          nl = TRUE), data = d4))
run("a * exp(b * x)", frm(bf(yp ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE),
                          data = d4, start = list(beta = c(1, 0))))
run("a + b, a ~ 1 + x, b ~ 0 + x + z (partial overlap)",
    frm(bf(y ~ a + b, a ~ 1 + x, b ~ 0 + x + z, nl = TRUE), data = d4))
run("a + b, a ~ 1 + x, b ~ 0 + z + (1 | g) (RE only in one)",
    frm(bf(y ~ a + b, a ~ 1 + x, b ~ 0 + z + (1 | g), nl = TRUE), data = d4))
run("a + b, a ~ 1 + (1 | g), b ~ 1 (shared intercept, RE in one)",
    frm(bf(y ~ a + b, a ~ 1 + (1 | g), b ~ 1, nl = TRUE), data = d4))
run("a + exp(b): different links, a ~ 1, b ~ 1",
    frm(bf(y ~ a + exp(b), a ~ 1 + x, b ~ 1, nl = TRUE), data = d4,
        start = list(beta = c(1, 0, -1))))
# 5. unidentified shapes the symbolic test can miss
run("MISS? a + 2 * b, shared intercepts",
    frm(bf(y ~ a + 2 * b, a ~ 1 + x, b ~ 1, nl = TRUE), data = d4))
run("MISS? a - b, shared intercepts",
    frm(bf(y ~ a - b, a ~ 1 + x, b ~ 1, nl = TRUE), data = d4))
run("b + a, shared intercepts",
    frm(bf(y ~ b + a, a ~ 1 + x, b ~ 1, nl = TRUE), data = d4))
run("(a + b) * x",
    frm(bf(y ~ (a + b) * x, a ~ 1, b ~ 1, nl = TRUE), data = d4))
run("MISS? inv_logit(a + b) * 3",
    frm(bf(y ~ 3 * inv_logit(a + b), a ~ 1, b ~ 1, nl = TRUE), data = d4))
run("a + b + 5 (constant)",
    frm(bf(y ~ a + b + 5, a ~ 1, b ~ 1, nl = TRUE), data = d4))
run("log(exp(a) * exp(b))",
    frm(bf(y ~ log(exp(a) * exp(b)), a ~ 1, b ~ 1, nl = TRUE), data = d4))
run("MISS? a * b (scale ridge, not a sum)",
    frm(bf(y ~ a * b, a ~ 1, b ~ 1, nl = TRUE), data = d4,
        start = list(beta = c(1, 1))))
# 6. the weak-prior exemption: a normal(0, 1000) on one parameter
f6 <- tryCatch(suppressMessages(suppressWarnings(
  frm(bf(y ~ a + b, a ~ 1 + x, b ~ 1 + x, nl = TRUE), data = d4,
      prior = set_prior("normal(0, 1000)", nlpar = "b")))),
  error = function(e) e)
if (inherits(f6, "error")) {
  cat("weak prior: ERROR", conditionMessage(f6), "\n")
} else {
  H <- optimHess(f6$opt$par, f6$obj$fn, f6$obj$gr)
  ev <- eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values
  cat("weak prior normal(0, 1000) on b: fits, conv", f6$opt$convergence,
      " Hessian eigenvalues", format(ev, digits = 3), " condition",
      format(max(ev) / min(ev), digits = 3), "\n")
  print(fixef(f6)[, 1:2])
}

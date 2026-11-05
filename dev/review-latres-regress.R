# Attack 5: gaussian and student paths, plain models. Run on both builds, output
# compared digit for digit (pin both to one core, dev/rtmb-pitfalls.md item 21).
# Usage: REV_ARM=lane|base Rscript review-latres-regress.R <outfile>
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
outf <- commandArgs(TRUE)[[1]]
set.seed(20261102)
d <- expand.grid(week = 1:6, subj = factor(1:30))
d <- d[order(d$subj, d$week), ]
n <- nrow(d)
d$x <- rnorm(n)
d$z <- rnorm(n)
d$f <- factor(sample(c("a", "b", "c"), n, TRUE))
e <- as.vector(replicate(30, arima.sim(list(ar = 0.5), 6, sd = 0.7)))
d$y <- 1 + 0.5 * d$x + e + rep(rnorm(30, 0, 0.3), each = 6)
d$t3 <- 1 + 0.5 * d$x + e + rt(n, 4) * 0.5
d$cnt <- rpois(n, exp(0.7 + 0.3 * d$x + rep(rnorm(30, 0, 0.3), each = 6)))
d$bin <- rbinom(n, 1, plogis(0.3 * d$x + rep(rnorm(30, 0, 0.3), each = 6)))
d$pos <- rgamma(n, 3, 3 / exp(0.3 * d$x))
specs <- list(
  ar1T = quote(frm(y ~ x + ar(week, subj, cov = TRUE), data = d)),
  ar2T = quote(frm(y ~ x + ar(week, subj, p = 2, cov = TRUE), data = d)),
  ma1T = quote(frm(y ~ x + ma(week, subj, cov = TRUE), data = d)),
  armaT = quote(frm(y ~ x + arma(week, subj, cov = TRUE), data = d)),
  arRE = quote(frm(y ~ x + (1 | subj) + ar(week, subj, cov = TRUE), data = d)),
  cosy = quote(frm(y ~ x + cosy(week, subj), data = d)),
  unstr = quote(frm(y ~ x + unstr(week, subj), data = d)),
  arF = quote(frm(y ~ x + ar(week, subj), data = d)),
  maF = quote(frm(y ~ x + ma(week, subj), data = d)),
  armaF = quote(frm(y ~ x + arma(week, subj, p = 2, q = 1), data = d)),
  stuT = quote(frm(t3 ~ x + ar(week, subj, cov = TRUE), family = student(), data = d)),
  stuF = quote(frm(t3 ~ x + ar(week, subj), family = student(), data = d)),
  sigAR = quote(frm(bf(y ~ x + ar(week, subj, cov = TRUE), sigma ~ x), data = d)),
  fcor = quote(frm(y ~ x + fcor(M), data = d, data2 = list(M = diag(n)))),
  lm = quote(frm(y ~ x * f + z, data = d)),
  lmer = quote(frm(y ~ x + (x | subj), data = d)),
  glmer = quote(frm(cnt ~ x + (1 | subj), family = poisson(), data = d)),
  bern = quote(frm(bin ~ x + (1 | subj), family = bernoulli(), data = d)),
  gam = quote(frm(pos ~ x + (1 | subj), family = Gamma(link = "log"), data = d)),
  smooth = quote(frm(y ~ s(x) + z, data = d)),
  ar1re = quote(frm(cnt ~ x + ar1(factor(week) + 0 | subj), family = poisson(), data = d)),
  mv = quote(frm(bf(y ~ x + (1 | subj)) + bf(cnt ~ x, family = poisson()), data = d)),
  quad = quote(frm(cnt ~ x + (1 | subj), family = poisson(), data = d, quadrature = TRUE)),
  imp = quote(frm(cnt ~ x + (1 | subj), family = poisson(), data = d, importance = 200, seed = 1)),
  gmm = quote(frm(y ~ x + (1 | subj), data = d, REML = TRUE)))
fmt <- function(v) paste(sprintf("%.15g", as.numeric(v)), collapse = " ")
con <- file(outf, "w")
for (nm in names(specs)) {
  m <- msg_of(eval(specs[[nm]]))
  writeLines(paste0("### ", nm), con)
  if (!is.null(m$err)) {
    writeLines(paste("ERR", substr(m$err, 1, 300)), con)
    next
  }
  f <- m$val
  writeLines(paste("warn", paste(m$warn, collapse = " | ")), con)
  writeLines(paste("par", fmt(f$opt$par)), con)
  writeLines(paste("names", paste(names(f$opt$par), collapse = " ")), con)
  writeLines(paste("ll", fmt(logLik(f))), con)
  writeLines(paste("code", f$opt$convergence, f$opt$iterations, f$opt$evaluations[1], f$opt$evaluations[2]), con)
  fx <- try(fitted(f), silent = TRUE)
  if (!inherits(fx, "try-error")) writeLines(paste("fit", fmt(as.matrix(fx)[1:5, 1])), con)
  writeLines(paste("vcov", fmt(diag(as.matrix(vcov(f, full = TRUE))))), con)
  r <- msg_of(residuals(f, type = "pearson"))
  if (is.null(r$err)) writeLines(paste("res", fmt(as.matrix(r$val)[1:5, 1])), con) else writeLines(paste("resERR", r$err), con)
  s <- msg_of(capture.output(print(summary(f))))
  writeLines(s$val %||% s$err, con)
  writeLines(paste("varnames", paste(variables(f), collapse = " ")), con)
  pr <- msg_of(default_prior(f))
  if (is.null(pr$err)) writeLines(capture.output(print(as.data.frame(pr$val))), con)
  sim <- msg_of(simulate(f, nsim = 1, seed = 5))
  if (is.null(sim$err)) writeLines(paste("sim", fmt(as.matrix(sim$val)[1:5, 1])), con) else writeLines(paste("simERR", sim$err), con)
  cv <- msg_of(confint_varcorr(f))
  if (is.null(cv$err)) writeLines(capture.output(print(cv$val)), con)
  a <- msg_of(autocor_matrix(f))
  if (is.null(a$err) && !is.null(a$val)) writeLines(paste("acm", fmt(a$val[1:2, 1:2])), con)
  writeLines(paste("tape", fmt(c(length(f$obj$env$par), f$obj$env$random %||% 0)[1])), con)
}
close(con)
# The error texts in the areas the lane rewrote
w <- function(e) tryCatch({e; "ok"}, error = function(z) conditionMessage(z))
cat(w(frm(y ~ x + ar(week, subj, cov = TRUE), family = student(), data = d)), "\n")

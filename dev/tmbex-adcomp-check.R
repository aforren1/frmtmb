# Audit harness for the adcomp-only TMB examples, at their published size.
#
# dev/tmb-examples-check.R covers the examples RTMB has ported: it runs
# each upstream R script as its own reference. The adcomp originals are
# C++ and cannot be run without a compiler, so the reference likelihood
# here is hand-written against RTMB from the .cpp template, then fitted
# and compared with frm() in the same process.
#
# The rows are the full-size ones. tests/testthat/test-tmb-examples.R
# carries the same models at test speed and needs no vendored data.
#
#   Rscript dev/tmbex-adcomp-check.R
#
# Needs dev/tmb-examples-adcomp/ (gitignored; see
# dev/tmb-examples-audit.md for how to refetch it).

self <- normalizePath(sub("^--file=", "",
                          grep("^--file=", commandArgs(), value = TRUE)[1]),
                      winslash = "/")
root <- normalizePath(file.path(dirname(self), ".."), winslash = "/")
AD <- file.path(root, "dev", "tmb-examples-adcomp")
suppressMessages(pkgload::load_all(root, quiet = TRUE))
library(RTMB)

say <- function(...) cat(sprintf(...), "\n", sep = "")
report <- function(nm, ref_ll, fit, t_ref, t_frm) {
  ll <- as.numeric(logLik(fit))
  say("%-12s ref=%.10f  frm=%.10f  |diff|=%.3e  ref %.1fs  frm %.1fs",
      nm, ref_ll, ll, abs(ref_ll - ll), t_ref, t_frm)
}

## ------------------------------------------------------------- orange_big
## Scaled-up Orange Tree: 5000 latent random effects, 35000 rows. The
## 192/726/356 constants in the template are starting-value offsets, so
## the fitted curve is the same one frm() reaches from `start`.
source(file.path(AD, "orange_data.R"))
d0 <- data_orange
mult <- d0$multiply
dd <- data.frame(y = rep(d0$y, mult), t = rep(d0$t, mult),
                 tree = factor(rep(seq_len(d0$M * mult),
                                   rep(d0$ngroup, mult))))
nll_orange <- function(p) {
  a0 <- 192 + p$beta[1] + p$u[as.integer(dd$tree)]
  f <- a0 / (1 + exp(-(dd$t - (726 + p$beta[2])) / (356 + p$beta[3])))
  -sum(dnorm(p$u, 0, exp(p$log_sigma_u), log = TRUE)) -
    sum(dnorm(dd$y, f, exp(p$log_sigma), log = TRUE))
}
t_ref <- system.time({
  o <- MakeADFun(nll_orange,
                 list(beta = c(0, 0, 0), log_sigma = 1, log_sigma_u = 2,
                      u = rep(0, d0$M * mult)),
                 random = "u", silent = TRUE)
  op <- nlminb(o$par, o$fn, o$gr, lower = c(-10, -10, -10, -5, -5),
               upper = c(10, 10, 10, 5, 5))})[3]
t_frm <- system.time(
  fit <- suppressWarnings(
    frm(bf(y ~ a0 / (1 + exp(-(t - a1) / a2)),
           a0 ~ 1 + (1 | tree), a1 ~ 1, a2 ~ 1, nl = TRUE),
        family = gaussian(), data = dd,
        start = list(beta = c(192, 726, 356)))))[3]
report("orange_big", -op$objective, fit, t_ref, t_frm)
say("             rows=%d latent=%d", nrow(dd), nlevels(dd$tree))

## ----------------------------------------------------------------- socatt
## Cumulative logit, 7 categories, 264 groups. The ADMB template writes
## the group effect as sigma * u with u ~ N(0, 1); (1 | g) is the same
## block in its own parameterization, so no pin is needed.
source(file.path(AD, "tools", "readdat.R"))
d <- readadmb(file.path(AD, "socatt.dat"))
X <- matrix(d$X, ncol = d$p, byrow = TRUE)
grp <- factor(rep(seq_len(d$M), each = 4))
yv <- d$y
S <- d$S
nll_soc <- function(p) {
  "[<-" <- ADoverload("[<-")
  alpha <- p$tmpk
  for (s in 2:length(p$tmpk)) alpha[s] <- alpha[s - 1] + exp(p$tmpk[s])
  eta <- as.vector(X %*% p$b) + exp(p$logsigma) * p$u[as.integer(grp)]
  P <- rep(0, length(yv))
  hi <- yv < S
  P[hi] <- 1 / (1 + exp(-(alpha[yv[hi]] - eta[hi])))
  P[!hi] <- 1
  lo <- yv > 1
  P[lo] <- P[lo] - 1 / (1 + exp(-(alpha[yv[lo] - 1] - eta[lo])))
  -sum(dnorm(p$u, 0, 1, log = TRUE)) - sum(log(1e-20 + P))
}
t_ref <- system.time({
  o <- MakeADFun(nll_soc,
                 list(b = rep(0, d$p), logsigma = 1, tmpk = rep(0, S - 1),
                      u = rep(0, nlevels(grp))),
                 random = "u", silent = TRUE)
  op <- nlminb(o$par, o$fn, o$gr)})[3]
ds <- data.frame(y = ordered(yv), g = grp)
for (j in seq_len(d$p)) ds[[paste0("x", j)]] <- X[, j]
fml <- stats::as.formula(paste("y ~",
                               paste(paste0("x", seq_len(d$p)),
                                     collapse = " + "), "+ (1 | g)"))
t_frm <- system.time(
  fit <- frm(bf(fml), family = cumulative(), data = ds))[3]
report("socatt", -op$objective, fit, t_ref, t_frm)

## ---------------------------------------------------------------- lr_test
## The example's point is TMB's map=. Here the restriction is a formula,
## and the parameter counts match, so anova() reproduces its table.
ngroup <- 5
nrep <- c(5, 8, 11, 13, 2)
set.seed(123)
ra <- lapply(seq_len(ngroup), function(i) rnorm(nrep[i], 0, c(1, 1, 1, 2, 2)[i]))
dl <- data.frame(obs = unlist(ra),
                 g = factor(rep(seq_len(ngroup), lengths(ra))))
nll_lr <- function(p) {
  gi <- as.integer(dl$g)
  -sum(dnorm(dl$obs, p$mu[gi], p$sd[gi], log = TRUE))
}
maps <- list(NULL, list(mu = factor(rep(1, 5))),
             list(mu = factor(rep(1, 5)), sd = factor(rep(1, 5))))
forms <- list(bf(obs ~ 0 + g, sigma ~ 0 + g), bf(obs ~ 1, sigma ~ 0 + g),
              bf(obs ~ 1, sigma ~ 1))
for (i in seq_along(maps)) {
  o <- MakeADFun(nll_lr, list(mu = rep(0, 5), sd = rep(1, 5)),
                 map = maps[[i]], silent = TRUE)
  op <- suppressWarnings(nlminb(o$par, o$fn, o$gr))
  f <- frm(forms[[i]], family = gaussian(), data = dl)
  say("lr_test[%d]  ref=%.10f  frm=%.10f  |diff|=%.3e  npar %d/%d",
      i, -op$objective, as.numeric(logLik(f)),
      abs(-op$objective - as.numeric(logLik(f))),
      length(o$par), length(f$opt$par))
}

## -------------------------------------------------------------- longlinreg
## 10^6 rows, three parameters. The example is about scale.
set.seed(123)
nobs <- 1e6
xl <- seq(0, 10, length = nobs)
dg <- data.frame(Y = 2 * xl + 1 + rnorm(nobs), x = xl)
nll_ll <- function(p) {
  -sum(dnorm(dg$Y, p$a + p$b * dg$x, exp(p$logSigma), log = TRUE))
}
t_ref <- system.time({
  o <- MakeADFun(nll_ll, list(a = 0, b = 0, logSigma = 0), silent = TRUE)
  op <- nlminb(o$par, o$fn, o$gr, o$he)})[3]
t_frm <- system.time(
  fit <- suppressWarnings(frm(bf(Y ~ x), family = gaussian(), data = dg)))[3]
report("longlinreg", -op$objective, fit, t_ref, t_frm)

## ------------------------------------- thetalog and mvrw: reference numbers
## Both are OUT OF SCOPE for frm() (no latent state-process term). Their
## reference likelihoods are printed so the audit's out-of-scope rows can
## be checked against something.
Y <- scan(file.path(AD, "thetalog.dat"), skip = 3, quiet = TRUE)
nll_tl <- function(p) {
  Xp <- head(p$X, -1)
  m <- Xp + exp(p$logr0) * (1 - (exp(Xp) / exp(p$logK))^exp(p$logtheta))
  -sum(dnorm(tail(p$X, -1), m, sqrt(exp(p$logQ)), log = TRUE)) -
    sum(dnorm(Y, p$X, sqrt(exp(p$logR)), log = TRUE))
}
o <- MakeADFun(nll_tl, list(X = Y * 0, logr0 = 0, logtheta = 0, logK = 6,
                            logQ = 0, logR = 0), random = "X", silent = TRUE)
TMB::newtonOption(o, smartsearch = FALSE)
op <- nlminb(o$par, o$fn, o$gr)
say("thetalog     reference nll=%.10f (out of scope: nonlinear state process)",
    op$objective)

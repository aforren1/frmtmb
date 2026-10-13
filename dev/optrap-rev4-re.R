# Reviewer of lane optrap, consolidation check: newton_polish() on fits
# with random effects, where its Hessian is optimHess() of the Laplace
# objective. Every design is fitted twice in one process: with the
# steps forced on (predictor_design_kappa() patched to Inf, so `ill`
# holds for every fit at code 0 with at most 100 parameters) and off
# (patched to 1). Per fit: objective change, steps, warnings, codes and
# the seconds newton_polish() took. Then the designs that trigger the
# steps by themselves (an uncentered year with year^2, a covariate near
# 1e5), as built, against lmer.
#   Rscript dev/optrap-rev4-re.R <lib> <label>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
lab <- a[2]
ns <- asNamespace("frmtmb")
cat(lab, "lib", dirname(find.package("frmtmb")), "\n")
kap0 <- get("predictor_design_kappa", ns)
np0 <- get("newton_polish", ns)
pol_secs <- 0
np_timed <- function(obj, opt, bounds, par_units = NULL, control = NULL,
                     ill = FALSE) {
  t0 <- proc.time()[[3]]
  r <- np0(obj, opt, bounds, par_units, control, ill = ill)
  pol_secs <<- proc.time()[[3]] - t0
  r
}
environment(np_timed) <- environment()
utils::assignInNamespace("newton_polish", np_timed, "frmtmb")
set_kappa <- function(v) {
  utils::assignInNamespace("predictor_design_kappa",
                           if (is.null(v)) kap0 else function(frame, ...) v,
                           "frmtmb")
}
fit1 <- function(fo, d, fam) {
  w <- character()
  pol_secs <<- 0
  f <- withCallingHandlers(suppressMessages(frm(fo, family = fam, data = d)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  list(obj = f$opt$objective, code = f$opt$convergence,
       steps = f$opt$polish[["steps"]] %||% 0,
       nwarn = sum(grepl("did not report|maximum absolute gradient|not available",
                         w)),
       secs = pol_secs, par = f$opt$par)
}
mk <- function(type, s) {
  set.seed(s)
  n <- 200
  d <- data.frame(g = factor(rep(1:20, each = 10)), x = rnorm(n),
                  h = factor(rep(1:10, 20)))
  b <- rnorm(20, 0, 0.6)[d$g]
  bs <- rnorm(20, 0, 0.3)[d$g]
  bh <- rnorm(10, 0, 0.4)[d$h]
  eta <- 0.2 + 0.5 * d$x + b
  switch(type,
    gauss_ri = {d$y <- eta + rnorm(n); list(y ~ x + (1 | g), d, gaussian())},
    gauss_rs = {d$y <- eta + bs * d$x + rnorm(n)
                list(y ~ x + (1 + x | g), d, gaussian())},
    crossed = {d$y <- eta + bh + rnorm(n)
               list(y ~ x + (1 | g) + (1 | h), d, gaussian())},
    binom = {d$k <- rbinom(n, 10, plogis(eta)); d$nt <- 10
             list(k | trials(nt) ~ x + (1 | g), d, binomial())},
    bern_rs = {d$y <- rbinom(n, 1, plogis(eta + bs * d$x))
               list(y ~ x + (1 + x | g), d, bernoulli())},
    pois = {d$y <- rpois(n, exp(eta)); list(y ~ x + (1 | g), d, poisson())},
    negbin = {d$y <- rnbinom(n, mu = exp(eta), size = 3)
              list(y ~ x + (1 | g), d, negbinomial())},
    student = {d$y <- eta + rt(n, 4); list(y ~ x + (1 | g), d, student())})
}
X <- list()
for (type in c("gauss_ri", "gauss_rs", "crossed", "binom", "bern_rs", "pois",
               "negbin", "student")) for (s in 1:5) {
  m <- mk(type, s)
  set_kappa(1)
  off <- fit1(m[[1]], m[[2]], m[[3]])
  set_kappa(Inf)
  on <- fit1(m[[1]], m[[2]], m[[3]])
  set_kappa(NULL)
  X[[length(X) + 1L]] <- data.frame(
    type = type, seed = s, dobj = on$obj - off$obj,
    steps = on$steps, code_off = off$code, code_on = on$code,
    warn_off = off$nwarn, warn_on = on$nwarn, secs = on$secs,
    dpar = max(abs(on$par - off$par)))
}
X <- do.call(rbind, X)
write.table(X, sprintf("dev/optrap-rev4-log/re-forced-%s.tsv", lab),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("forced steps on 40 fits with random effects:\n")
cat(sprintf("  fits with a step %d | objective lower %d, higher %d (max +%.3g) | warnings off -> on %d -> %d | codes != 0 off -> on %d -> %d | max polish secs %.3f | max |dpar| %.3g\n",
            sum(X$steps > 0), sum(X$dobj < -1e-12), sum(X$dobj > 1e-12),
            max(X$dobj), sum(X$warn_off > 0), sum(X$warn_on > 0),
            sum(X$code_off != 0), sum(X$code_on != 0), max(X$secs),
            max(X$dpar)))
print(X[X$steps > 0 | X$warn_off != X$warn_on, ], digits = 3)
# designs that trigger the steps as built
Y <- list()
for (s in 1:10) {
  set.seed(s)
  n <- 400
  d <- data.frame(g = factor(rep(1:40, each = 10)),
                  year = sample(1990:2020, n, TRUE), big = 1e5 + runif(n))
  d$y <- 0.002 * (d$year - 2005)^2 + 0.1 * (d$year - 2005) +
    rnorm(40, 0, 0.5)[d$g] + rnorm(n)
  d$y2 <- 0.5 * (d$big - 1e5) + rnorm(40, 0, 0.5)[d$g] + rnorm(n)
  for (des in list(list("year_sq", y ~ year + I(year^2) + (1 | g)),
                   list("big", y2 ~ big + (1 | g)))) {
    set_kappa(NULL)
    on <- fit1(des[[2]], d, gaussian())
    set_kappa(1)
    off <- fit1(des[[2]], d, gaussian())
    set_kappa(NULL)
    l4 <- suppressMessages(suppressWarnings(lmer(des[[2]], data = d,
                                                 REML = FALSE)))
    Y[[length(Y) + 1L]] <- data.frame(
      design = des[[1]], seed = s, kappa = kap0(frm(des[[2]], data = d,
                                                   dry_run = "frame")),
      steps = on$steps, secs = on$secs,
      dll_on = -on$obj - as.numeric(logLik(l4)),
      dll_off = -off$obj - as.numeric(logLik(l4)),
      warn_on = on$nwarn, warn_off = off$nwarn)
  }
}
Y <- do.call(rbind, Y)
write.table(Y, sprintf("dev/optrap-rev4-log/re-trigger-%s.tsv", lab),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(Y, digits = 3)

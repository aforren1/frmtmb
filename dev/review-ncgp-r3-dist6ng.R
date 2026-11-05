# Short check: distinct = 6L in gp_marginal_opt() and gp_marginal_max(), on the non-gaussian rows of review-ncgp-r3-ng.R, against the reference those runs read (each row's
# best of 15 starts). The start is the only thing that moves. A fit above the reference raises the reference.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
ns <- asNamespace("frmtmb")
for (nm in c("gp_marginal_opt", "gp_marginal_max")) {
  f <- ns[[nm]]; formals(f)$distinct <- 6L; environment(f) <- ns
  unlockBinding(nm, ns); assign(nm, f, envir = ns); lockBinding(nm, ns)
}
ref <- c(readLines(file.path(REV_OUT, "r3-ng-lane-a.txt")), readLines(file.path(REV_OUT, "r3-ng-lane-b.txt")))
refb <- c(readLines(file.path(REV_OUT, "r3-ng-base-a.txt")), readLines(file.path(REV_OUT, "r3-ng-base-b.txt")))
curves <- list(mix = function(x) 2 * sin(3 * x) + 0.5 * x, wig = function(x) sin(7 * x) + 0.2 * x, two = function(x) sin(x) + 0.4 * sin(10 * x), smooth = function(x) sin(x), step = function(x) 2 * (x > 5) - 1)
mk <- function(case, x, f, n) {
  eta <- f(x)
  switch(case,
    pois_low = list(d = data.frame(x = x, y = rpois(n, exp(eta / 2))), fml = y ~ gp(x), fam = poisson()),
    pois_mid = list(d = data.frame(x = x, y = rpois(n, exp(1.5 + eta / 2))), fml = y ~ gp(x), fam = poisson()),
    bern = list(d = data.frame(x = x, y = rbinom(n, 1, plogis(eta - mean(eta)))), fml = y ~ gp(x), fam = bernoulli()),
    binom5 = list(d = data.frame(x = x, tr = 5, y = rbinom(n, 5, plogis(eta - mean(eta)))), fml = y | trials(tr) ~ gp(x), fam = binomial()),
    negbin = list(d = data.frame(x = x, y = rnbinom(n, size = 20, mu = exp(4 + eta))), fml = y ~ gp(x), fam = negbinomial()),
    gamma = list(d = data.frame(x = x, y = rgamma(n, 400, 400 / exp(eta))), fml = y ~ gp(x), fam = Gamma(link = "log")),
    lognorm = list(d = data.frame(x = x, y = exp(eta + rnorm(n, sd = 0.02))), fml = y ~ gp(x), fam = lognormal()),
    student = list(d = data.frame(x = x, y = eta + 0.02 * rt(n, 5)), fml = y ~ gp(x), fam = student()),
    gsigma = list(d = data.frame(x = x, z = rnorm(n), y = eta + 0.02 * rnorm(n)), fml = bf(y ~ gp(x), sigma ~ z), fam = gaussian()),
    gar = list(d = data.frame(x = x, t = seq_len(n), y = eta + as.numeric(arima.sim(list(ar = 0.4), n, sd = 0.03))), fml = y ~ gp(x) + ar(time = t), fam = gaussian()))
}
tb <- list()
for (cs in c("pois_low", "pois_mid", "bern", "binom5", "negbin", "gamma", "lognorm", "student", "gsigma", "gar")) for (cn in names(curves)) for (n in c(40L, 80L, 160L)) {
  set.seed(74000 + 1 * 31 + n); x <- sort(round(runif(n, 0, 10), 4)); m <- mk(cs, x, curves[[cn]], n)
  r <- do.call(rev_fit, list(m$fml, m$d, family = m$fam))
  if (inherits(r$fit, "error")) next
  key <- sprintf("^RN3 lane %s %s %d 1 ", cs, cn, n); lr <- grep(key, ref, value = TRUE)
  best <- as.numeric(sub(".* best ([-0-9.]+).*", "\\1", lr)); ll4 <- as.numeric(sub(".* ll ([-0-9.]+) best.*", "\\1", lr))
  lb <- grep(sub("lane", "base", key), refb, value = TRUE); gb <- as.numeric(sub(".* gap ([-0-9.e+]+) conv.*", "\\1", lb))
  ll <- as.numeric(logLik(r$fit)); top <- max(ll, best)
  tb[[cs]] <- rbind(tb[[cs]], c(ll - top, ll4 - top, gb))
  if (ll - top < -1e-3) cat(sprintf("RD6N %s %s %d gap %.3g (distinct 4 %.3g, base %.3g) conv %d w %d\n", cs, cn, n, ll - top, ll4 - top, gb, r$fit$opt$convergence, length(r$warn)))
  flush.console()
}
for (k in names(tb)) { a <- tb[[k]]; cat(sprintf("RD6NSUM %-9s rows %d | distinct 6 short %d worst %.3g | distinct 4 short %d worst %.3g | base short %d worst %.3g | lane-only at 6: %d\n", k, nrow(a), sum(a[, 1] < -1e-3), min(a[, 1]), sum(a[, 2] < -1e-3), min(a[, 2]), sum(a[, 3] < -1e-3), min(a[, 3]), sum(a[, 1] < -1e-3 & a[, 3] >= -1e-3))) }
cat("DONE\n")

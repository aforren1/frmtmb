# Lane gpstop: the brief's silent rows, default fit against theta starts.
# Args: arm (base | lane) rows (comma list of E120w2, L120w2, F120w1,
#   st80, st160, or "all"). Row data: dev/review-fitedge-designs.R
#   (set.seed(52000 + 100 seed + n + 10 design + curve)) and dev/arstep-ng.R
#   (set.seed(74000 + 31 seed + n)).
a <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/gpstop-lib.R")
gs_libs(a[1])
rows <- strsplit(a[2], ",")[[1]]
source(file.path(GS_WT, "dev/review-fitedge-designs.R"))
ng_mk <- function(case, curve, n, seed) {
  set.seed(74000 + seed * 31 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  f <- switch(curve, step = function(x) 2 * (x > 5) - 1,
              two = function(x) sin(x) + 0.4 * sin(10 * x))
  eta <- f(x)
  switch(case,
    student = list(d = data.frame(x = x, y = eta + 0.02 * rt(n, 5)),
                   fml = y ~ gp(x), fam = student()),
    gsigma = list(d = data.frame(x = x, z = rnorm(n),
                                 y = eta + 0.02 * rnorm(n)),
                  fml = bf(y ~ gp(x), sigma ~ z), fam = gaussian()),
    pois_mid = list(d = data.frame(x = x, y = rpois(n, exp(1.5 + eta / 2))),
                    fml = y ~ gp(x), fam = poisson()))
}
defs <- list(
  E120w2 = function() list(d = rv_mk("E", 2, 120L, 2), fml = rv_spell$E$nl,
                           fam = gaussian(), th = function(s, r) c(s, r)),
  E60w2 = function() list(d = rv_mk("E", 2, 60L, 2), fml = rv_spell$E$nl,
                          fam = gaussian(), th = function(s, r) c(s, r)),
  L120w2 = function() list(d = rv_mk("L", 2, 120L, 2), fml = rv_spell$L$nl,
                           fam = gaussian(), th = function(s, r) c(s, r)),
  F120w1 = function() list(d = rv_mk("F", 1, 120L, 2), fml = rv_spell$F$nl,
                           fam = gaussian(), th = function(s, r) c(s, r)),
  st40 = function() c(ng_mk("student", "step", 40L, 1), th = NULL),
  st80 = function() c(ng_mk("student", "step", 80L, 1), th = NULL),
  st160 = function() c(ng_mk("student", "step", 160L, 1), th = NULL),
  gs40 = function() c(ng_mk("gsigma", "step", 40L, 1), th = NULL),
  pm80 = function() c(ng_mk("pois_mid", "two", 80L, 1), th = NULL))
if (identical(rows, "all")) rows <- names(defs)
show <- function(fit) {
  e <- fit$estimates
  paste(names(fit$opt$par), signif(fit$opt$par, 5), sep = "=", collapse = " ")
}
for (id in rows) {
  m <- defs[[id]]()
  t0 <- proc.time()[[3]]
  q <- fe_fit(m$fml, data = m$d, family = m$fam)
  el <- proc.time()[[3]] - t0
  if (inherits(q$fit, "error")) {
    cat("PROBE", id, "ERROR", conditionMessage(q$fit), "\n")
    next
  }
  ll0 <- as.numeric(logLik(q$fit))
  cat(sprintf("PROBE %s default ll %.6f conv %d nwarn %d wall %.2f\n", id, ll0,
              q$fit$opt$convergence, length(q$warn), el))
  if (length(q$warn)) cat("PROBEW", id, paste(q$warn, collapse = " || "), "\n")
  cat("PROBEP", id, "default", show(q$fit), "\n")
  best <- -Inf
  bp <- ""
  for (rho in c(0.1, 0.3, 1, 3, 8)) for (sd0 in c(0.3, 1.5, 6)) {
    th <- q$fit$opt$par[names(q$fit$opt$par) == "theta"]
    th[] <- c(log(sd0), log(rho))
    r <- fe_fit(m$fml, data = m$d, family = m$fam,
                start = list(theta = unname(th)))
    if (inherits(r$fit, "error")) next
    v <- as.numeric(logLik(r$fit))
    if (is.finite(v) && v > best) {
      best <- v
      bp <- sprintf("rho0 %g sd0 %g conv %d nwarn %d | %s", rho, sd0,
                    r$fit$opt$convergence, length(r$warn), show(r$fit))
    }
  }
  cat(sprintf("PROBE %s best of 15 theta starts %.6f gap %.4f\n", id, best,
              ll0 - best))
  cat("PROBEP", id, "best", bp, "\n")
  flush.console()
}
cat("DONE\n")

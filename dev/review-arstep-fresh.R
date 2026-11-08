# Review of lane arstep: fresh designs, not the lane's sets.
# Seeds 815000 + 97 * row index (one seed per row, new). x on [0, 10].
# Part G: gaussian + ar(time = t): curves step, smooth, periodic;
#   n 40, 100, 200, 300; AR coefficient 0, 0.4, 0.8 (innovation sd 0.1);
#   layout "ord" (t in x order) everywhere, and "rnd" (t unrelated to x)
#   at n 100. 45 rows.
# Part N: the working route. Families poisson, bernoulli, binomial(10),
#   student, gaussian with sigma ~ z, on the same curves, a latent AR(phi)
#   field of sd 0.3 in t = x order added to eta, model y ~ gp(x); n 40 and
#   120 by phi 0, 0.4, 0.8, and n 240 at phi 0.4. 105 rows. Part P:
#   poisson with y ~ gp(x) + ar(time = t), n 40 and 100, 18 rows.
# With ref 1 each row also runs multi-start references (start = theta, and
# thetaac for G), which bypass every start heuristic on both builds.
# Usage: Rscript review-arstep-fresh.R arm part ref [rowfrom rowto]
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/review-arstep-lib.R")
a <- commandArgs(TRUE)
arm <- a[1]; part <- a[2]; ref <- as.integer(a[3])
rv_libs(arm)
curves <- list(step = function(x) 2 * (x > 5) - 1,
               smooth = function(x) sin(x),
               periodic = function(x) sin(2 * pi * x / 2.5))
rows <- list()
if (part == "G") {
  for (cn in names(curves)) for (n in c(40L, 100L, 200L, 300L))
    for (phi in c(0, 0.4, 0.8)) {
      rows[[length(rows) + 1L]] <- list(fam = "gar", cn = cn, n = n, phi = phi, lay = "ord")
      if (n == 100L) rows[[length(rows) + 1L]] <- list(fam = "gar", cn = cn, n = n, phi = phi, lay = "rnd")
    }
} else if (part == "N") {
  for (fam in c("pois", "bern", "binom10", "student", "gsigma"))
    for (cn in names(curves)) {
      for (n in c(40L, 120L)) for (phi in c(0, 0.4, 0.8))
        rows[[length(rows) + 1L]] <- list(fam = fam, cn = cn, n = n, phi = phi, lay = "ord")
      rows[[length(rows) + 1L]] <- list(fam = fam, cn = cn, n = 240L, phi = 0.4, lay = "ord")
    }
} else if (part == "P") {
  for (cn in names(curves)) for (n in c(40L, 100L)) for (phi in c(0, 0.4, 0.8))
    rows[[length(rows) + 1L]] <- list(fam = "pois_ar", cn = cn, n = n, phi = phi, lay = "ord")
}
from <- if (length(a) > 3) as.integer(a[4]) else 1L
to <- if (length(a) > 4) min(as.integer(a[5]), length(rows)) else length(rows)
base_seed <- c(G = 815000L, N = 825000L, P = 835000L)[[part]]
arsim <- function(n, phi, sd) {
  if (phi == 0) rnorm(n, 0, sd) else as.numeric(arima.sim(list(ar = phi), n, sd = sd))
}
mk <- function(rw, seed) {
  set.seed(seed)
  n <- rw$n
  x <- sort(round(runif(n, 0, 10), 4))
  t <- if (rw$lay == "rnd") sample.int(n) else seq_len(n)
  f <- curves[[rw$cn]](x)
  switch(rw$fam,
    gar = {
      e <- numeric(n); e[order(t)] <- arsim(n, rw$phi, 0.1)
      list(d = data.frame(x = x, t = t, y = f + e), fml = y ~ gp(x) + ar(time = t), fam = gaussian())
    },
    pois = list(d = data.frame(x = x, y = rpois(n, exp(1 + f / 2 + arsim(n, rw$phi, 0.3)))), fml = y ~ gp(x), fam = poisson()),
    bern = list(d = data.frame(x = x, y = rbinom(n, 1, plogis(f + arsim(n, rw$phi, 0.3)))), fml = y ~ gp(x), fam = bernoulli()),
    binom10 = list(d = data.frame(x = x, tr = 10, y = rbinom(n, 10, plogis(f + arsim(n, rw$phi, 0.3)))), fml = y | trials(tr) ~ gp(x), fam = binomial()),
    student = list(d = data.frame(x = x, y = f + arsim(n, rw$phi, 0.1) + 0.05 * rt(n, 4)), fml = y ~ gp(x), fam = student()),
    gsigma = {
      z <- rnorm(n)
      list(d = data.frame(x = x, z = z, y = f + arsim(n, rw$phi, 0.1) * exp(0.3 * z)), fml = bf(y ~ gp(x), sigma ~ z), fam = gaussian())
    },
    pois_ar = list(d = data.frame(x = x, t = t, y = rpois(n, exp(1 + f / 2 + arsim(n, rw$phi, 0.3)))), fml = y ~ gp(x) + ar(time = t), fam = poisson()))
}
fitit <- function(m, start = NULL) {
  args <- list(m$fml, data = m$d, family = m$fam)
  if (!is.null(start)) args$start <- start
  t0 <- proc.time()[[3]]
  r <- rv_quiet(do.call(frm, args))
  r$wall <- proc.time()[[3]] - t0
  r
}
for (i in seq(from, to)) {
  rw <- rows[[i]]
  seed <- base_seed + 97L * i
  m <- mk(rw, seed)
  lbl <- sprintf("%s %s %s %d %.1f %s", part, rw$fam, rw$cn, rw$n, rw$phi, rw$lay)
  r <- fitit(m)
  if (inherits(r$value, "error")) {
    cat(sprintf("FR %s %d %s ERROR %s\n", arm, i, lbl, substr(conditionMessage(r$value), 1, 120)))
    next
  }
  fit <- r$value
  ll <- as.numeric(logLik(fit))
  th <- fit$estimates$theta
  ac <- fit$estimates$thetaac
  cat(sprintf("FR %s %d %s ll %.6f conv %d w %d wall %.2f th %s ac %s second %s\n",
              arm, i, lbl, ll, fit$opt$convergence, length(r$warnings), r$wall,
              paste(signif(th, 5), collapse = ","),
              if (length(ac)) paste(signif(ac, 5), collapse = ",") else "-",
              if (is.null(fit$opt$gp_second)) "-" else
                paste(signif(fit$opt$gp_second, 4), collapse = ",")))
  if (length(r$warnings)) {
    cat(sprintf("FRW %s %d %s %s\n", arm, i, lbl,
                paste(unique(substr(r$warnings, 1, 80)), collapse = " || ")))
  }
  if (ref) {
    best <- ll; bst <- "default"
    gpb <- Filter(function(b) identical(b$covstruct, "gp"), fit$frame$re_blocks)
    sds <- if (rw$n >= 200L) c(0.5, 2) else c(0.3, 1.5, 6)
    rhos <- if (rw$n >= 200L) c(0.1, 0.6, 3) else c(0.1, 0.5, 3)
    acs <- if (length(ac)) c(0, 0.95) else NA
    for (sd0 in sds) for (rho0 in rhos) for (ac0 in acs) {
      th0 <- th
      for (b in gpb) { th0[b$theta_idx[1]] <- log(sd0); th0[b$theta_idx[-1]] <- log(rho0) }
      st <- list(theta = th0)
      if (!is.na(ac0)) st$thetaac <- rep(ac0, length(ac)) * c(1, numeric(length(ac) - 1))
      r2 <- fitit(m, st)
      if (inherits(r2$value, "error")) next
      v <- as.numeric(logLik(r2$value))
      if (is.finite(v) && v > best) {
        best <- v
        bst <- sprintf("sd %g rho %g ac %g", sd0, rho0, ac0)
      }
    }
    cat(sprintf("FRR %s %d %s best %.6f from %s\n", arm, i, lbl, best, bst))
  }
  flush.console()
}
cat("DONE\n")

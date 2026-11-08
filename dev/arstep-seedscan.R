# Which seed of dev/review-ncgp-r3-ng.R is the review's 25.9 row? Same data
# construction (seed 74000 + 31 seed + n), step curve, `gar` and `student`;
# default fit against a few short-length-scale starts. Args: arm seeds
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/arstep-lib.R")
arm <- args[1]; seeds <- eval(parse(text = args[2]))
rev_libs(arm)
f <- function(x) 2 * (x > 5) - 1
mk <- function(case, x, n) {
  eta <- f(x)
  switch(case,
    student = list(d = data.frame(x = x, y = eta + 0.02 * rt(n, 5)), fml = y ~ gp(x), fam = student()),
    gar = list(d = data.frame(x = x, t = seq_len(n), y = eta + as.numeric(arima.sim(list(ar = 0.4), n, sd = 0.03))), fml = y ~ gp(x) + ar(time = t), fam = gaussian()))
}
fitit <- function(m, start = NULL) { a <- list(m$fml, m$d, family = m$fam); if (!is.null(start)) a$start <- start; do.call(rev_fit, a) }
for (seed in seeds) for (cs in c("gar", "student")) for (n in c(40L, 80L, 160L)) {
  set.seed(74000 + seed * 31 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  m <- mk(cs, x, n)
  r <- fitit(m)
  if (inherits(r$fit, "error")) { cat(sprintf("ASS %s %s %d %d ERROR\n", arm, cs, n, seed)); next }
  ll0 <- as.numeric(logLik(r$fit)); th0 <- r$fit$estimates$theta; best <- ll0
  for (rho0 in c(0.1, 0.3, 1)) for (sd0 in c(1.5)) {
    th <- th0
    for (b in r$fit$frame$re_blocks) if (identical(b$covstruct, "gp")) { th[b$theta_idx[1]] <- log(sd0); th[b$theta_idx[-1]] <- log(rho0) }
    r2 <- fitit(m, list(theta = th))
    if (!inherits(r2$fit, "error")) { v <- as.numeric(logLik(r2$fit)); if (is.finite(v) && v > best) best <- v }
  }
  cat(sprintf("ASS %s %s %d %d gap %.3e conv %d w %d ll %.4f best %.4f th %s\n", arm, cs, n, seed, ll0 - best, r$fit$opt$convergence, length(r$warn), ll0, best, paste(round(th0, 3), collapse = ",")))
  flush.console()
}
cat("DONE\n")

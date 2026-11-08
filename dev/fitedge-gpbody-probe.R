# Item 1 probe: gp() in a nonlinear body against the linear spelling of the
# same model. Args: arm. Seeds 1811 (the release review's x18.R D seed) and
# 58001 (the ncgp review's edge2 row).
arm <- commandArgs(TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
mk <- function(seed, n, curve) {
  set.seed(seed)
  x <- sort(runif(n, 0, 10)); z <- rnorm(n)
  f <- switch(curve, sin = sin(x), wig = 2 * sin(3 * x) + x,
              exp = exp(0.3 * x), step = sign(sin(1.5 * x)))
  data.frame(x, z, y = f + 0.3 * z + 0.1 * rnorm(n))
}
ns <- asNamespace("frmtmb")
for (cs in list(list(1811, 60, "sin"), list(1811, 60, "wig"),
                list(58001, 60, "wig"))) {
  d <- mk(cs[[1]], cs[[2]], cs[[3]])
  t0 <- proc.time()[[3]]
  lin <- fe_fit(y ~ z + gp(x), data = d)
  t1 <- proc.time()[[3]]
  nl <- fe_fit(bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE), data = d)
  t2 <- proc.time()[[3]]
  cat(sprintf("P1 %s seed %d %s lin %.4f conv %d w %d (%.1fs) | nl %.4f conv %d w %d (%.1fs) | gap %.3f\n",
              arm, cs[[1]], cs[[3]], logLik(lin$fit), lin$fit$opt$convergence,
              length(lin$warn), t1 - t0, logLik(nl$fit), nl$fit$opt$convergence,
              length(nl$warn), t2 - t1, logLik(nl$fit) - logLik(lin$fit)))
  if (length(nl$warn)) cat("  W:", substr(nl$warn, 1, 150), sep = "\n  ")
  cat("  lin theta", signif(lin$fit$estimates$theta, 4), " nl theta",
      signif(nl$fit$estimates$theta, 4), "\n")
}
fr <- nl$fit$frame
for (k in names(fr$linpreds)) {
  l <- fr$linpreds[[k]]
  cat("LP", k, "dpar", l$dpar %||% "-", "par", l$par %||% "-", "nlpar",
      l$nlpar %||% "-", "body", !is.null(l$nl_body), "Zcol", NCOL(l$Z), "\n")
}
cat("blocks", vapply(fr$re_blocks, `[[`, "", "covstruct"), "\n")
cat("DONE\n")

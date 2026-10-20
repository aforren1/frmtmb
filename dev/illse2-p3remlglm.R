# Lane illse2, item 3: REML (Laplace over beta and the group effects)
# on non-gaussian raw polynomials with (1 | g): does the fit converge,
# lose nothing, warn nothing, and agree with the same model written with
# poly(x, deg), whose design is orthogonal (fitted SEs on the response
# scale at x = 1, 1.5, 2; the variance parameters, which a change of
# basis of the fixed effects leaves where they are)? No exact reference:
# the Laplace REML is not
# a closed form.
#   ILLSE2_ARM=lane|base Rscript dev/illse2-p3remlglm.R <family> <deg> <seeds>
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
a <- commandArgs(TRUE)
fam <- a[1]
deg <- as.integer(a[2])
seeds <- eval(parse(text = a[3]))
xn <- paste0("x", seq_len(deg))
fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "),
                              "+ (1 | g)"))
fop <- stats::as.formula(paste("y ~ poly(x,", deg, ") + (1 | g)"))
famo <- switch(fam, bernoulli = bernoulli(), poisson = poisson())
for (seed in seeds) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(400, 1, 2),
                  g = factor(rep(1:20, length.out = 400)))
  eta <- sin(3 * d$x) + stats::rnorm(20, 0, 0.5)[d$g]
  d$y <- switch(fam, bernoulli = stats::rbinom(400, 1, stats::plogis(eta)),
                poisson = stats::rpois(400, exp(eta)))
  for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
  r <- cap(frm(fo, data = d, family = famo, REML = TRUE))
  rp <- cap(frm(fop, data = d, family = famo, REML = TRUE))
  if (is_err(r$v) || is_err(rp$v)) {
    cat("seed", seed, "ERROR", if (is_err(r$v)) r$v else rp$v, "\n")
    next
  }
  nd <- data.frame(x = c(1, 1.5, 2), g = d$g[1])
  for (k in seq_len(deg)) nd[[xn[k]]] <- nd$x^k
  f1 <- cap(fitted(r$v, newdata = nd, re_formula = NA))
  f2 <- cap(fitted(rp$v, newdata = nd, re_formula = NA))
  cat(sprintf(paste0("%s deg %d seed %d | code %d (poly %d) | theta - poly ",
                     "%+.1e | lost %s | whitened %s | fitted SE / poly %s | ",
                     "warn %s\n"),
              fam, deg, seed, r$v$opt$convergence, rp$v$opt$convergence,
              max(abs(r$v$opt$par - rp$v$opt$par)),
              lost_of(r$v), !is.null(r$v$beta_whiten),
              if (is_err(f1$v) || is_err(f2$v)) "ERROR" else
                ratio_str(f1$v[, "Est.Error"] / f2$v[, "Est.Error"]),
              wtag(r$w)))
}

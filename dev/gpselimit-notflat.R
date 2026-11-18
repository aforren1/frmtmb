# Lane gpselimit: approximate gp() fits whose length scale sits under a
# tenth of the basis's spacing (frmtmb_gp_basis fires) at an interior
# maximum, so the ridge to the white limit is NOT flat. The basis
# coefficients are set to their prior sd at a chosen length scale ls0 =
# frac * spacing, with alternating signs, so that the variance profile's
# maximum sits at ls0; the noise is small so that the coefficients are
# determined. For each fit: the length scale over the spacing, the
# objective along the ridge (log sdgp up D/2, log lscale down 1 per unit t)
# at t = -2, 2 and 12 against the fit, the SEs and lost reasons, and the
# conditions.
# Usage: Rscript gpselimit-notflat.R <lane|base> [jitter]
source("C:/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/gpselimit-lib.R")
a <- commandArgs(TRUE)
arm <- a[1]
jit <- if (length(a) >= 2) as.numeric(a[2]) else 0
gs_libs(arm)
ns <- asNamespace("frmtmb")
ridge <- function(fit, bk) {
  om <- ns$outer_par_map(fit)
  thpos <- which(om$comp == "theta")[bk$theta_idx]
  p <- fit$opt$par
  r <- numeric(length(p))
  r[thpos[1]] <- bk$gp_D / 2
  r[thpos[2]] <- -1
  f0 <- fit$obj$fn(p)
  v <- vapply(c(-2, 2, 12), function(t) fit$obj$fn(p + t * r) - f0, 0)
  fit$obj$fn(p)
  v
}
n <- 200
set.seed(4100)
x <- sort(runif(n, 0, 10))
for (K in c(10, 30, 60)) for (frac in c(0.02, 0.05, 0.08)) {
  d0 <- data.frame(x = x, y = rnorm(n))
  f0 <- suppressWarnings(suppressMessages(
    frm(bf(y ~ gp(x, k = K)), data = d0, family = gaussian())))
  bk <- Filter(function(b) b$covstruct == "hsgp", f0$frame$re_blocks)[[1]]
  lp <- f0$frame$linpreds[[ns$gp_block_lp(bk, f0$frame)]]
  Z <- as.matrix(ns$lp_Z(f0, lp))[, bk$c_idx, drop = FALSE]
  sp <- bk$gp_spacing
  th0 <- c(0, log(frac * sp))
  sds <- as.vector(ns$hsgp_sds(th0, bk$aux_omega, bk$gp_iso,
                               ns$gp_cov_of(bk)))
  set.seed(4200 + K)
  b <- sds * rep(c(1, -1), length.out = length(sds))
  y <- as.vector(Z %*% b) + 0.01 * rnorm(n)
  if (jit) y <- y * (1 + 2e-15 * jit)
  d <- data.frame(x = x, y = y)
  q <- gs_fit(bf(y ~ gp(x, k = K)), data = d, family = gaussian())
  fit <- q$fit
  if (inherits(fit, "error")) {
    cat("NOTFLAT K", K, "frac", frac, "ERROR", conditionMessage(fit), "\n")
    next
  }
  bk <- Filter(function(b) b$covstruct == "hsgp", fit$frame$re_blocks)[[1]]
  th <- fit$estimates$theta[bk$theta_idx]
  rg <- ridge(fit, bk)
  ss <- gs_ses(fit)
  tn <- grepl("^theta", names(ss$se))
  cat(sprintf(paste0("NOTFLAT %s K %2d frac %.2f jit %g | ll %.6f conv %d | ",
                     "ls/spacing %.4g sd %.4g | ridge f(t)-f0 t=-2 %.3g ",
                     "t=2 %.3g t=12 %.3g | se(theta) %s | lost [%s] | ",
                     "warn [%s] | msg [%s]\n"),
              arm, K, frac, jit, as.numeric(logLik(fit)),
              fit$opt$convergence, exp(th[2]) / sp, exp(th[1]),
              rg[1], rg[2], rg[3],
              paste(format(ss$se[tn], digits = 4), collapse = ","),
              paste(paste0(names(ss$lost), "=", ss$lost), collapse = ","),
              paste(unique(q$wclass), collapse = ";"),
              paste(unique(q$mclass), collapse = ";")))
}

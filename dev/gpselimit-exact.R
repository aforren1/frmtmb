# Lane gpselimit: the exact gp() (no k) at its white limit. Data: an iid
# effect per position (sd 1) plus noise (sd 0.3), so the field the data
# hold is white at every length scale above none. `rep` rows per position:
# with 1 the field and the noise are confounded (only sdgp^2 + sigma^2 is
# determined); with 3 they are not. Per fit: the length scale over the
# smallest gap, the SEs and lost reasons, the boundary kinds, the
# conditions, and the objective with the length scale moved -2 and -12
# (rest held) against the fit.
# Usage: Rscript gpselimit-exact.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/gpselimit-lib.R")
a <- commandArgs(TRUE)
arm <- a[1]
gs_libs(arm)
ns <- asNamespace("frmtmb")
for (rep in c(1, 3)) for (seed in 1:4) {
  set.seed(7700 + 10 * rep + seed)
  np <- 40
  xs <- sort(runif(np, 0, 10))
  u <- rnorm(np)
  x <- rep(xs, each = rep)
  y <- rep(u, each = rep) + 0.3 * rnorm(np * rep)
  q <- gs_fit(bf(y ~ gp(x)), data = data.frame(x = x, y = y),
              family = gaussian())
  fit <- q$fit
  if (inherits(fit, "error")) {
    cat("EXACT rep", rep, "seed", seed, "ERROR", conditionMessage(fit), "\n")
    next
  }
  bk <- Filter(function(b) b$covstruct == "gp", fit$frame$re_blocks)[[1]]
  th <- fit$estimates$theta[bk$theta_idx]
  om <- ns$outer_par_map(fit)
  tl <- which(om$comp == "theta")[bk$theta_idx[2]]
  p <- fit$opt$par
  f0 <- fit$obj$fn(p)
  mv <- vapply(c(-2, -12, 2), function(t) {
    pp <- p
    pp[tl] <- pp[tl] + t
    fit$obj$fn(pp) - f0
  }, 0)
  fit$obj$fn(p)
  ss <- gs_ses(fit)
  gap <- min(diff(sort(unique(x))))
  cat(sprintf(paste0("EXACT %s rep %d seed %d | ll %.6f conv %d | ",
                     "ls/mingap %.3g sd %.4g sigma %.4g | f(ls t)-f0 ",
                     "t=-2 %.3g t=-12 %.3g t=2 %.3g | se %s | lost [%s] | ",
                     "kind [%s] | warn [%s] | msg [%s]\n"),
              arm, rep, seed, as.numeric(logLik(fit)), fit$opt$convergence,
              exp(th[2]) * (bk$gp_dmax %||% 1) / gap, exp(th[1]),
              exp(fit$estimates$betad[1]), mv[1], mv[2], mv[3],
              paste(format(ss$se, digits = 3), collapse = ","),
              paste(paste0(names(ss$lost), "=", ss$lost), collapse = ","),
              paste(paste0(names(fit$cache$se_boundary_kind), "=",
                           fit$cache$se_boundary_kind), collapse = ","),
              paste(unique(q$wclass), collapse = ";"),
              paste(unique(q$mclass), collapse = ";")))
}

# Review of lane gpselimit, check 7: the exact gp() MEDIUM. Design of
# dev/gpselimit-exact2.R (iid effect per position plus noise, 40
# positions, 3 rows each), seeds 1 to 4: the length scale over the
# smallest gap, the objective with log lscale moved by -2 and -12, the
# log-scale SEs and the lost reasons, in a fresh process.
# Usage: Rscript gpselimit-review-exact.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/reviews/gpselimit-review-lib.R")
arm <- commandArgs(TRUE)[1]
rv_libs(arm)
ns <- asNamespace("frmtmb")
for (seed in 1:4) {
  rep <- 3
  set.seed(7700 + 10 * rep + seed)
  np <- 40
  xs <- sort(runif(np, 0, 10))
  u <- rnorm(np)
  x <- rep(xs, each = rep)
  y <- rep(u, each = rep) + 0.3 * rnorm(np * rep)
  q <- rv_fit(frm(bf(y ~ gp(x)), data = data.frame(x = x, y = y),
                  family = gaussian()))
  fit <- q$fit
  bk <- Filter(function(b) identical(b$covstruct, "gp"),
               fit$frame$re_blocks)[[1]]
  th <- fit$estimates$theta[bk$theta_idx]
  gap <- min(diff(xs))
  om <- ns$outer_par_map(fit)
  pos <- which(om$comp == "theta")[bk$theta_idx]
  p <- fit$opt$par
  pf <- ns$se_probe_fn(fit)
  f0 <- pf(p)
  mv <- vapply(c(-2, -12, 2), function(t) {
    qq <- p; qq[pos[2]] <- qq[pos[2]] + t; pf(qq) - f0
  }, 0)
  pf(p)
  ss <- rv_ses(fit)
  cat(sprintf(paste0("EXACT %s seed %d | ll %.6f | sdgp %.4g lscale/gap %.3g",
                     " | f(log ls -2, -12, +2) - f0 %s | se %s | lost [%s]",
                     " | conds [%s]\n"),
              arm, seed, as.numeric(logLik(fit)), exp(th[1]),
              exp(th[2]) / gap, paste(format(mv, digits = 3), collapse = ","),
              paste(format(ss$se, digits = 4), collapse = ","),
              paste(paste0(names(ss$lost), "=", ss$lost), collapse = ","),
              rv_cstr(q$conds)))
}

# rev-ndt: is "sd(ndt) is recovered" met by the MODEL or by the FLOORS,
# and does the answer move with the number of trials?
#
# Three arms on the SAME data, at the tier's own design and seed:
#
#   re      ndt ~ 1 + (1 | s) with ndt_group(s)   what the lane ships
#   nore    ndt ~ 1           with ndt_group(s)   the floors alone
#   global  ndt ~ 1 + (1 | s) with no grouping    the 0.6.0 model
#
# `nore` is the estimator the acceptance criterion has to beat: it has
# NO random effect on ndt at all, so every scrap of between-subject
# spread it reports comes from the per-subject floors. If it matches
# `re`, the criterion "recovers sd(ndt)" is met by the data.
#
# Trial counts 100, 200 and 400, because the fastest of n draws is a
# biased ceiling whose bias shrinks with n, and an estimator whose
# answer moves with n for that reason is worth knowing about.
#
# Seed 20260908, the tier's own. Worktree build.

.libPaths(c("C:/Users/adf44/source/r/rev-ndt-lib",
            "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})

tr <- list(mu0 = 0.4, mu_cond = 0.9, bs = 1.4, ndt = 0.25,
           sd_mu = 0.35, sd_lbs = 0.20, sd_lndt = 0.12)

mk <- function(nt, ns = 30L, seed = 20260908L) {
  set.seed(seed)
  u_mu <- rnorm(ns, 0, tr$sd_mu)
  u_bs <- rnorm(ns, 0, tr$sd_lbs)
  u_nd <- rnorm(ns, 0, tr$sd_lndt)
  s <- rep(seq_len(ns), each = nt)
  cond <- rep(rep(0:1, each = nt / 2L), times = ns)
  d <- ddm_simulate(ns * nt,
                    mu = tr$mu0 + tr$mu_cond * cond + u_mu[s],
                    bs = tr$bs * exp(u_bs[s]),
                    ndt = tr$ndt * exp(u_nd[s]), bias = 0.5, sv = 0)
  d$s <- factor(s)
  d$cond <- factor(cond, labels = c("a", "b"))
  attr(d, "truth") <- tr$ndt * exp(u_nd)
  d
}

el <- function(e) { t0 <- Sys.time()
  force(e); as.numeric(difftime(Sys.time(), t0, units = "secs")) }

row <- function(nt, arm) {
  d <- mk(nt)
  truth <- attr(d, "truth")
  one <- d[match(levels(d$s), as.character(d$s)), , drop = FALSE]
  own <- as.numeric(tapply(d$rt, d$s, min))
  grp <- arm != "global"
  re <- arm != "nore"
  lhs <- if (grp) quote(rt | dec(upper) + ndt_group(s)) else
    quote(rt | dec(upper))
  form <- if (grp && re) {
    bf(rt | dec(upper) + ndt_group(s) ~ cond + (1 | s), bs ~ 1 + (1 | s),
       ndt ~ 1 + (1 | s), bias = 0.5)
  } else if (grp && !re) {
    bf(rt | dec(upper) + ndt_group(s) ~ cond + (1 | s), bs ~ 1 + (1 | s),
       ndt ~ 1, bias = 0.5)
  } else {
    bf(rt | dec(upper) ~ cond + (1 | s), bs ~ 1 + (1 | s),
       ndt ~ 1 + (1 | s), bias = 0.5)
  }
  fit <- NULL
  secs <- el(fit <- frm(form, family = wiener(), data = d))
  # ndt_time() and not predict(dpar = "ndt"): every arm here runs on
  # the worktree build, where predict() reports a FRACTION even for a
  # model with no ndt_group(). The first version of this script read
  # the global arm's fraction as a time and reported a 740 ms bias.
  hat <- as.numeric(ndt_time(fit, newdata = one))
  cat(sprintf(paste0("nt=%-4d arm=%-7s conv=%d logLik=%12.3f ",
                     "ndt_mean=%.5f bias_ms=%+7.2f sd_hat=%.5f ",
                     "sd_true=%.5f rmse_ms=%6.2f cor=%.4f ",
                     "floor_gap_ms=%6.2f below=%d/%d s=%.1f\n"),
              nt, arm, fit$opt$convergence, as.numeric(logLik(fit)),
              mean(hat), 1000 * (mean(hat) - mean(truth)), sd(hat),
              sd(truth), 1000 * sqrt(mean((hat - truth)^2)),
              suppressWarnings(cor(hat, truth)),
              1000 * mean(own - truth), sum(hat < own), length(hat),
              secs))
  utils::flush.console()
  invisible(NULL)
}

# what the floors alone say about the design, before any fitting
for (nt in c(100L, 200L, 400L)) {
  d <- mk(nt)
  own <- as.numeric(tapply(d$rt, d$s, min))
  truth <- attr(d, "truth")
  cat(sprintf(paste0("nt=%-4d floors: overshoot mean %.2f ms sd %.2f ",
                     "max %.2f | sd(floor) %.5f sd(truth) %.5f | ",
                     "above GLOBAL floor %d/30\n"),
              nt, 1000 * mean(own - truth), 1000 * sd(own - truth),
              1000 * max(own - truth), sd(own), sd(truth),
              sum(truth > min(d$rt))))
}
cat("\n")

arms <- strsplit(Sys.getenv("REV_ARMS", "re,nore,global"), ",")[[1L]]
for (nt in c(100L, 200L, 400L)) {
  for (arm in arms) row(nt, arm)
}

# Lane nlstop2, item 5079: the widest piece between the walls of K + xp
# in median gaps, on the fresh rows' xp (seeds 11 to 55); and a body whose
# walls lie on both sides of its fit, a / (1 + b * x) with x of both
# signs (gaussian, n 300, a 2, b 0.1, seeds 1 to 5): the distance to its
# nearest wall in median gaps, and whether the fit restarts.
# Usage: NLSTOP2_ARM=before|after Rscript dev/nlstop2-gaps.R
source("dev/nlstop2-common.R")
fresh_defs()
ns <- asNamespace("frmtmb")
for (seed in c(11, 22, 33, 44, 55)) {
  d <- gen("gaussian", "mm", seed * 1000 + 1)
  w <- sort(unique(-d$xp))
  g <- diff(w)
  cat(sprintf("xp seed %d: widest piece %.2f median gaps\n", seed,
              max(g) / stats::median(g)))
}
for (seed in 1:5) {
  set.seed(seed)
  d <- data.frame(x = rnorm(300))
  d$y <- 2 / (1 + 0.1 * d$x) + rnorm(300, 0, 0.3)
  r <- tt(frm(bf(y ~ a / (1 + b * x), a + b ~ 1, nl = TRUE), data = d))
  lp <- Filter(function(q) !is.null(q$nl_body), r$frame$linpreds)[[1]]
  wl <- ns$nl_pole_walls(r$frame, lp, r$estimates, "beta", 2L,
                         ns$nl_denoms(lp$nl_body))[[1]]
  cat(sprintf("1 + b x seed %d: b %.4g walls %.3g to %.3g dist %.3g gap %.3g (%.1f gaps) w %s m %s\n",
              seed, r$estimates$beta[2], wl$lo, wl$hi, wl$dist, wl$gap,
              wl$dist / wl$gap, wsum(r), msum(r)))
}
cat("DONE\n")

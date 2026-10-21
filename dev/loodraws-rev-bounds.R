# Reviewer of lane loodraws, fix round 2: the gated same-draws test's
# bounds, evaluated and printed beside the measured difference, on the
# test's own construction (data seed 3, brms chains = 2, iter = 1000,
# seed 1). Also the same bounds against the 0.20.0 r_eff = 1 rule on the
# subset, and the weighted mean's tight bound (2 dl + dx) max|x| without
# the factor S, to show how much headroom S gives.
#   Rscript dev/loodraws-rev-bounds.R -> dev/loodraws-rev-out/fix2-bounds.txt
.libPaths(c("C:/Users/adf44/source/r/wt-loodraws-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
options(mc.cores = 1)
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", dirname(system.file(package = "frmtmb.sample")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(3)
n <- 80
dd <- data.frame(x = stats::rnorm(n), g = factor(rep(1:8, 10)))
dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x + stats::rnorm(8, 0, 0.7)[dd$g])
fitb <- function(f) q(brms::brm(f, data = dd,
  family = brms::brmsfamily("gaussian"), chains = 2, iter = 1000, seed = 1,
  refresh = 0, silent = 2))
mk <- function(bfit, f) {
  fit <- q(frm(f, family = gaussian(), data = dd))
  labs <- frmtmb::brms_par_labels(fit)
  M <- as.matrix(bfit)
  if ("theta_1" %in% labs) M <- cbind(M, theta_1 = M[, "sd_g__Intercept"])
  structure(list(stanfit = bfit$fit,
                 draws = cbind(M[, labs, drop = FALSE], lp__ = 0),
                 fit = fit), class = "frmtmb_draws")
}
bg <- fitb(y ~ x)
br <- fitb(y ~ x + (1 | g))
dg <- mk(bg, bf(y ~ x))
dr <- mk(br, bf(y ~ x + (1 | g)))
eps <- .Machine$double.eps
N <- nrow(dd)
row <- function(lab, d, tol) {
  cat(sprintf("  %-28s diff %-10.3g bound %-10.3g bound/eps %-10.3g\n",
              lab, d, tol, tol / eps))
}
for (k in list(list("gauss", bg, dg), list("ri", br, dr))) {
  cat("==", k[[1]], "\n")
  llf <- unname(unclass(q(log_lik(k[[3]])))[, ])
  llb <- unname(q(brms::log_lik(k[[2]])))
  S <- nrow(llb)
  dl <- max(abs(llf - llb))
  xb <- q(brms::posterior_epred(k[[2]]))
  dx <- max(abs(posterior_epred(k[[3]]) - xb)) / max(abs(xb))
  per <- 2 * dl + dx + eps
  ps <- q(psis(k[[3]]))
  w <- stats::weights(ps, log = FALSE)
  cat(sprintf("  S %d dl %.3g dx %.3g max|x| %.3g\n", S, dl, dx,
              max(abs(xb))))
  row("log_lik guard", dl, S * eps * max(abs(llb)))
  mb <- q(brms::loo_epred(k[[2]]))
  row("loo_epred mean", max(abs(q(loo_epred(k[[3]])) - mb)),
      S * per * max(abs(xb)))
  vb <- q(brms::loo_epred(k[[2]], type = "var"))
  row("loo_epred var", max(abs(q(loo_epred(k[[3]], type = "var")) - vb)),
      2 * S * (2 * dl + 2 * dx + eps) * max(xb^2) /
        (1 - max(colSums(w^2))))
  pr <- c(0.025, 0.5, 0.975)
  L <- max(vapply(seq_len(ncol(xb)), function(i) {
    o <- order(xb[, i]); max(diff(xb[o, i]) / w[o[-1L], i])
  }, 1))
  qb <- q(brms::loo_epred(k[[2]], type = "quantile", probs = pr))
  row("loo_epred quantile",
      max(abs(q(loo_epred(k[[3]], type = "quantile", probs = pr)) - qb)),
      S * (2 * dl + eps) * L + (dx + eps) * max(abs(xb)))
  cat(sprintf("  L (largest gap over weight) %.3g\n", L))
  ids <- seq(1L, 1000L, by = 3L)
  sb <- q(brms::loo_epred(k[[2]], draw_ids = ids))
  sf <- q(loo_epred(k[[3]], draw_ids = ids))
  row("loo_epred(draw_ids)", max(abs(sf - sb)),
      length(ids) * per * max(abs(xb)))
  # the 0.20.0 rule on the same subset
  ll <- log_lik(k[[3]], draw_ids = ids)
  old <- loo::E_loo(posterior_epred(k[[3]], draw_ids = ids),
                    q(loo::psis(-ll, r_eff = 1)))$value
  row("old r_eff = 1 on the subset", max(abs(old - sb[, 1])),
      length(ids) * per * max(abs(xb)))
  cat(sprintf("  tight mean bound without S: %.3g\n",
              (2 * dl + dx + eps) * max(abs(xb))))
}
cat("done\n")

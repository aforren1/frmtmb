# The gated same-draws comparison of test-loo-predict.R, as numbers:
# frmtmb.sample against brms 2.23.0 on brms's own draws, the largest
# difference in units of eps times the largest value. Run under R's
# reference BLAS and under OpenBLAS 0.3.26 with LAPACK, where brms's
# X %*% b and frmtmb's predictor round differently.
#   <Rscript> dev/loodraws-blas-same.R
# Logs: dev/loodraws-log/blas-same-<ref|openblas>.txt
.libPaths(c("C:/Users/adf44/source/r/wt-loodraws-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("R home:", R.home(), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(3)
n <- 80
dd <- data.frame(x = stats::rnorm(n), g = factor(rep(1:8, 10)))
dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x + stats::rnorm(8, 0, 0.7)[dd$g])
fitb <- function(f) q(brms::brm(f, data = dd,
  family = brms::brmsfamily("gaussian"), chains = 2, iter = 1000,
  seed = 1, refresh = 0, silent = 2))
mk <- function(bfit, f) {
  fit <- q(frm(f, family = gaussian(), data = dd))
  labs <- frmtmb::brms_par_labels(fit)
  M <- as.matrix(bfit)
  if ("theta_1" %in% labs) M <- cbind(M, theta_1 = M[, "sd_g__Intercept"])
  structure(list(stanfit = bfit$fit,
                 draws = cbind(M[, labs, drop = FALSE], lp__ = 0),
                 fit = fit), class = "frmtmb_draws")
}
ulps <- function(a, b) {
  max(abs(as.numeric(a) - as.numeric(b))) /
    (.Machine$double.eps * max(abs(as.numeric(b))))
}
bg <- fitb(y ~ x)
br <- fitb(y ~ x + (1 | g))
ks <- list(gauss = list(bg, mk(bg, bf(y ~ x))),
           ri = list(br, mk(br, bf(y ~ x + (1 | g)))))
for (nm in names(ks)) {
  k <- ks[[nm]]
  cat(nm, "log_lik:", ulps(log_lik(k[[2]]), q(brms::log_lik(k[[1]]))), "\n")
  for (ty in c("mean", "var", "quantile")) {
    pr <- if (ty == "quantile") c(0.025, 0.5, 0.975) else 0.5
    cat(nm, "loo_epred", ty, ":",
        ulps(q(loo_epred(k[[2]], type = ty, probs = pr)),
             q(brms::loo_epred(k[[1]], type = ty, probs = pr))), "\n")
  }
  ids <- seq(1L, 1000L, by = 3L)
  cat(nm, "loo_epred draw_ids:", ulps(q(loo_epred(k[[2]], draw_ids = ids)),
      q(brms::loo_epred(k[[1]], draw_ids = ids))), "\n")
  cat(nm, "loo_R2 draws:", ulps(q(loo_R2(k[[2]], seed = 7, summary = FALSE)),
      q(brms::loo_R2(k[[1]], seed = 7, summary = FALSE))), "\n")
}
set.seed(5)
wf <- q(loo_model_weights(ks$gauss[[2]], ks$ri[[2]]))
set.seed(5)
wb <- q(brms::loo_model_weights(bg, br))
cat("stacking:", format(as.numeric(wf), digits = 10), "/",
    format(as.numeric(wb), digits = 10), "ulps", ulps(wf, wb), "\n")

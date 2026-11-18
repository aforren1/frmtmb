# rel084-kt copy of dev/squant-qgam.R: lib = rellib-r23 then user lib (lane lib and rellib-r22 removed). No other change.
# Lane squant, item 4(b): quantile_elf() against qgam 2.0.0 at the same
# smoothing. qgam's elf(qu, co, theta) has the absolute kernel width
# co = lambda sigma and the learning rate sigma = exp(theta); with sigma
# held at exp(theta) and smoothing = co / exp(theta), quantile_elf() is
# the same density. qgam is called with lsig given, so it does not tune
# sigma, and its co comes from its own err rule (.getErrParam()).
#
# Run: Rscript dev/squant-qgam.R > dev/squant-out/qgam.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/rellib-r23",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(qgam)
})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")),
    "qgam", format(packageVersion("qgam")), "mgcv",
    format(packageVersion("mgcv")), "\n")

same_model <- function(qf) {
  co <- unique(qf$family$getCo())
  th <- qf$family$getTheta()
  stopifnot(length(co) == 1L)
  list(sigma = exp(th), smoothing = co / exp(th), co = co, theta = th)
}

# --- 1. parametric models: qgam's gamSim(1) data, linear terms --------
set.seed(2)
dat <- mgcv::gamSim(1, n = 400, dist = "normal", scale = 2, verbose = FALSE)
for (qu in c(0.1, 0.5, 0.8)) {
  for (lsig in c(-2, 0)) {
    qf <- qgam(y ~ x0 + x1 + x2 + x3, data = dat, qu = qu, lsig = lsig)
    m <- same_model(qf)
    ff <- frm(bf(y ~ x0 + x1 + x2 + x3, quantile = qu, sigma = m$sigma),
              family = quantile_elf(smoothing = m$smoothing), data = dat)
    b_q <- stats::coef(qf)
    b_f <- fixef(ff)[, "Estimate"]
    # qgam's log-likelihood: the elf family's aic() is -2 sum log f
    ll_q <- -qf$family$aic(dat$y, stats::fitted(qf), qf$family$getTheta(),
                           rep(1, nrow(dat)), NULL) / 2
    # the score of the coefficients, X' (plogis(r / co) - (1 - qu)) /
    # sigma, at each package's estimate: which one stopped nearer the
    # maximum where the coefficients differ
    X <- stats::model.matrix(~ x0 + x1 + x2 + x3, dat)
    score <- function(b) {
      r <- dat$y - drop(X %*% b)
      max(abs(crossprod(X, stats::plogis(r / m$co) - (1 - qu)))) / m$sigma
    }
    cat(sprintf(paste0("parametric qu %.1f lsig %+.0f: smoothing %.5f (co",
                       " %.5f) | max |coef diff| %.3g, max relative %.3g |",
                       " logLik qgam %.8f frm %.8f diff %.3g | max |score|",
                       " qgam %.2g frm %.2g | conv %d\n"),
                qu, lsig, m$smoothing, m$co, max(abs(b_q - b_f)),
                max(abs(b_q - b_f) / abs(b_q)), ll_q,
                as.numeric(logLik(ff)), ll_q - as.numeric(logLik(ff)),
                score(b_q), score(b_f), ff$opt$convergence))
    # the same fit with nlminb's relative tolerance at 1e-14
    ft <- frm(bf(y ~ x0 + x1 + x2 + x3, quantile = qu, sigma = m$sigma),
              family = quantile_elf(smoothing = m$smoothing), data = dat,
              control = frmtmb_control(optCtrl = list(
                iter.max = 1000, eval.max = 1000, rel.tol = 1e-14)))
    b_t <- fixef(ft)[, "Estimate"]
    cat(sprintf(paste0("  rel.tol 1e-14: max |coef diff| %.3g, max relative",
                       " %.3g, max |score| frm %.2g, conv %d\n"),
                max(abs(b_q - b_t)), max(abs(b_q - b_t) / abs(b_q)),
                score(b_t), ft$opt$convergence))
  }
}

# --- 2. qgam's own examples with smooths ------------------------------
# ?qgam: gamSim(1) with s(x0) + s(x1) + s(x2) + s(x3) at qu = 0.5, and
# MASS::mcycle with s(times, k = 20, bs = "ad") at qu = 0.8. The
# adaptive smooth has no frmtmb equivalent, so mcycle uses bs = "tp".
# qgam (mgcv) selects the smoothing parameters by its REML/LAML; frmtmb
# fits the smooths as random effects by the Laplace approximation with
# REML = TRUE. The sigma is held as above.
cases <- list(
  list(name = "gamSim s(x0..x3) qu 0.5", data = dat, qu = 0.5,
       fq = y ~ s(x0) + s(x1) + s(x2) + s(x3),
       ff = y ~ s(x0) + s(x1) + s(x2) + s(x3)),
  list(name = "mcycle s(times, k = 20) qu 0.8", data = MASS::mcycle,
       qu = 0.8, fq = accel ~ s(times, k = 20),
       ff = accel ~ s(times, k = 20)))
for (cs in cases) {
  qf <- qgam(cs$fq, data = cs$data, qu = cs$qu, lsig = 0)
  m <- same_model(qf)
  fo <- bf(cs$ff, quantile = cs$qu, sigma = m$sigma)
  ws <- character()
  ff <- withCallingHandlers(
    frm(fo, family = quantile_elf(smoothing = m$smoothing), data = cs$data,
        REML = TRUE),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  fq_fit <- as.numeric(stats::fitted(qf))
  ff_fit <- as.numeric(fitted(ff)[, "Estimate"])
  # fitted() is the mean; the location mu is predict's linear predictor
  ff_mu <- as.numeric(frm_linpred(ff, dpar = "mu"))
  cat(sprintf(paste0("%s: smoothing %.5f | max |mu diff| %.3g, relative",
                     " to the range of mu %.3g | edf qgam %.3f | conv %d",
                     " warn %d %s\n"),
              cs$name, m$smoothing, max(abs(fq_fit - ff_mu)),
              max(abs(fq_fit - ff_mu)) / diff(range(fq_fit)),
              sum(qf$edf), ff$opt$convergence, length(ws),
              paste(substr(ws, 1, 60), collapse = " | ")))
}
cat("DONE\n")

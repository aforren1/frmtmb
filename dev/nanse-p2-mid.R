# Punch round 2, RB1: the counted-work rule on mid-size mixed models.
# Per model: outer parameters, the optimizer's counted evaluations, the
# decision (fit time or wait), and the fit-time cost against
# check_se = "ignore" (arms interleaved, each grown past 1 s, minimum of
# 3 rounds; the control re-times the check arm).
#   Rscript dev/nanse-p2-mid.R
.libPaths(c("C:/Users/adf44/source/r/wt-nanse-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
mk <- function(p, n = 1500, G = 40, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("x", 1:p)))
  g <- factor(sample(G, n, TRUE))
  eta <- 0.3 + X %*% rnorm(p, 0, 0.1) + rnorm(G, 0, 0.5)[g]
  data.frame(X, g, y = rpois(n, exp(eta)), yg = eta + rnorm(n))
}
models <- list(
  pois_p12 = list(p = 12, fam = poisson(), y = "y"),
  pois_p20 = list(p = 20, fam = poisson(), y = "y"),
  pois_p40 = list(p = 40, fam = poisson(), y = "y"),
  gaus_p20 = list(p = 20, fam = gaussian(), y = "yg"),
  gaus_p40 = list(p = 40, fam = gaussian(), y = "yg")
)
time_arm <- function(f) {
  k <- 1L
  repeat {
    t0 <- proc.time()[["elapsed"]]
    for (i in seq_len(k)) f()
    el <- proc.time()[["elapsed"]] - t0
    if (el > 1) return(el / k)
    k <- k * 2L
  }
}
for (nm in names(models)) {
  m <- models[[nm]]
  d <- mk(m$p)
  fo <- stats::as.formula(paste(m$y, "~", paste0("x", 1:m$p, collapse = " + "),
                                "+ (1 | g)"))
  fit <- suppressWarnings(frm(fo, data = d, family = m$fam))
  np <- length(fit$opt$par)
  at <- ns$se_check_at_fit(fit)
  chk <- function() suppressWarnings(frm(fo, data = d, family = m$fam))
  ign <- function() suppressWarnings(frm(fo, data = d, family = m$fam,
    control = frmtmb_control(check_se = "ignore")))
  r <- matrix(NA_real_, 3, 3)
  for (i in 1:3) r[i, ] <- c(time_arm(chk), time_arm(ign), time_arm(chk))
  t <- apply(r, 2, min)
  cat(sprintf("%-9s outer pars %3d  evals %4d  2*np/evals %.2f  at fit %-5s | check %.3fs ignore %.3fs ratio %.3f control %.3f\n",
              nm, np, fit$opt$evals, 2 * np / fit$opt$evals, at, t[1], t[2],
              t[1] / t[2], t[3] / t[1]))
}

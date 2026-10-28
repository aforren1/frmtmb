# Models without a cov = TRUE ar(), ma() or arma() term must not move:
# logLik, every estimate and standard error, the resolved priors of a
# MAP fit and frm_sample() draws under a fixed seed, printed to 12
# digits on rellib-r14 ("before") and on the lane library ("after").
# Usage: Rscript dev/sigmainn-unchanged.R before|after
# Log: dev/sigmainn-out/unchanged-<arm>.txt; compare with diff past the
# first line.
arm <- commandArgs(TRUE)[1]
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sigmainn"
.libPaths(c(if (identical(arm, "after")) "C:/Users/adf44/source/r/wt-sigmainn-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
Sys.setenv(FRMTMB_STAN_CACHE = file.path(wt, "dev", "stan-cache"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(11)
n <- 240
d <- data.frame(t = rep(1:6, 40), g = rep(1:40, each = 6),
                x = rnorm(n), z = rnorm(n))
d$y <- 1 + 0.5 * d$x + rnorm(40, 0, 0.5)[d$g] + rnorm(n, 0, exp(0.2 * d$z))
d$cnt <- rpois(n, exp(0.3 + 0.4 * d$x))
d$y2 <- 0.3 * d$y + rnorm(n)
M <- 0.4^abs(outer(1:n, 1:n, "-"))
p12 <- function(lab, v) {
  cat(lab, ":", paste(formatC(unname(v), digits = 12, format = "g"),
                      collapse = " "), "\n")
}
models <- list(
  gauss = list(f = bf(y ~ x), fam = gaussian()),
  re = list(f = bf(y ~ x + (1 | g)), fam = gaussian()),
  sigz = list(f = bf(y ~ x, sigma ~ z), fam = gaussian()),
  student = list(f = bf(y ~ x + (1 | g)), fam = student()),
  pois = list(f = bf(cnt ~ x + (1 | g)), fam = poisson()),
  cosy = list(f = bf(y ~ x + cosy(t, g)), fam = gaussian()),
  cosy_sigz = list(f = bf(y ~ x + cosy(t, g), sigma ~ z), fam = gaussian()),
  unstr = list(f = bf(y ~ x + unstr(t, g)), fam = gaussian()),
  arcond = list(f = bf(y ~ x + ar(t, g)), fam = gaussian()),
  armacond = list(f = bf(y ~ x + arma(t, g)), fam = student()),
  fcor = list(f = bf(y ~ x + fcor(M)), fam = gaussian(), data2 = list(M = M)),
  mv = list(f = mvbf(bf(y ~ x), bf(y2 ~ x), rescor = TRUE), fam = NULL))
for (m in names(models)) {
  mo <- models[[m]]
  fit <- suppressWarnings(if (is.null(mo$fam)) {
    frm(mo$f, data = d)
  } else if (is.null(mo$data2)) {
    frm(mo$f, family = mo$fam, data = d)
  } else frm(mo$f, family = mo$fam, data = d, data2 = mo$data2))
  p12(paste(m, "logLik"), logLik(fit))
  p12(paste(m, "par"), fit$obj$env$last.par.best)
  p12(paste(m, "se"), sqrt(diag(vcov(fit, full = TRUE))))
  if ("sigma" %in% names(fit$family %||% list()) || TRUE) {
    s <- tryCatch(sigma(fit), error = function(e) NA)
    p12(paste(m, "sigma"), s)
  }
}
# MAP fits with priors on sigma, b and the correlation classes
pr <- set_prior("student_t(3, 0, 2.5)", class = "sigma") +
  set_prior("normal(0, 1)", class = "b")
for (f in list(bf(y ~ x), bf(y ~ x + cosy(t, g)), bf(y ~ x + ar(t, g)))) {
  fit <- frm(f, family = gaussian(), data = d, prior = pr)
  p12(paste("MAP", deparse(f$formula), "objective"), fit$opt$objective)
  p12(paste("MAP", deparse(f$formula), "par"), fit$obj$env$last.par.best)
}
fx <- frm(bf(y ~ x, sigma ~ z), family = gaussian(), data = d,
          prior = set_prior("normal(0, 1)", class = "Intercept",
                            dpar = "sigma"))
p12("MAP sigma ~ z objective", fx$opt$objective)
# draws under a fixed seed, on two models
for (f in list(bf(y ~ x), bf(y ~ x + cosy(t, g)))) {
  ds <- suppressWarnings(suppressMessages(frm_sample(
    f, family = gaussian(), data = d, chains = 1, iter = 300, seed = 5,
    refresh = 0)))
  mm <- as.matrix(ds)
  for (j in seq_len(ncol(mm))) {
    p12(paste("draws", deparse(f$formula), colnames(mm)[j]),
        mm[c(1, 50, 150), j])
  }
}

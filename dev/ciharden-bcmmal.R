# test-bcm-latent-mixtures.R "Malingering_2 matches its Stan program"
# fails with OpenBLAS (stated constant 0, measured -1.3e-4 against a
# bound of 5.9e-5). Where does the fit stop, how big are the log-gamma
# terms the two programs cancel there, and how big is the disagreement
# against their rounding scale?
# Usage: Rscript dev/ciharden-bcmmal.R <lib or base>   (FRMTMB_BRMS_FIT_TESTS
# must be true; the Stan program comes from dev/stan-cache)
a <- commandArgs(TRUE)
.libPaths(c(if (!identical(a[1], "base")) a[1],
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
root <- "C:/Users/adf44/source/r/frmtmb-wt-ciharden"
if (!nzchar(Sys.getenv("FRMTMB_STAN_CACHE"))) {
  Sys.setenv(FRMTMB_STAN_CACHE = file.path(root, "dev/stan-cache"))
}
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
cat("R:", R.home(), " threads:", Sys.getenv("OPENBLAS_NUM_THREADS"), "\n")
e <- new.env(parent = asNamespace("frmtmb"))
td <- file.path(root, "tests/testthat")
for (f in list.files(td, "^helper-.*[.]R$", full.names = TRUE)) {
  sys.source(f, envir = e)
}
src <- parse(file.path(td, "test-bcm-latent-mixtures.R"))
for (x in src) {
  if (is.call(x) && identical(x[[1]], as.name("<-")) &&
      grepl("^bcm_", as.character(x[[2]]))) eval(x, e)
}
d <- e$bcm_malingering_data()
fit <- suppressWarnings(
  frm(bf(k | trials(n) ~ 1),
      family = mixture(beta_binomial, beta_binomial), data = d,
      start = e$bcm_mix_start(c(0.55, 0.95), phi = c(15, 15), theta = 0)))
p <- e$bcm_betamix_pars(fit)

res <- e$stan_lp_check(e$bcm_betamix_code(),
                       data = list(p = nrow(d), k = as.integer(d$k), n = 45L),
                       fit = fit, pars = e$bcm_betamix_pars, const = 0,
                       tol = Inf, tol_grad = Inf)
terms <- function(mu, phi) {
  a <- mu * phi; b <- (1 - mu) * phi
  sum(abs(lgamma(d$k + a)) + abs(lgamma(45 - d$k + b)) +
        abs(lgamma(a + b + 45)) + abs(lgamma(a)) + abs(lgamma(b)) +
        abs(lgamma(a + b)))
}
S <- terms(p$mu1, p$phi1) + terms(p$mu2, p$phi2)
cat(sprintf(paste0("phi1 %.6g phi2 %.6g ours %.10g measured %.6g\n",
                   "log-gamma terms S %.6g, eps S %.6g, ",
                   "|measured| / (eps S) %.4g\n"),
            p$phi1, p$phi2, res$ours, res$measured_const, S,
            .Machine$double.eps * S,
            abs(res$measured_const) / (.Machine$double.eps * S)))

# Lane sampleprior: which small negbinomial data set brms 2.23.0 samples
# cleanly under its default inv_gamma(0.4, 0.3) on shape, so that the
# agreement case is a check of the prior and not of a stuck brms chain.
#   Rscript dev/sampleprior-nbprobe.R <lane|r9>
arm <- commandArgs(TRUE)[1]
lib <- if (identical(arm, "r9")) character(0) else
  "C:/Users/adf44/source/r/wt-sampleprior-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sampleprior"
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
source(file.path(wt, "dev/sampleprior-agree-common.R"))
cat("arm", arm, "lib:", find.package("frmtmb.sample"), "\n")
for (n in c(15L, 25L, 40L)) {
  for (size in c(1, 3)) {
    set.seed(20261014)
    d <- data.frame(x = stats::rnorm(n))
    d$y <- stats::rnbinom(n, mu = exp(1 + 0.3 * d$x), size = size)
    cs <- list(form = "y ~ x", fam = "negbinomial", data = d,
               par = "shape")
    b <- suppressWarnings(agree_brms(cs, 3L, 2L, 2000L,
                                     file.path(wt, "dev/stan-cache")))
    f <- agree_frm(cs, 3L, 2L, 2000L)
    lb <- log(b)
    lf <- log(f)
    rh <- posterior::rhat(matrix(lb, ncol = 2))
    z <- (mean(lf) - mean(lb)) /
      sqrt(agree_mcse(lf, 2)^2 + agree_mcse(lb, 2)^2)
    cat(sprintf(paste0("n %2d size %g: brms rhat %.3f mean %.3f q95 %.3f | ",
                       "frmtmb mean %.3f q95 %.3f | z %.2f\n"), n, size, rh,
                mean(lb), stats::quantile(lb, 0.95), mean(lf),
                stats::quantile(lf, 0.95), z))
  }
}

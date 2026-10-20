# Lane sampleprior: the second negbinomial probe. The first
# (dev/sampleprior-nbprobe.R) found brms 2.23.0's chains stuck near
# log(shape) = 50 whenever the counts are close to Poisson, a numerical
# failure of the program rather than a property of the prior; here the
# counts are strongly overdispersed, where inv_gamma(0.4, 0.3) acts on
# small shapes.
#   Rscript dev/sampleprior-nbprobe2.R <lane|r9>
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
for (n in c(10L, 15L)) {
  for (size in c(0.3, 0.5, 0.8)) {
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

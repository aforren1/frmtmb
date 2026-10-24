# Reviewer: does registering a structure named after an R function
# break that function in a formula? Arg 1: library path.
.libPaths(c(strsplit(commandArgs(TRUE)[1], ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
set.seed(1)
d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
d$y <- d$x + rnorm(20)[d$g] + rnorm(100)
r <- function(lab, e) cat(lab, "->", tryCatch({
  f <- e; format(as.numeric(logLik(f)), digits = 10)
}, error = function(e) paste("ERROR:", conditionMessage(e))), "\n")
r("poly(x, 2) unregistered", frm(y ~ poly(x, 2) + (1 | g), data = d))
r("ns via splines", frm(y ~ splines::ns(x, 3) + (1 | g), data = d))
if (exists("frmtmb_register_covstruct")) {
  lp <- function(b, theta, blk) sum(RTMB::dnorm(b, 0, exp(theta[1]), log = TRUE))
  vc <- function(theta, blk) matrix(exp(2 * theta[1]), 1, 1,
                                    dimnames = list(blk$cnms, blk$cnms))
  frmtmb_register_covstruct("poly", lpdf = lp, vcov = vc, npar = 1L,
                            sd_idx = 1L)
  r("poly(x, 2) after registering poly", frm(y ~ poly(x, 2) + (1 | g),
                                             data = d))
  frmtmb_register_covstruct("exp2", lpdf = lp, vcov = vc, npar = 1L,
                            sd_idx = 1L)
  r("I(x^2) control", frm(y ~ I(x^2) + (1 | g), data = d))
  frmtmb_register_covstruct("scale", lpdf = lp, vcov = vc, npar = 1L,
                            sd_idx = 1L)
  r("scale(x) after registering scale", frm(y ~ scale(x) + (1 | g),
                                           data = d))
  frmtmb_register_covstruct("factor", lpdf = lp, vcov = vc, npar = 1L,
                            sd_idx = 1L)
  d$h <- rep(1:4, 25)
  r("factor(h) after registering factor", frm(y ~ factor(h) + (1 | g),
                                             data = d))
}

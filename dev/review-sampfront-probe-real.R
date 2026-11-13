# Review of lane sampfront: the probe/watch hole on real sampler output.
#   Rscript dev/review-sampfront-probe-real.R [base]
# frm_sample(laplace = TRUE) and frm_sample() of y ~ x + (1 + x | g);
# posterior_epred() at newdata x = -1 in group 1. Data seed 77, sampler
# seed 1, 1 chain, 600 iterations.
args <- commandArgs(TRUE)
base <- length(args) && identical(args[1], "base")
.libPaths(c(if (!base) "C:/Users/adf44/source/r/wt-sampfront-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-sampfront/dev/stan-cache")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", if (base) "base" else "lane", find.package("frmtmb.sample"), "\n")
set.seed(77)
G <- 8; n <- 12
dd <- data.frame(g = factor(rep(seq_len(G), each = n)), x = rnorm(G * n))
u0 <- rnorm(G, 0, 0.8); u1 <- rnorm(G, 0, 0.8)
dd$y <- 0.5 + 0.4 * dd$x + u0[dd$g] + u1[dd$g] * dd$x +
  rnorm(nrow(dd), 0, 0.5)
cat("true group-1 value at x = -1:", 0.5 - 0.4 + u0[1] - u1[1], "\n")
smp <- function(...) {
  suppressWarnings(suppressMessages(
    frm_sample(..., chains = 1, iter = 600, refresh = 0, seed = 1)))
}
fit <- suppressWarnings(frm(y ~ x + (1 + x | g), data = dd))
lap <- smp(fit, laplace = TRUE)
ful <- smp(fit)
nd <- data.frame(x = -1, g = factor(1, levels = 1:G))
el <- tryCatch(posterior_epred(lap, newdata = nd), error = function(e) e)
ef <- posterior_epred(ful, newdata = nd)
en <- posterior_epred(ful, newdata = nd, re_formula = NA)
cat("full draws, group 1, x = -1: mean", round(mean(ef), 4), "\n")
cat("full draws, re_formula = NA, x = -1: mean", round(mean(en), 4), "\n")
if (inherits(el, "error")) {
  cat("laplace draws: REFUSED:", substr(conditionMessage(el), 1, 100), "\n")
} else {
  cat("laplace draws, group 1, x = -1: mean", round(mean(el), 4),
      "(returned without a word)\n")
}
cat("DONE\n")

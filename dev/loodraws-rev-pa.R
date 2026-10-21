# Reviewer of lane loodraws: posterior_average(weights = "stacking") on
# draws that carry a stored loo criterion of other draws, base and lane.
#   Rscript dev/loodraws-rev-pa.R <rev|base> -> dev/loodraws-rev-out/pa-<arm>.txt
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(if (a[1] == "rev") "C:/Users/adf44/source/r/wt-loodraws-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", dirname(system.file(package = "frmtmb.sample")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(21)
n <- 60
dd <- data.frame(x = rnorm(n), g = factor(rep(1:6, 10)))
dd$y <- rnorm(n, 1 + 0.5 * dd$x + rnorm(6, 0, 0.6)[dd$g])
mk <- function(fit, nd = 400, sd = 0.05, seed = 5) {
  tpl <- fit$frame$par_template
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  set.seed(seed)
  M <- matrix(rep(est, each = nd) + rnorm(nd * length(est), 0, sd), nd,
              dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  structure(list(stanfit = NULL, draws = M, fit = fit), class = "frmtmb_draws")
}
d1 <- mk(q(frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)))
d2 <- mk(q(frm(bf(y ~ x), family = gaussian(), data = dd)))
w0 <- attr(q(posterior_average(d1, d2, variable = "b_x", seed = 1)), "weights")
s1 <- q(add_criterion(d1, "loo", ndraws = 40))
s2 <- q(add_criterion(d2, "loo", ndraws = 40))
w1 <- attr(q(posterior_average(s1, s2, variable = "b_x", seed = 1)), "weights")
cat("no stored loo:", format(unname(w0), digits = 10), "\n")
cat("stored loo(ndraws = 40):", format(unname(w1), digits = 10), "\n")
cat("identical:", identical(unname(w0), unname(w1)), "\n")

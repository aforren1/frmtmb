# Reviewer probe: wiener(link_bias = "identity") fits that end at a bias
# above 1, against the default logit fit on the same data, and the
# tape value of the guard. Seeds 20261101 to 20261105, n = 80, bias 0.97.
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
for (s in 1:5) {
  set.seed(20261100 + s)
  d2 <- frmtmb.eam:::ddm_simulate(80, mu = 0.5, bs = 1.0, ndt = 0.25,
                                  bias = 0.97)
  q <- function(expr) withCallingHandlers(expr,
    warning = function(w) invokeRestart("muffleWarning"))
  fi <- q(frm(bf(rt | dec(upper) ~ 1), family = wiener(link_bias = "identity"),
              data = d2))
  f0 <- q(frm(bf(rt | dec(upper) ~ 1), family = wiener(), data = d2))
  bi <- summary(fi)$spec_pars["bias", 1]
  b0 <- summary(f0)$spec_pars["bias", 1]
  # the same density, evaluated off the tape at the identity fit's
  # parameters, as the lane's guard intends
  dp <- function(f, nm) frm_linpred(f, dpar = nm, type = "response")
  off <- sum(family(fi)$lpdf(d2$rt, list(mu = dp(fi, "mu"), bs = dp(fi, "bs"),
                                          ndt = dp(fi, "ndt"),
                                          bias = dp(fi, "bias")),
                              list(dec = d2$upper)))
  cat(sprintf("seed %d: identity bias %.5f logLik %.4f (off-tape density sum %s) | logit bias %.5f logLik %.4f\n",
              s, bi, as.numeric(logLik(fi)), format(off), b0,
              as.numeric(logLik(f0))))
}

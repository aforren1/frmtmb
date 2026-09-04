# Probe 6: the lapse overlay.
#
# The paper's Eq. 14 (M1) and Eq. 13 (M5) both carry a uniform contaminant:
#   p*_i(t) = 0.95 p_i(t - t_nd) + 0.025
# so a fit of "the paper's 3-parameter DDM" without a lapse component is not the paper's
# model. It also matters practically: PyDDM's tutorial puts non-decision time at 0.211 s
# when the fastest observed RT is 0.203 s, which is only admissible because the
# contaminant carries the trials the diffusion cannot reach. wiener()'s ndt link ceiling
# is min(y), so without a lapse component frmtmb cannot go there by construction.
#
# Question: does core's mixture() reach this today? The survey says mixture() accepts any
# component with a mu dpar, no extra_pars and no drop_intercept, and that a dpar can be
# pinned to a constant by passing a named scalar to bf(). wiener() satisfies all of those.

sp <- Sys.getenv("EL_SCRATCH")
.libPaths(c(file.path(sp, "el-lib"), .libPaths()))
suppressPackageStartupMessages({
  library(frmtmb); library(frmtmb.ddm); library(RWiener)
})

T_DUR <- 2
COH <- c(0, 0.032, 0.064, 0.128, 0.256, 0.512)
LAPSE <- 0.05

set.seed(42)
truth <- list(mu0 = 10.5, bs = 1.6, ndt = 0.28)
dat <- do.call(rbind, lapply(COH, function(cc) {
  x <- RWiener::rwiener(440, alpha = truth$bs, tau = truth$ndt, beta = 0.5,
                        delta = truth$mu0 * cc)
  data.frame(rt = x$q, upper = as.integer(x$resp == "upper"), coh = cc)
}))
lap <- runif(nrow(dat)) < LAPSE
dat$rt[lap] <- runif(sum(lap), 0.05, T_DUR)          # fast guesses AND slow lapses
dat$upper[lap] <- rbinom(sum(lap), 1, 0.5)
cat(nrow(dat), "trials,", sum(lap), "contaminated; min RT", round(min(dat$rt), 4), "\n")
cat("(the fastest trial is now below the true ndt of", truth$ndt, ")\n\n")

# A uniform contaminant over [0, T_dur], split evenly across the two responses. It needs
# a mu dpar to satisfy mixture(); the density ignores it.
unif_contam <- custom_family(
  "unif_contam", dpars = "mu", links = list(mu = "identity"),
  lpdf = function(y, dpars, aterms) {
    0 * dpars$mu + log(0.5 / 2)
  }
)

cat("=== attempt 1: mixture(wiener(), unif_contam), lapse weight free ===\n")
a1 <- try(frm(bf(rt | vint(upper) ~ 0 + coh, bias1 = 0.5) +
                mixture(wiener(max_ndt = 0.4), unif_contam),
              data = dat), silent = TRUE)
if (inherits(a1, "try-error")) cat("REFUSED:\n", attr(a1, "condition")$message, "\n")

cat("\n=== attempt 2: same, lapse weight pinned at 0.05 (the paper's value) ===\n")
a2 <- try(frm(bf(rt | vint(upper) ~ 0 + coh, bias1 = 0.5, theta1 = qlogis(1 - LAPSE)) +
                mixture(wiener(max_ndt = 0.4), unif_contam),
              data = dat), silent = TRUE)
if (inherits(a2, "try-error")) {
  cat("REFUSED:\n", attr(a2, "condition")$message, "\n")
} else {
  print(summary(a2))
}

cat("\n=== attempt 3: no mixture, for contrast (ndt capped at min RT) ===\n")
a3 <- try(frm(bf(rt | vint(upper) ~ 0 + coh, bias = 0.5), family = wiener(), data = dat),
          silent = TRUE)
if (inherits(a3, "try-error")) {
  cat("REFUSED:\n", attr(a3, "condition")$message, "\n")
} else {
  fe <- fixef(a3)
  cat("mu0 =", round(fe$mu[1], 4), " (true", truth$mu0, ")\n")
  cat("bs  =", round(exp(fe$bs[1]), 4), " (true", truth$bs, ")\n")
  cat("ndt =", round(plogis(fe$ndt[1]) * min(dat$rt), 4), " (true", truth$ndt,
      "; ceiling", round(min(dat$rt), 4), ")\n")
  cat("\nThis is the cost of having no lapse component: a handful of contaminant trials\n",
      "drags the ndt ceiling down and the whole fit with it.\n", sep = "")
}

# Does the wiener lpdf return NaN or -Inf when rt < ndt? That decides whether a mixture
# can rescue those rows at all, since NaN would poison the log-sum-exp.
cat("\n=== behavior of the wiener density below ndt ===\n")
v <- frmtmb.ddm:::ddm_lpdf_both(c(-0.05, -1e-8, 1e-8, 0.05), rep(2, 4), rep(1.5, 4),
                                rep(0.5, 4), rep(1L, 4))
print(data.frame(t_minus_ndt = c(-0.05, -1e-8, 1e-8, 0.05), lpdf = as.numeric(v)))

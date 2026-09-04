# Probe 7: closing the lapse gap with a custom family.
#
# Probe 6 showed two things block the paper's overlay:
#   (a) wiener() refuses max_ndt above min(rt), so mixture() never gets built, and
#   (b) the density returns NaN below the non-decision time, which would poison a
#       log-sum-exp even if the guard were lifted.
# Both refusals are correct for a bare Wiener model and wrong for a mixture, because the
# contaminant is exactly what covers the rows the diffusion cannot reach.
#
# This probe asks whether a user can route around both today, without any change to core
# or to frmtmb.ddm, by writing one custom_family that carries the overlay itself. If yes,
# the paper's M1 and M2 go from BLOCKED to FITS WITH WORKAROUND and the cost is countable.

sp <- Sys.getenv("EL_SCRATCH")
.libPaths(c(file.path(sp, "el-lib"), .libPaths()))
suppressPackageStartupMessages({
  library(frmtmb); library(frmtmb.ddm); library(RWiener)
})

T_DUR <- 2
COH <- c(0, 0.032, 0.064, 0.128, 0.256, 0.512)
LAPSE <- 0.05

# The whole workaround. Note ndt gets a log link, not wiener()'s scaled_logit, so it is
# free to exceed the fastest observed RT. The zero-density region is handled without a
# branch: relu(y - ndt) + eps pins the argument to eps below the non-decision time, where
# the Wiener log-density is about -2.8e7, so exp() of it underflows to a clean zero. The
# lapse term keeps the sum strictly positive, so no row can ever return -Inf or NaN.
wiener_lapse <- function(lapse = 0.05, t_dur = 2) {
  custom_family(
    "wiener_lapse",
    dpars = c("mu", "bs", "ndt", "bias"),
    links = list(mu = "identity", bs = "log", ndt = "log", bias = "logit"),
    lpdf = function(y, dpars, aterms) {
      d <- y - dpars$ndt
      safe <- 0.5 * (d + abs(d)) + 1e-8
      lw <- frmtmb.ddm:::ddm_lpdf_both(safe, dpars$mu, dpars$bs, dpars$bias, aterms$vint1)
      log((1 - lapse) * exp(lw) + lapse * 0.5 / t_dur)
    },
    # init_dpars entries are function(y, aterms), NOT scalars. Passing scalars is accepted
    # silently by custom_family() and by check_custom_family(), then fails inside frm()
    # with "could not find function init_fn".
    init_dpars = list(mu = function(y, aterms) 0,
                      bs = function(y, aterms) 1.5,
                      ndt = function(y, aterms) 0.5 * min(y),
                      bias = function(y, aterms) 0.5)
  )
}

set.seed(42)
truth <- list(mu0 = 10.5, bs = 1.6, ndt = 0.28)
dat <- do.call(rbind, lapply(COH, function(cc) {
  x <- RWiener::rwiener(440, alpha = truth$bs, tau = truth$ndt, beta = 0.5,
                        delta = truth$mu0 * cc)
  data.frame(rt = x$q, upper = as.integer(x$resp == "upper"), coh = cc)
}))
lap <- runif(nrow(dat)) < LAPSE
dat$rt[lap] <- runif(sum(lap), 0.05, T_DUR)
dat$upper[lap] <- rbinom(sum(lap), 1, 0.5)
dat$cohf <- factor(dat$coh, levels = COH)
cat(nrow(dat), "trials,", sum(lap), "contaminated; min RT", round(min(dat$rt), 4),
    "vs true ndt", truth$ndt, "\n\n")

fam <- wiener_lapse(LAPSE, T_DUR)

cat("=== check_custom_family ===\n")
chk <- try(check_custom_family(fam, y = dat$rt,
                               dpars = list(mu = rep(2, nrow(dat)), bs = rep(1.5, nrow(dat)),
                                            ndt = rep(0.28, nrow(dat)),
                                            bias = rep(0.5, nrow(dat))),
                               aterms = list(vint1 = dat$upper)), silent = TRUE)
if (inherits(chk, "try-error")) cat("FAILED:", attr(chk, "condition")$message, "\n") else
  cat("passed\n")

cat("\n=== M1 with the paper's lapse overlay ===\n")
f1 <- try(frm(bf(rt | vint(upper) ~ 0 + coh, bias = 0.5) + fam,
              data = dat), silent = TRUE)
if (inherits(f1, "try-error")) {
  cat("REFUSED:\n", attr(f1, "condition")$message, "\n")
} else {
  fe <- fixef(f1)
  est <- c(mu0 = unname(fe$mu[1]), bs = unname(exp(fe$bs[1])), ndt = unname(exp(fe$ndt[1])))
  cat("logLik", round(logLik(f1), 2), "\n")
  print(data.frame(truth = unlist(truth), estimate = round(est, 4)))
  cat("\ncompare probe 6's no-lapse fit on these same data:",
      "mu0 5.37, bs 2.08, ndt 0.045\n")
}

cat("\n=== M2 (free drift per coherence) with the lapse overlay ===\n")
f2 <- try(frm(bf(rt | vint(upper) ~ 0 + cohf, bias = 0.5) + fam,
              data = dat), silent = TRUE)
if (inherits(f2, "try-error")) {
  cat("REFUSED:\n", attr(f2, "condition")$message, "\n")
} else {
  cat("drift per coherence vs truth:\n")
  print(round(data.frame(coh = COH, est = unname(fixef(f2)$mu),
                         truth = truth$mu0 * COH), 3))
}

# ---- and on the real monkey data, where the ndt ceiling actually bound -------------
csv <- file.path(sp, "el-pyddm", "roitman_rts.csv")
if (file.exists(csv)) {
  cat("\n=== real Roitman monkey 1, M1 with lapse ===\n")
  rd <- read.csv(csv)
  rd <- subset(rd, monkey == 1 & rt > 0.1 & rt < 1.65)
  rd$upper <- as.integer(rd$correct == 1)
  f3 <- try(frm(bf(rt | vint(upper) ~ 0 + coh, bias = 0.5) + fam, data = rd),
            silent = TRUE)
  if (inherits(f3, "try-error")) {
    cat("REFUSED:\n", attr(f3, "condition")$message, "\n")
  } else {
    fe <- fixef(f3)
    cat("mu0 =", round(fe$mu[1], 4), "  bs =", round(exp(fe$bs[1]), 4),
        "  ndt =", round(exp(fe$ndt[1]), 4), "\n")
    cat("min RT =", min(rd$rt), "-- ndt is now free to sit above it if the data ask.\n")
    cat("PyDDM tutorial GDDM for this subject had nondectime = 0.211.\n")
  }
}

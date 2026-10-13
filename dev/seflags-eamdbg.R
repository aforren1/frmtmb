# Lane seflags, punch round 2: frmtmb.eam's test-defects.R:92 (a
# wiener plus lognormal mixture whose Wiener parameters are flat):
# which parameters are lost and why, the null basis, and which of
# residuals(type = "response")'s rows get NaN, on a given library.
#   Rscript dev/seflags-eamdbg.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
ns <- asNamespace("frmtmb")
cat("lib", find.package("frmtmb"), "\n")
set.seed(3)
dat <- frmtmb.eam:::ddm_simulate(300, mu = 0.9, bs = 1.4, ndt = 0.30)
k <- sample(300, 18)
dat$rt[k] <- stats::runif(18, 0.12, 2.5)
fit <- suppressWarnings(frm(bf(rt | dec(upper) ~ 1, bias1 = 0.5),
  family = mixture(wiener(max_ndt = 0.4, allow_unreachable = TRUE),
                   lognormal()), data = dat))
sdr <- suppressWarnings(ns$sdr_of(fit))
print(sdr$se_lost)
print(round(fit$opt$par, 3))
fl <- sdr$se_null_flat
cat("flat cols", length(fl$cols), "tol", fl$tol, "ill", !is.null(fl$ill), "\n")
cat("dscale", signif(fl$dscale, 3), "\nunits", signif(fit$par_units %||% 1, 3), "\n")
print(signif(sdr$se_null, 3))
w <- character()
r <- withCallingHandlers(stats::residuals(fit, type = "response"),
  warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
cat("residual Est.Error NaN:", sum(is.nan(r[, "Est.Error"])), "of", nrow(r), "| warnings", length(w), "\n")
print(summary(r[, "Est.Error"]))
# the cosines jc_nonest() reads, in both metrics, for those rows
trace("jc_nonest", where = ns, print = FALSE, tracer = quote({
  N0 <- jc$null; G0 <- as.matrix(G); n0 <- NROW(G0)
  cosm <- function(s) {
    Nq <- N0 * s; Gq <- G0 / rep(s[pos], each = n0)
    max(abs(Gq %*% Nq[pos, , drop = FALSE]) /
          outer(sqrt(rowSums(Gq^2)), sqrt(colSums(Nq^2))))
  }
  s1 <- jc$flat$dscale; s1[!is.finite(s1) | !(s1 > 0)] <- 1
  u <- jc$units %||% rep(1, NROW(N0))
  cat(sprintf("largest cosine: unit-diagonal %.3g | optimizer units %.3g\n",
              cosm(s1), cosm(1 / u)))
}))
invisible(suppressWarnings(stats::residuals(fit, type = "response")))

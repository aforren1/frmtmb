# Release 0.70.0: the "poly" part of dev/optrap-rev-att/scripts/
# setier-rev-probe.R (same data) with (1 | g), against an exact
# reference: the GLS covariance (X' V^-1 X)^-1 from a Householder QR of
# the whitened design, at frmtmb's own variance estimates and at lmer's.
# The probe compares with lmer()'s vcov(), which carries lmer's own
# estimates and arithmetic, so a gap there does not say which side is
# off.
#   Rscript dev/rel070-optrap-diag-gls.R <lib> <label>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
cat("lib", dirname(find.package("frmtmb")), "frmtmb",
    as.character(packageVersion("frmtmb")), "label", a[2], "\n")
gls_se <- function(X, Z, sig, tau) {
  V <- sig^2 * diag(nrow(X)) + tau^2 * tcrossprod(Z)
  Xw <- backsolve(chol(V), X, transpose = TRUE)
  q <- qr(Xw, LAPACK = TRUE)
  Ri <- backsolve(qr.R(q), diag(ncol(X)))
  se <- sqrt(rowSums(Ri^2))
  out <- numeric(ncol(X))
  out[q$pivot] <- se
  out
}
for (scale in c(1, 1e5)) for (deg in 5:8) {
  set.seed(deg)
  d <- data.frame(x = runif(200, 1, 2) * scale, g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x / scale) + rnorm(20, 0, 0.3)[d$g] +
    rnorm(200, 0, 0.2)
  for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
  rhs <- paste0("x", 1:deg, collapse = " + ")
  fo <- as.formula(paste("y ~", rhs, "+ (1 | g)"))
  w <- character()
  f <- withCallingHandlers(suppressMessages(frm(fo, data = d)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    })
  se <- suppressWarnings(fixef(f))[, "Est.Error"]
  X <- cbind(1, as.matrix(d[paste0("x", 1:deg)]))
  Z <- model.matrix(~ 0 + g, d)
  vc <- as.data.frame(VarCorr(f))
  sig_f <- exp(f$estimates$betad[[1]])
  tau_f <- exp(f$estimates$theta[[1]])
  ref_f <- gls_se(X, Z, sig_f, tau_f)
  m <- tryCatch(suppressMessages(suppressWarnings(
    lmer(fo, data = d, REML = FALSE))), error = function(e) e)
  if (inherits(m, "error")) {
    ref_l <- rep(NA, ncol(X)); se_l <- ref_l; dll <- NA
    sig_l <- NA; tau_l <- NA
  } else {
    vl <- as.data.frame(VarCorr(m))
    sig_l <- sigma(m); tau_l <- vl$sdcor[vl$grp == "g"]
    ref_l <- gls_se(X, Z, sig_l, tau_l)
    se_l <- sqrt(diag(as.matrix(vcov(m))))
    dll <- as.numeric(logLik(f)) - as.numeric(logLik(m))
  }
  r <- function(x, y) {
    v <- abs(x / y - 1)
    if (all(!is.finite(v))) NA else max(v[is.finite(v)])
  }
  cat(sprintf(paste0("scale %g deg %d code %d | sd f %.6g/%.6g lmer ",
                     "%.6g/%.6g | logLik f - lmer %.3g | finite SE %d of %d",
                     " | max |SE/QR(own)-1| %.3g | |SE/QR(lmer's sds)-1| ",
                     "%.3g | lmer vcov/QR(lmer's sds) %.3g | warn %s\n"),
              scale, deg, f$opt$convergence, sig_f, tau_f, sig_l, tau_l,
              dll, sum(is.finite(se)), length(se), r(se, ref_f),
              r(se, ref_l), r(se_l, ref_l),
              if (length(w)) substr(paste(w, collapse = " || "), 1, 90)
              else "none"))
}

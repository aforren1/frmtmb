# Reviewer of lane illse, attack 9: frmtmb.spline's curves on root fits.
#   Rscript dev/illse-rev-spline.R <old|new>
# old: the lane core with the 0.70.0 frmtmb.spline (rellib-r8); new: the
# lane core with the lane's frmtmb.spline. Raw polynomials of degree 5
# (not lossy) and 8 (lossy) through I(x^k); every curve route against
# lm()'s QR.
ext <- commandArgs(TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb, lib.loc = "C:/Users/adf44/source/r/wt-illse-lib")
  library(frmtmb.spline, lib.loc = if (ext == "old") {
    "C:/Users/adf44/source/r/rellib-r8"
  } else "C:/Users/adf44/source/r/wt-illse-lib")
})
cat("frmtmb", find.package("frmtmb"), "| spline", find.package("frmtmb.spline"),
    format(packageVersion("frmtmb.spline")), "\n")
cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(tryCatch(expr, error = function(e)
    structure(conditionMessage(e), class = "rev_err")),
    warning = function(x) { w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning") })
  list(v = v, w = w)
}
for (deg in c(5, 8)) for (seed in 1:3) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(200, 1, 2))
  d$y <- sin(3 * d$x) + stats::rnorm(200, 0, 0.3)
  fo <- stats::as.formula(paste("y ~ x +",
                                paste0("I(x^", 2:deg, ")", collapse = " + ")))
  fit <- suppressWarnings(frm(fo, data = d))
  lmf <- stats::lm(fo, data = d)
  n <- nrow(d); k <- deg + 1
  Ri <- backsolve(qr.R(lmf$qr), diag(k))
  s <- sqrt(sum(stats::residuals(lmf)^2) / n)
  g <- data.frame(x = seq(1.05, 1.95, length.out = 7))
  X1 <- outer(g$x, 0:deg, `^`)
  ref_c <- sqrt(rowSums((X1 %*% Ri)^2)) * s
  Dd <- outer(g$x, 0:deg, function(x, j) j * x^pmax(j - 1, 0))
  ref_d <- sqrt(rowSums((Dd %*% Ri)^2)) * s
  g2 <- data.frame(x = g$x + 0.05)
  X2 <- outer(g2$x, 0:deg, `^`)
  ref_k <- sqrt(rowSums(((X1 - X2) %*% Ri)^2)) * s
  one <- function(lab, r, ref) {
    if (inherits(r$v, "rev_err")) return(sprintf("%s ERROR(%s)", lab, substr(r$v, 1, 70)))
    se <- r$v$.se
    sprintf("%s %.2e%s", lab, max(abs(se / ref - 1)),
            if (length(r$w)) paste0(" [warn: ", substr(r$w[1], 1, 50), "]") else "")
  }
  cv <- cap(frm_curve(fit, newdata = g, simultaneous = FALSE))
  ck <- cap(frm_curve(fit, newdata = g, contrast = g2, simultaneous = FALSE))
  cd <- cap(frm_curve_deriv(fit, var = "x", order = 1, newdata = g,
                            simultaneous = FALSE))
  cf <- cap(frm_curve_feature(fit, var = "x", type = "minimum", newdata = data.frame(x = seq(1, 2, length.out = 41))))
  cfs <- if (inherits(cf$v, "rev_err")) sprintf("feature ERROR(%s)", substr(cf$v, 1, 70)) else
    paste("feature ok se", paste(signif(unlist(as.data.frame(cf$v)[1, grep("se", names(as.data.frame(cf$v)))]), 4), collapse = " "))
  cat(sprintf("deg %d s%d | root %s | %s | %s | %s | %s\n", deg, seed,
              !is.null(frmtmb:::sdr_of(fit)$cov_fixed_root),
              one("curve", cv, ref_c), one("contrast", ck, ref_k),
              one("deriv", cd, ref_d), cfs))
}

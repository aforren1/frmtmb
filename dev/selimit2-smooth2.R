# Lane selimit2, backlog 2757: `s(x) + s(x2)` with x2 = x was named
# under OpenBLAS LAPACK and not under the reference BLAS
# (dev/setier-rev-smooth2.R, on 0.69.0). Rerun on 0.82.0 and the lane;
# run it with the R of dev/ciharden-out/ for the OpenBLAS arm.
#   Rscript dev/selimit2-smooth2.R <base|lane>
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-selimit2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("arm", arm, dirname(find.package("frmtmb")), "| BLAS",
    basename(extSoftVersion()[["BLAS"]]), La_library(), "\n")
for (s in 1:5) {
  set.seed(s)
  d <- data.frame(x = runif(300))
  d$x2 <- d$x
  d$y <- sin(2 * pi * d$x) + rnorm(300, 0, 0.3)
  w <- character(); m <- character()
  f <- tryCatch(withCallingHandlers(frm(y ~ s(x) + s(x2), data = d),
    warning = function(z) {w <<- c(w, conditionMessage(z))
      invokeRestart("muffleWarning")},
    message = function(z) {m <<- c(m, conditionMessage(z))
      invokeRestart("muffleMessage")}), error = function(e) e)
  if (inherits(f, "error")) {
    cat("seed", s, "ERROR", conditionMessage(f), "\n")
    next
  }
  nm <- ns$outer_par_names(f)
  th <- grep("^theta", nm)
  p <- f$opt$par
  f0 <- f$obj$fn(p)
  # the trade-off: raise one sd and lower the other
  q <- p; q[th[1]] <- q[th[1]] + 0.5; q[th[2]] <- q[th[2]] - 0.5
  lost <- ns$sdr_of(f)$se_lost
  se <- suppressWarnings(sqrt(diag(vcov(f, full = TRUE))))
  cat(sprintf(paste0("seed %d code %d | theta %s | SE %s | lost %s | ",
                     "nll trade +-0.5 %.3g | warn %d msg %d %s\n"),
              s, f$opt$convergence, paste(signif(p[th], 4), collapse = ","),
              paste(signif(se[th], 3), collapse = ","),
              if (length(lost)) paste(names(lost), lost, sep = ":",
                                      collapse = ",") else "none",
              f$obj$fn(q) - f0, length(w), length(m),
              if (length(m)) substr(m[1], 1, 60) else ""))
}

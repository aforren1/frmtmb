# Reviewer of lane setier: a smoothing sd that is flat because it trades
# off with another (two smooths of the same covariate), not because it is
# at its limit. Does se_boundary_names() call it a boundary (silent)?
#   Rscript dev/setier-rev-smooth2.R lane|base
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-setier-lib",
           "C:/Users/adf44/source/r/rellib-r6"),
  base = "C:/Users/adf44/source/r/rellib-r6")
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("arm", arm, "\n")
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
  if (inherits(f, "error")) {cat("seed", s, "ERROR", conditionMessage(f), "\n"); next}
  nm <- ns$outer_par_names(f)
  th <- grep("^theta", nm)
  p <- f$opt$par
  f0 <- f$obj$fn(p)
  # the trade-off: raise one sd and lower the other
  q <- p; q[th[1]] <- q[th[1]] + 0.5; q[th[2]] <- q[th[2]] - 0.5
  ed <- vapply(th, function(j) {q2 <- p; q2[j] <- q2[j] - 2
    f$obj$fn(q2) - f0}, 0)
  lost <- ns$sdr_of(f)$se_lost
  se <- suppressWarnings(sqrt(diag(vcov(f, full = TRUE))))
  cat(sprintf("seed %d code %d | theta %s | SE %s | lost %s | nll: trade +-0.5 %.3g, each sd -2: %s | warn %d msg %d\n",
              s, f$opt$convergence, paste(signif(p[th], 4), collapse = ","),
              paste(signif(se[th], 3), collapse = ","),
              if (length(lost)) paste(names(lost), lost, sep = ":",
                                      collapse = ",") else "none",
              f$obj$fn(q) - f0, paste(signif(ed, 3), collapse = ","),
              length(w), length(m)))
}

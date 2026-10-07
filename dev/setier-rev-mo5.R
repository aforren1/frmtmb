# Reviewer of lane setier: the mo() seeds whose verdict the lane changed
# (3, 15, 103, 142, 171 of brms_monotonic's data code, interaction).
# Are zeta2_1 and zeta2_2 flat at the estimate? Exact Hessian spectrum,
# the likelihood along each coordinate, and base's SEs.
#   Rscript dev/setier-rev-mo5.R lane|base
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-setier-lib",
           "C:/Users/adf44/source/r/rellib-r6"),
  base = "C:/Users/adf44/source/r/rellib-r6")
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, find.package("frmtmb"), "\n")
mk <- function(s) {
  set.seed(s)
  lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
  ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
  d <- data.frame(income, ls)
  d$age <- rnorm(100, mean = 40, sd = 10)
  d
}
for (s in c(3, 15, 103, 142, 171)) {
  f <- suppressWarnings(frm(ls ~ mo(income) * age, data = mk(s)))
  nm <- frmtmb:::outer_par_names(f)
  p <- f$opt$par
  H <- f$obj$he(p)
  D <- sqrt(abs(diag(H)))
  ev <- eigen(H / outer(D, D), symmetric = TRUE)
  se <- suppressWarnings(sqrt(diag(vcov(f, full = TRUE))))
  cat(sprintf("seed %d code %d logLik %.10f\n", s, f$opt$convergence,
              as.numeric(logLik(f))))
  cat("  par:", paste0(nm, "=", signif(p, 4), collapse = " "), "\n")
  cat("  SE :", paste0(nm, "=", signif(se, 4), collapse = " "), "\n")
  cat("  unit-diag ev / max:", signif(ev$values / max(abs(ev$values)), 3),
      "\n")
  f0 <- f$obj$fn(p)
  for (z in c("zeta2_1", "zeta2_2")) {
    j <- match(z, nm)
    dl <- vapply(c(-5, -1, 1, 5), function(h) {
      q <- p; q[j] <- q[j] + h; f$obj$fn(q) - f0
    }, 0)
    cat("  ", z, "diag H", signif(H[j, j], 3), "nll change at -5,-1,+1,+5:",
        signif(dl, 3), "\n")
  }
  w <- exp(c(0, f$estimates$zeta2)); w <- w / sum(w)
  cat("  simplex 2:", signif(w, 4), "\n")
}

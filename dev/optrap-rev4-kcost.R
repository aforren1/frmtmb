# Reviewer of lane optrap, consolidation check: the time of
# predictor_design_kappa() on designs at and above the 5e6-cell cap.
.libPaths(c("C:/Users/adf44/source/r/wt-optrap-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
kap <- frmtmb:::predictor_design_kappa
for (cfg in list(c(1e6, 5), c(1e5, 50), c(1e6, 50), c(2e5, 100))) {
  n <- cfg[1]; p <- cfg[2]
  set.seed(1)
  Z <- matrix(rnorm(n * (p - 1)), n, p - 1)
  colnames(Z) <- paste0("z", seq_len(p - 1))
  d <- data.frame(Z, g = factor(rep(seq_len(1000), length.out = n)))
  d$y <- rnorm(n)
  fo <- as.formula(paste("y ~", paste(colnames(Z), collapse = " + "), "+ (1 | g)"))
  fr <- frm(fo, data = d, dry_run = "frame")
  t <- system.time(k <- kap(fr))[["elapsed"]]
  cat(sprintf("n %g p %d | kappa %.3g | %.2f s\n", n, p, k, t))
}

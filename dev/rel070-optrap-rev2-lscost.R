# Reviewer of lane optrap, re-check: the unit predictor_units() gives,
# and the time it takes, on high-R^2 designs on either side of
# predictor_resid_sd()'s cap of 5e6 design cells (n * p), from the frame
# alone (dry_run = "frame").
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
pu <- frmtmb:::predictor_units
for (cfg in list(c(9.9e4, 50), c(1.01e5, 50), c(9.9e5, 5), c(1.01e6, 5),
                 c(1e6, 4))) {
  n <- cfg[1]; p <- cfg[2]
  set.seed(1)
  Z <- matrix(rnorm(n * (p - 1)), n, p - 1)
  colnames(Z) <- paste0("z", seq_len(p - 1))
  d <- data.frame(Z, g = factor(rep(seq_len(1000), length.out = n)))
  d$y <- 1 + 1e4 * d$z1 + rnorm(1000, 0, 0.5)[d$g] + rnorm(n)
  fo <- as.formula(paste("y ~", paste(colnames(Z), collapse = " + "),
                         "+ (1 | g)"))
  fr <- frm(fo, data = d, dry_run = "frame")
  t <- system.time(u <- pu(fr))[["elapsed"]]
  cat(sprintf("n %g p %d cells %.3g | unit %s | predictor_units() %.2f s\n",
              n, p, n * p, paste(u, collapse = ","), t))
}
